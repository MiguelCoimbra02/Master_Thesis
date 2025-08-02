from flask import Flask, render_template, request, jsonify, session, send_file, after_this_request
from pathlib import Path
from flask_cors import CORS, cross_origin
from pdfReport import PDFReport
import time
import psutil
import io
import csv


#Using Pickle to serialize the NetworkX graph to save it in the session
import pickle
#Using Redis to store the session data
import redis

# --- Neo4j Service and utils ---
from neo4j_graph_service import Neo4jGraphService

from utils import *

# --- PDF generation ---
from io import BytesIO
import base64
from PIL import Image
import os

#Initiate Neo4j object to interact with the db
neo4_service = Neo4jGraphService()
base_dir = Path('../')

app = Flask(__name__)
CORS(app)

# Using the default Redis configuration
redis_host = os.getenv('REDIS_HOST', 'localhost')
redis_port = int(os.getenv('REDIS_PORT', 6380))

r = redis.StrictRedis(host=redis_host, port=redis_port, db=0)


@app.route('/')
def main():
    return render_template('index.html')


# Search route to handle POST requests from the JS
@app.route('/search', methods=['POST'])
def search():
    print('Received POST request on /search')
    try:
        # Retrieve the JSON data from the request
        data = request.get_json(force=False)
        
        if not data or 'nodes' not in data:
            return {'error': 'No nodes or data provided.'}, 400
        
        query_nodes = {node.upper() for node in data['nodes']} # requested_genes
        tf_ranks = list(data['ranks'])
        slideRange_co_exp = data['rangeSliderValue']
        rangeSliderValue_viz_nodes = data['rangeSliderValue_viz_nodes']
        
        #----- Check database for valid nodes -----
        not_valid_nodes, valid_nodes = neo4_service.search_gene_in_neo4j(query_nodes)
        
        #----- Extract subgraph data only for valid nodes -----
        nodes_data, edges_data, max_result = neo4_service.get_subgraph_data(valid_nodes)
        print('EDGES DATA HERE - ',edges_data)
    
        #----- Create a NetworkX graph from the extracted neo4j data -----
        subgraph_flask = create_graph(nodes_data, edges_data)
      
        #----- Filter edges and nodes based on the query -----
        filtered_sub_graph = filter_edges_nodes(subgraph_flask, valid_nodes, tf_ranks, slideRange_co_exp)
        graph_4_display = get_top_edges_by_irpscore(filtered_sub_graph, rangeSliderValue_viz_nodes) #aqui - filtered
        # Delete the existing graph_key if it exists
        r.delete('graph_key')
        #save the subgraph in the flask session
        r.set('graph_key', pickle.dumps(filtered_sub_graph.copy())) #aqui - filtered
        r.set('graph_key_exports', pickle.dumps(filtered_sub_graph.copy()))

        graph_json = graph2json(graph_4_display, valid_nodes) 
        
        return jsonify({"found": valid_nodes, "not_found": not_valid_nodes, "max_result": max_result, "graph": graph_json})
    except Exception as e:
        return {'error': f'Invalid query data: {str(e)}'}, 400
    
@app.route('/report', methods=['GET', 'POST'])
def generate_report():
    try:
        # Retrieve the JSON data from the request
        data = request.get_json(force=False)
        if not data or 'quried_nodes' not in data or 'tf_ranks' not in data or 'rangeSliderValue' not in data:
            return {'error': 'Not all data required provided.'}, 400
        quried_nodes = set(data['quried_nodes']) # validNodes
        tf_ranks = list(data['tf_ranks'])

        subgraph_flask = pickle.loads(r.get('graph_key_exports'))

        if 4 in tf_ranks:
            tf_ranks.remove(4)
        rangeSliderValue = data['rangeSliderValue']
        nodes_dict = dict(subgraph_flask.nodes(data=True))  
        edges_dict = subgraph_flask.edges(keys=True, data=True)
        network_image_base64 = data['network_image']


               # Handle image if provided
        image_path = None
        if network_image_base64:
            try:
                # Decode the Base64 string to bytes
                network_image_data = base64.b64decode(network_image_base64.split(",")[1])
                # Convert to PIL Image for use in the PDF
                image = Image.open(BytesIO(network_image_data))
                # Save the image temporarily
                image_path = '/tmp/network_image.png'
                image.save(image_path)
            except Exception as e:
                return jsonify({'error': f'Error decoding or saving image: {str(e)}'}), 400

        #generate pdf
        pdf = PDFReport(quried_nodes)
        pdf.add_page()

        pdf_file_path = 'gene_relation_report.pdf'
        
        # Add content to PDF
        pdf.add_summary(rangeSliderValue, tf_ranks)
        pdf.add_legend()
        pdf.add_edges_table(edges_dict)
        pdf.add_nodes_table(nodes_dict)
        pdf.addImage(image_path)
        pdf.output(pdf_file_path)
        byte_string = pdf.output(dest='S')  # Write PDF content to the BytesIO stream
        stream = BytesIO(byte_string)
        stream.seek(0)  # Move cursor to the start of the stream

        response = send_file(stream, as_attachment=True, download_name='gene_relation_report.pdf', mimetype='application/pdf')

        # Delete the file after sending
        @after_this_request
        def cleanup(response):
            try:
                os.remove(pdf_file_path)
            except Exception as e:
                print(f'Failed to delete {pdf_file_path}. Reason: {e}')
            return response
        
        return response
            
    except Exception as e:
        return {'error': f'Invalid data: {str(e)}'}, 500    


@app.route('/export_nodes')
def export_nodes():
    try:
        subgraph_flask = pickle.loads(r.get('graph_key_exports'))
        nodes = subgraph_flask.nodes(data=True)


        if not nodes:
            return "No nodes to export.", 400

        output = io.StringIO()
        writer = csv.writer(output)

        header = ['name', 'arabidopsis_gene', 'isTR', 'isTF', 'gene_annotation']
        writer.writerow(header)

        for node, data in nodes:
            writer.writerow([
                data.get('name', ''),
                data.get('arabidopsis_gene', ''),
                data.get('isTR', ''),
                data.get('isTF', ''),
                data.get('gene_annotation', '')
            ])

        output.seek(0)
        return send_file(
            io.BytesIO(output.getvalue().encode('utf-8')),
            mimetype='text/csv',
            as_attachment=True,
            download_name='nodes.csv'
        )
    except Exception as e:
        app.logger.error("CSV node export error: %s", str(e))
        return f"Error generating node CSV: {str(e)}", 500
    

@app.route('/export_edges')
def export_edges():
    try:
        
        subgraph_flask = pickle.loads(r.get('graph_key_exports'))
        edges = subgraph_flask.edges(keys=True, data=True)
        print('edges_export', edges)
        print(len(edges))

        if not edges:
            return "No edges to export.", 400

        output = io.StringIO()
        writer = csv.writer(output)

        header = ['label', 'source', 'target', 'ConnecTF', 'irp_score', 'edge_type', 'gene_name', 'cis_elements', 'dap_seq', 'tf_rank']
        writer.writerow(header)
        for u, v, k, data in edges:
            if data.get('interaction') == 'interacts with':
                print('ENTREI AQUI')
                writer.writerow([
                    data.get('id', '--'),
                    '--',
                    '--',
                    '--',
                    data.get('irp_score', '--'),
                    data.get('interaction', '--'),
                    '--',
                    '--',
                    '--',
                    '--'
                ])
            
            else:
                print('ENTREI AQUI 2')
                writer.writerow([
                    data.get('id', '--'),
                    data.get('source', '--'),
                    data.get('target', '--'),
                    data.get('connecTF', '--'),
                    data.get('irp_score', '--'),
                    data.get('interaction', '--'),
                    data.get('gene_name', '--'),
                    data.get('cis_elements', '--'),
                    data.get('dap_seq', '--'),
                    data.get('tf_rank', '--')
                ])

        output.seek(0)
        return send_file(
            io.BytesIO(output.getvalue().encode('utf-8')),
            mimetype='text/csv',
            as_attachment=True,
            download_name='edges.csv'
        )
    except Exception as e:
        return f"Error generating CSV: {str(e)}", 500
    

@app.route('/expand', methods=['GET', 'POST'])
def expand():
    try:
        data = request.get_json(force=False)
        query_node = data['node_name']
        tf_ranks = list(data['ranks'])
        slideRange_co_exp = data['rangeSliderValue']
    except Exception as e:
        return {'error': f'Invalid query data: {str(e)}'}, 400
    
    #loading the subgraph from redis
    subgraph_flask = pickle.loads(r.get('graph_key'))


    subgraph_expanded, new_edges, max_result = neo4_service.expand_node(subgraph_flask, query_node)

    subgraph_completed = autocomplete_graph_data(subgraph_expanded, new_edges)  # Autocomplete edges for the new ones, the others already have this information
    print(subgraph_completed.edges(data=True))
    filtered_sub_graph = filter_edges_nodes(subgraph_completed, [query_node], tf_ranks, slideRange_co_exp)
    
    graph4display = expanded_graph4display(subgraph_flask,filtered_sub_graph)
    
    #saving graph to redis
    r.set('graph_key_exports', pickle.dumps(filtered_sub_graph))
    # write potential edges in JSON
    json_data = graph2json(graph4display)
    # print('graph_json: ', json_data)

    return jsonify({"graph": json_data, "max_result": max_result })


#------------ TEMPLATES ------------

@app.route('/about')
def about():
    return render_template('about_template.html')

@app.route('/user_guide')
def user_guide():
    return render_template('user_guide_template.html')

#------------ Full data exports ------------

DOWNLOAD_FOLDER = os.path.join(base_dir, 'viz_tool', 'data')
app.config['DOWNLOAD_FOLDER'] = DOWNLOAD_FOLDER

@app.route('/downloads/<filename>', methods=['GET'])
def download_file(filename):
    try:
        # Ensure the filename is safe
        safe_filename = os.path.basename(filename)
        file_path = os.path.join(app.config['DOWNLOAD_FOLDER'], safe_filename)
        print('file_path', file_path)
        
        if not os.path.exists(file_path):
            return jsonify({"error": "File not found"}), 404
        
        return send_file(file_path, as_attachment=True)
    except Exception as e:
        return jsonify({"error": str(e)}), 500




#------------ Testing - app.run(debug=True) ------------

@app.before_request
def start_timer():
    request.start_time = time.time()


@app.after_request
def log_response_time(response):
    if hasattr(request, "start_time"):
        elapsed_time = time.time() - request.start_time
        app.logger.info(f"Endpoint {request.path} took {elapsed_time:.4f} sec")
    return response

@app.after_request
def log_memory_usage(response):
    process = psutil.Process()
    memory_used = process.memory_info().rss / (1024 * 1024)  # Convert to MB
    app.logger.info(f"Memory usage: {memory_used:.2f} MB")
    return response


# Run the app
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5001,debug=True)
