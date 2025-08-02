#from utils import autocomplete_graph_data #remove later

# Start the Neo4j server before running this script
# sudo systemctl start neo4j
# sudo systemctl stop neo4j

# Check the server status and option to access the Neo4j browser
# sudo systemctl status neo4j

# if prefered to work directly on cypher shell:
# cypher-shell -u neo4j -p password (miguel@itqb)
# if want do delete all data: run MATCH (n) DETACH DELETE n; in cypher shell
import os
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
from flask import jsonify
import csv

# class Neo4jGraphService:
#     def __init__(self, username = 'neo4j', password = 'miguel@itqb', uri = 'bolt://neo4j:7687'):
#         self.driver = GraphDatabase.driver(uri, auth=(username, password))
#         self.subgraph = None

#     #--- Connection Utils ---

#     def get_session(self):
#         return self.driver.session()

#     def close(self):
#         self.driver.close()^

class Neo4jGraphService:
    def __init__(self):
        username = os.getenv('NEO4J_USER', 'neo4j')
        password = os.getenv('NEO4J_PASSWORD', 'miguel@itqb')
        uri = os.getenv('NEO4J_URI', 'bolt://neo4j:7687')
        
        self.driver = None
        
        # Retry logic: wait until Neo4j is ready
        for attempt in range(10):
            try:
                self.driver = GraphDatabase.driver(uri, auth=(username, password))
                # Test connection
                with self.driver.session() as session:
                    session.run("RETURN 1")
                print("Connected to Neo4j successfully.")
                break
            except ServiceUnavailable:
                print(f"Neo4j not ready, retrying {attempt+1}/10...")
                import time
                time.sleep(3)
        else:
            raise Exception("Could not connect to Neo4j after 10 retries.")
        
        self.subgraph = None
    
    
#-------------- loading utils ----------------

# Function to load CSV into Neo4j
    def load_nodes_to_neo4j(self, data_dir, csv_file_name):
        '''
        Loads nodes from a CSV file into a Neo4j database, one by one using parameterized Cypher queries.
        '''
        with self.get_session() as session:
            print("Connected to Neo4j!")

            # Open the CSV file for reading
            with open(f'{data_dir}/{csv_file_name}', 'r') as file:
                reader = csv.DictReader(file)

                for row in reader:
                    cypher_query = """
                    MERGE (n:GeneNode {name: $name})
                    SET n.arabidopsis_gene = $arabidopsis_gene,
                        n.isTF = $isTF,
                        n.isTR = $isTR,
                        n.gene_annotation = $gene_annotation
                    """

                    parameters = {
                        "name": row["name"],
                        "arabidopsis_gene": row["arabidopsis_gene"],
                        "isTF": row["isTF"],
                        "isTR": row["isTR"],
                        "gene_annotation": row["gene_annotation"]
                    }

                    # ✅ Moved inside the loop
                    session.run(cypher_query, parameters)

            print(f"Nodes from {csv_file_name} have been successfully loaded into Neo4j!")


    # Function to load CSV into Neo4j
    def load_edges_to_neo4j(self, data_dir, csv_file_name):
        '''
        Loads edges from a CSV file into a Neo4j database.
        Creates an interaction from node A to node B, with all the associated attributes.
        '''

        with self.get_session() as session:
            print("Connected to Neo4j!")

            # Open the CSV file for reading
            with open(f'{data_dir}/{csv_file_name}', 'r') as file:
                reader = csv.DictReader(file) # Read the CSV as a dictionary
                batch = []

                for row in reader:
                    # Create a Cypher query using parameters for values
                    if int(row["directed"]) == 0:
                        cypher_query = """
                        MATCH (Source:GeneNode {name: $source_name}), (Target:GeneNode {name: $target_name})
                        MERGE (Source)-[e:Interaction {id: $id}]-(Target)
                        SET e.irp_score = $irp_score,
                            e.interaction = $interaction,
                            e.connecTF = $connecTF,
                            e.edge_type = $edge_type,
                            e.gene_name = $gene_name,
                            e.cis_elements = $cis_elements,
                            e.cis_value = $cis_value,
                            e.dap_seq = $dap_seq,
                            e.tf_rank = $tf_rank,
                            e.directed = $directed
                        """
                    else:
                        cypher_query = """
                        MATCH (Source:GeneNode {name: $source_name}), (Target:GeneNode {name: $target_name})
                        MERGE (Source)-[e:Interaction {id: $id}]-(Target)
                        SET e.irp_score = $irp_score,
                            e.interaction = $interaction,
                            e.connecTF = $connecTF,
                            e.edge_type = $edge_type,
                            e.gene_name = $gene_name,
                            e.cis_elements = $cis_elements,
                            e.cis_value = $cis_value,
                            e.dap_seq = $dap_seq,
                            e.tf_rank = $tf_rank,
                            e.directed = $directed,
                            e.source = $source_name,
                            e.target = $target_name
                        """
                    # Create the parameters dictionary
                    parameters = {
                        "source_name": row["Source"],
                        "target_name": row["Target"],
                        "id": row["id"],
                        "irp_score": row["irp_score"],
                        "interaction": row["interaction"],
                        "connecTF": row["connecTF"],
                        "edge_type": row["EDGE_TYPE"],
                        "gene_name": row["gene_name"],
                        "cis_elements": row["cis_elements"],
                        "cis_value": row["cis_value"],
                        "dap_seq": row["dap_seq"],
                        "tf_rank": row["tf_rank"],
                        "directed": row["directed"]
                    }

                    # Add the query with parameters to the batch
                    batch.append((cypher_query, parameters))

                    # When the batch reaches a certain size, run the queries in a single transaction
                    if len(batch) >= 1000:
                        for query, params in batch:
                            session.run(query, params)
                        print(f"Inserted {len(batch)} edges.")
                        batch = []  # Reset batch for next round

                # Run any remaining queries in the batch
                if batch:
                    for query, params in batch:
                        session.run(query, params)
                    print(f"Inserted remaining {len(batch)} edges.")

            print(f"Edges from {csv_file_name} have been successfully loaded into Neo4j!")


    
    def search_gene_in_neo4j(self, query_nodes):
        if not query_nodes:
            return jsonify({'error': 'No query nodes provided'}), 400

        query_nodes = list(query_nodes)  # Ensure it's a list

        with self.get_session() as session:
            # Query to find the nodes that exist in the database
            cypher_query = """
            MATCH (n:GeneNode)
            WHERE n.name IN $query_nodes
            RETURN n.name AS name
            """
            result = session.run(cypher_query, query_nodes=query_nodes)  # Execute query

            found_nodes = {record["name"] for record in result}  # Extract found nodes
            not_found_nodes = list(set(query_nodes) - found_nodes)  # Find missing nodes

        valid_nodes = list(found_nodes)  # Convert set to list

        print(f"Valid nodes: {valid_nodes}")
        print(f"Nodes not found: {not_found_nodes}")

        return not_found_nodes , valid_nodes

    def remove_edges_from_neo4j(self, data_dir, csv_file_name):
        '''
        Removes edges from the Neo4j database based on a CSV file.
        The CSV must contain an 'id' column that uniquely identifies relationships.
        '''

        with self.get_session() as session:
            print("Connected to Neo4j!")

            # Open the CSV file for reading
            with open(f'{data_dir}/{csv_file_name}', 'r') as file:
                reader = csv.DictReader(file)
                batch = []

                for row in reader:
                    edge_id = row.get("id")
                    if edge_id:
                        cypher_query = """
                        MATCH ()-[e:Interaction {id: $id}]->()
                        DELETE e
                        """
                        parameters = {"id": edge_id}
                        batch.append((cypher_query, parameters))

                    # Execute in batches
                    if len(batch) >= 1000:
                        for query, params in batch:
                            session.run(query, params)
                        print(f"Deleted {len(batch)} edges.")
                        batch = []

                # Delete remaining edges
                if batch:
                    for query, params in batch:
                        session.run(query, params)
                    print(f"Deleted remaining {len(batch)} edges.")

            print(f"Edges from {csv_file_name} have been successfully removed from Neo4j.")

    def get_subgraph_data(self, valid_nodes):
        max = False
        with self.get_session() as session:
            gene_nodes = {}  # Dictionary to store node attributes
            edges = set()  # Set to store edges and avoid duplicates

            if len(valid_nodes) == 1:
                # Case 1: Single gene → Get direct neighbors
                cypher_query = """
                MATCH (n:GeneNode)-[r]-(m:GeneNode)
                WHERE (n.name = $query_node OR m.name = $query_node)
                RETURN n, m, r;
                """
                # Execute the query with the first valid node"""
                result = session.run(cypher_query, query_node=valid_nodes[0])           

                for record in result:
                    n, m, r = record["n"], record["m"], record["r"]

                    # Store node attributes
                    gene_nodes[n["name"]] = dict(n)
                    gene_nodes[m["name"]] = dict(m)

                    # Store edge attributes (converted to a tuple to ensure hashability)
                    edges.add((n["name"], m["name"], tuple(sorted(r.items()))))  

            else:
                # Case 2: Multiple genes → Get direct neighbors and shortest paths
                for gene in valid_nodes:
                    cypher_query = """
                    MATCH (n:GeneNode)-[r]-(m:GeneNode)
                    WHERE (n.name = $query_node or m.name = $query_node)
                    RETURN n, m, r;
                    """
                    result = session.run(cypher_query, query_node=gene)

                    for record in result:
                        n, m, r = record["n"], record["m"], record["r"]
                        
                        # Store node attributes
                        gene_nodes[n["name"]] = dict(n)
                        gene_nodes[m["name"]] = dict(m)

                        # Store edge attributes (converted to a tuple to ensure hashability)
                        edges.add((n["name"], m["name"], frozenset(r.items())))  # Use frozenset for edge attributes

                # Find shortest paths between all pairs of valid nodes
                for i in range(len(valid_nodes)):
                    for j in range(i + 1, len(valid_nodes)):
                        cypher_query = """
                        MATCH (n1:GeneNode), (n2:GeneNode),
                            p = shortestPath((n1)-[*]-(n2))
                        WHERE n1.name = $node1 AND n2.name = $node2
                        RETURN nodes(p) AS path_nodes, relationships(p) AS path_edges; 
                        """
                        result = session.run(cypher_query, node1=valid_nodes[i], node2=valid_nodes[j])
                                
                        for record in result:
                            path_nodes = record["path_nodes"]
                            path_edges = record["path_edges"]

                            # Store node attributes
                            for node in path_nodes:
                                gene_nodes[node["name"]] = dict(node)

                            # Store edge attributes (converted to a tuple to ensure hashability)
                            for rel in path_edges:
                                source = rel.start_node["name"]
                                target = rel.end_node["name"]
                                edges.add((source, target, frozenset(rel.items())))  # Use frozenset for edge attributes
        return gene_nodes, edges, max

    
    def expand_node(self, graph, node_name):
        g = graph.copy()
        max = False

        """Expands a node by retrieving its direct neighbors and additional interactions with existing nodes."""
        
        existing_nodes = set(g.nodes())  # Nodes already in the subgraph
        new_edges = []  # Store new edges
        new_nodes = set()  # Track new nodes added
        
        query = """
            MATCH (n:GeneNode)-[r:Interaction]-(m:GeneNode)
            WHERE (n.name = $node_name OR m.name = $node_name)
            RETURN n, m, r;
            """
        
        with self.get_session() as session:
            results = session.run(query, node_name=node_name)            

            for record in results:
                node_data = record["n"]
                neighbor_data = record["m"]
                edge_data = record["r"]

                source = node_data["name"]
                target = neighbor_data["name"]
                edge_attrs = dict(edge_data)

                # Add new node with all properties
                if source not in g:
                    g.add_node(source, **dict(node_data))
                if target not in g:
                    g.add_node(target, **dict(neighbor_data))
                    new_nodes.add(target)

                # Know if this edge is directed or not
                
                if not g.has_edge(source, target):
                    g.add_edge(source, target, **edge_attrs)
                    new_edges.append((source, target, edge_attrs))
                

        # Second step: Get interactions between new nodes and existing nodes
        if new_nodes:
            query = """
            MATCH (a)-[r]-(b)
            WHERE ((a.name IN $new_nodes AND b.name IN $existing_nodes) 
                OR (b.name IN $new_nodes AND a.name IN $existing_nodes))
            RETURN a, b, r;
            """
            with self.get_session() as session:
                results = session.run(query, new_nodes=list(new_nodes), existing_nodes=list(existing_nodes))
            

                for record in results:
                    node_data = record["a"]
                    neighbor_data = record["b"]
                    edge_data = record["r"]

                    source = node_data["name"]
                    target = neighbor_data["name"]
                    edge_attrs = dict(edge_data)

                    # Ensure nodes are added with attributes
                    if source not in g:
                        g.add_node(source, **dict(node_data))
                    if target not in g:
                        g.add_node(target, **dict(neighbor_data))

                    if not g.has_edge(source, target):
                        g.add_edge(source, target, **edge_attrs)
                        new_edges.append((source, target, edge_attrs))
              
        return g, new_edges, max 