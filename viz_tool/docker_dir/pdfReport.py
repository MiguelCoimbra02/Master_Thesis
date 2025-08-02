from fpdf import *
from datetime import datetime

#cell:
# width - height - text - border - linebreak - alignement
#   0       10      ...      0/1      0/1         L C R

from fpdf import FPDF, XPos, YPos
from datetime import datetime

class PDFReport(FPDF):

    def __init__(self, queried_genes):
        super().__init__()
        self.queried_genes = queried_genes

    def header(self):
        queried_genes_strings = [str(gene) for gene in self.queried_genes]
        self.set_font('Helvetica', 'B', 14)  # Use Helvetica
        self.cell(0, 10, f'Gene Interaction Report - {", ".join(queried_genes_strings)}', new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')
        self.cell(0, 10, f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')
        self.ln(10)

    def footer(self):
        self.set_y(-15)
        self.set_font('Helvetica', 'I', 8)  # Use Helvetica
        self.cell(0, 10, f'Page {self.page_no()}', new_x=XPos.RIGHT, new_y=YPos.TOP)

    def add_summary(self, irp_score, tf_rank):
        if len(tf_rank) == 0:
            tf_rank = [0, 1, 2, 3]            
        self.set_font('Helvetica', 'B', 12)  # Use Helvetica
        self.cell(0, 6, 'Analysis Summary', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font('Helvetica', '', 10)  # Use Helvetica
        self.cell(0, 6, f'- Minimum IRP Score: {irp_score}', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        tf_rank_strings = [str(rank) for rank in tf_rank]
        self.cell(0, 6, f'- Trancription factor Ranks searched: {", ".join(tf_rank_strings)}', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(10)

    def add_legend(self):
        self.set_font('Helvetica', 'B', 10)  # Use Helvetica
        self.cell(0, 4, 'TF Ranks Legend:', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font('Helvetica', '', 8)  # Use Helvetica
        legend_items = [
            '- 0 : Co-Expression Relations predicted using Seidr Software',
            '- 1 : Regulation predicted in ConnecTF software',
            '- 2 : Presence of cis-elements in the promoter region of the target gene',
            '- 3 : Data from the DAP-seq experiments (Only for LOC111997151, LOC112008346 and LOC112030452 )'
        ] 
        for item in legend_items:
            self.cell(0, 4, item, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(5)
    

    def add_edges_table(self, edges_dict):
        col_widths = [60, 13, 12, 15, 50, 15, 10]  
        total_width = sum(col_widths)  # Total width of the table
        # Calculate X position to center the table
        page_width = 210  
        x_offset = (page_width - total_width) / 2  # Offset to center the table

        # Move to the calculated X position
        self.set_x(x_offset)
        self.set_font('Helvetica', 'B', 12)  # Use Helvetica
        self.cell(0, 6, 'Edge Table', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font('Helvetica', 'B', 6)
        # Updated Column headers with 'id'
        headers = ['ID','Direction', 'IRP Score', 'TF Rank', 'ConnecTF Method', 'Product Gene', 'cis Value']
        

        # Header row
        self.set_x(x_offset)
        for header, width in zip(headers, col_widths):
            self.cell(width, 8, header, border=1, align='C')
        self.ln()

        # Set font for the table content
        self.set_font('Helvetica', '', 6)

        # Loop through the dictionary and add rows to the table
        for source, target, key, edge_data in edges_dict:
            if edge_data['directed'] in ('0', 0):
                irp_score = edge_data['irp_score']
                directed = 'undirected'
                cis_value = '--'
                edge_type = '--'
                gene_name = '--'
            else:
                irp_score = '--'
                directed = 'directed'
                edge_type = edge_data['edge_type']
                gene_name = edge_data['gene_name']
                cis_value = edge_data['cis_value']
            edge_id = edge_data['id']  # The key itself is the edge ID
            tf_rank = str(edge_data['tf_rank'])


            
            
            self.set_x(x_offset)

            # Add a row for each edge, starting with the id
            self.cell(col_widths[0], 7, edge_id, border=1)  
            self.cell(col_widths[1], 7, directed, border=1)
            self.cell(col_widths[2], 7, irp_score, border=1)
            self.cell(col_widths[3], 7, tf_rank, border=1)
            self.cell(col_widths[4], 7, edge_type, border=1)
            self.cell(col_widths[5], 7, gene_name, border=1)
            self.cell(col_widths[6], 7, cis_value, border=1)
            self.ln()

        self.ln(8)

    def add_nodes_table(self, nodes_dict):
        col_widths = [30, 40, 8, 8, 89]  # Adjust the widths as needed
        total_width = sum(col_widths)  # Total width of the table

        # Calculate X position to center the table
        page_width = 210  
        x_offset = (page_width - total_width) / 2  # Offset to center the table


        # Move to the calculated X position
        self.set_x(x_offset)
        self.set_font('Helvetica', 'B', 12)  # Use Helvetica
        self.cell(0, 6, 'Node Table', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font('Helvetica', 'B', 8)
        
        # Column headers
        headers = ['ID', 'Arabidopsis Homolog', 'TR', 'TF', 'Annotations']
        

        # Header row
        self.set_x(x_offset)
        for header, width in zip(headers, col_widths):
            self.cell(width, 8, header, border=1, align='C')
        self.ln()

        # Set font for the table content
        self.set_font('Helvetica', '', 8)

        # Loop through the dictionary and add rows to the table
        for key, node_data in nodes_dict.items():
            node_id = str(key)  # The key itself is the node ID
            homolog = str(node_data.get('arabidopsis_gene'))  # Get homolog or default to '--'
            is_tr = str(node_data.get('isTR'))  # Get isTR or default to '--'
            is_tf = str(node_data.get('isTF'))  # Get isTF or default to '--'
            annotations = node_data.get('gene_annotation', '--')  # Get annotations or default to '--'
            
            is_tr = 'no' if is_tr == '0' else 'yes'
            is_tf = 'no' if is_tf == '0' else 'yes'
            homolog = '--' if homolog == 'None' else homolog
            
            self.set_x(x_offset)

            # Add a row for each node
            self.cell(col_widths[0], 7, node_id, border=1)  
            self.cell(col_widths[1], 7, homolog, border=1)
            self.cell(col_widths[2], 7, is_tr, border=1)
            self.cell(col_widths[3], 7, is_tf, border=1)
            self.cell(col_widths[4], 7, annotations, border=1)
            self.ln()

        self.add_page() # Add some space after the table

    def addImage(self, path, h=150):
        self.set_font('Helvetica', 'B', 12)  # Use Helvetica
        self.cell(0, 2, 'Network Image', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        # Add your image, scale it down to fit
        self.image(path, h=h, x=0, alt_text='Network screenshot')
        self.ln(2)
        self.set_font('Helvetica', '', 8)  # Use Helvetica
        self.cell(0, 2, 'For better image quality procede to download the Network via export Network in the tools EXPORT option ', new_x=XPos.LMARGIN, new_y=YPos.NEXT)


#generate PDF
'''node_dict = {'LOC112007793': {'Arabidopsis_concise': 'AT2G32300', 'Arabidopsis_gene': 'AT2G32300.1', 'AverageShortestPathLength': 7.496770025839793, 'BetweennessCentrality': 0.0017220678142381, 'ClosenessCentrality': 0.1333907798362774, 'Degree': 4, 'Eccentricity': 18, 'Flow': 0.000145994, 'IsSingleNode': False, 'NeighborhoodConnectivity': 6.75, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 4, 'Radiality': 0.9532606472961166, 'SelfLoops': 0, 'Stress': 5886036, 'TopologicalCoefficient': 0.2613636363636363, 'id': 'LOC112007793', 'isDEG': None, 'isTF': None, 'isTR': None, 'label': 'LOC112007793', 'title': {}}, 'LOC112015362': {'Arabidopsis_concise': None, 'Arabidopsis_gene': None, 'AverageShortestPathLength': 7.18249354005168, 'BetweennessCentrality': 0.010042615019298, 'ClosenessCentrality': 0.139227413769843, 'Degree': 3, 'Eccentricity': 18, 'Flow': 0.000103692, 'IsSingleNode': False, 'NeighborhoodConnectivity': 11, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 3, 'Radiality': 0.9555216292082614, 'SelfLoops': 0, 'Stress': 53733712, 'TopologicalCoefficient': 0.3333333333333333, 'id': 'LOC112015362', 'isDEG': None, 'isTF': None, 'isTR': None, 'label': 'LOC112015362', 'title': {}}, 'LOC111998585': {'Arabidopsis_concise': 'AT4G35160', 'Arabidopsis_gene': 'AT4G35160.1', 'AverageShortestPathLength': 6.55765503875969, 'BetweennessCentrality': 0.0143566830629939, 'ClosenessCentrality': 0.1524935352789065, 'Degree': 17, 'Eccentricity': 19, 'Flow': 0.000641682, 'IsSingleNode': False, 'NeighborhoodConnectivity': 3.176470588235294, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 17, 'Radiality': 0.9600168702247504, 'SelfLoops': 0, 'Stress': 76887364, 'TopologicalCoefficient': 0.0806100217864923, 'id': 'LOC111998585', 'isDEG': 'DEG', 'isTF': None, 'isTR': None, 'label': 'LOC111998585', 'title': {}}, 'LOC112009257': {'Arabidopsis_concise': 'AT2G38970', 'Arabidopsis_gene': 'AT2G38970.1', 'AverageShortestPathLength': 8.394056847545219, 'BetweennessCentrality': 0, 'ClosenessCentrality': 0.1191319070340157, 'Degree': 1, 'Eccentricity': 20, 'Flow': 3.39888e-05, 'IsSingleNode': False, 'NeighborhoodConnectivity': 5, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 1, 'Radiality': 0.946805346420538, 'SelfLoops': 0, 'Stress': 0, 'TopologicalCoefficient': 0, 'borderWidth': 2, 'color': {'background': '#f09d62', 'border': '#e3530b'}, 'id': 'LOC112009257', 'isDEG': None, 'isTF': None, 'isTR': None, 'label': 'LOC112009257', 'title': {}}, 'LOC111984382': {'Arabidopsis_concise': None, 'Arabidopsis_gene': None, 'AverageShortestPathLength': 8.064114987080103, 'BetweennessCentrality': 0.0013472972787204, 'ClosenessCentrality': 0.1240061682654757, 'Degree': 5, 'Eccentricity': 19, 'Flow': 0.00019686, 'IsSingleNode': False, 'NeighborhoodConnectivity': 3.8, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 5, 'Radiality': 0.9491790288699272, 'SelfLoops': 0, 'Stress': 4379248, 'TopologicalCoefficient': 0.2153846153846154, 'id': 'LOC111984382', 'isDEG': 'DEG', 'isTF': None, 'isTR': None, 'label': 'LOC111984382', 'title': {}}, 'LOC112015374': {'Arabidopsis_concise': 'AT1G03870', 'Arabidopsis_gene': 'AT1G03870.1', 'AverageShortestPathLength': 6.887596899224806, 'BetweennessCentrality': 0.0014857313016808, 'ClosenessCentrality': 0.1451885199774901, 'Degree': 4, 'Eccentricity': 18, 'Flow': 0.000169592, 'IsSingleNode': False, 'NeighborhoodConnectivity': 13.75, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 4, 'Radiality': 0.9576431877753612, 'SelfLoops': 0, 'Stress': 7981286, 'TopologicalCoefficient': 0.2833333333333333, 'id': 'LOC112015374', 'isDEG': None, 'isTF': None, 'isTR': None, 'label': 'LOC112015374', 'title': {}}, 'LOC111997687': {'Arabidopsis_concise': 'AT5G62280', 'Arabidopsis_gene': 'AT5G62280.1', 'AverageShortestPathLength': 6.41311369509044, 'BetweennessCentrality': 0.0059972226516229, 'ClosenessCentrality': 0.155930496096701, 'Degree': 6, 'Eccentricity': 18, 'Flow': 0.00021304, 'IsSingleNode': False, 'NeighborhoodConnectivity': 13.333333333333334, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 6, 'Radiality': 0.9610567360065436, 'SelfLoops': 0, 'Stress': 59711796, 'TopologicalCoefficient': 0.1813725490196078, 'id': 'LOC111997687', 'isDEG': None, 'isTF': None, 'isTR': None, 'label': 'LOC111997687', 'title': {}}, 'LOC112015478': {'Arabidopsis_concise': 'AT3G12460', 'Arabidopsis_gene': 'AT3G12460.1', 'AverageShortestPathLength': 7.39421834625323, 'BetweennessCentrality': 0.0062121842199872, 'ClosenessCentrality': 0.135240799388446, 'Degree': 5, 'Eccentricity': 19, 'Flow': 0.000178488, 'IsSingleNode': False, 'NeighborhoodConnectivity': 2, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 5, 'Radiality': 0.9539984291636456, 'SelfLoops': 0, 'Stress': 20717248, 'TopologicalCoefficient': 0.2, 'id': 'LOC112015478', 'isDEG': 'DEG', 'isTF': None, 'isTR': None, 'label': 'LOC112015478', 'title': {}}, 'LOC111987461': {'Arabidopsis_concise': 'AT2G32030', 'Arabidopsis_gene': 'AT2G32030.1', 'AverageShortestPathLength': 6.935239018087855, 'BetweennessCentrality': 0.0072612772364046, 'ClosenessCentrality': 0.144191137088699, 'Degree': 9, 'Eccentricity': 17, 'Flow': 0.000336021, 'IsSingleNode': False, 'NeighborhoodConnectivity': 3.7777777777777777, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 9, 'Radiality': 0.9573004387187924, 'SelfLoops': 0, 'Stress': 64738426, 'TopologicalCoefficient': 0.1322751322751322, 'id': 'LOC111987461', 'isDEG': 'DEG', 'isTF': None, 'isTR': 'TR', 'label': 'LOC111987461', 'title': {}}, 'LOC111995138': {'Arabidopsis_concise': 'AT3G13790', 'Arabidopsis_gene': 'AT3G13790.1', 'AverageShortestPathLength': 8.749354005167959, 'BetweennessCentrality': 4.1402725288644974e-05, 'ClosenessCentrality': 0.1142941523922031, 'Degree': 2, 'Eccentricity': 20, 'Flow': 9.8771e-05, 'IsSingleNode': False, 'NeighborhoodConnectivity': 4.5, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 2, 'Radiality': 0.9442492517613816, 'SelfLoops': 0, 'Stress': 228026, 'TopologicalCoefficient': 0.5833333333333334, 'borderWidth': 2, 'color': {'background': '#f09d62', 'border': '#e3530b'}, 'id': 'LOC111995138', 'isDEG': None, 'isTF': None, 'isTR': None, 'label': 'LOC111995138', 'title': {}}, 'LOC111997034': {'Arabidopsis_concise': None, 'Arabidopsis_gene': None, 'AverageShortestPathLength': 6.966892764857882, 'BetweennessCentrality': 0.0021349176544175, 'ClosenessCentrality': 0.1435360114977166, 'Degree': 7, 'Eccentricity': 17, 'Flow': 0.000263602, 'IsSingleNode': False, 'NeighborhoodConnectivity': 4, 'NumberOfDirectedEdges': 0, 'NumberOfUndirectedEdges': 7, 'Radiality': 0.9570727139218858, 'SelfLoops': 0, 'Stress': 8016142, 'TopologicalCoefficient': 0.1666666666666666, 'id': 'LOC111997034', 'isDEG': 'DEG', 'isTF': None, 'isTR': None, 'label': 'LOC111997034', 'title': {}}}
pdf = PDFReport(['LOC1','loc2'])
pdf.add_page()

# Add content to PDF
pdf.add_summary(irp_score="0.8", tf_rank=["0", "1"])
pdf.add_legend()
#pdf.add_edges_table()
pdf.add_nodes_table(node_dict)

pdf.output('gene_relation_report.pdf')
'''
