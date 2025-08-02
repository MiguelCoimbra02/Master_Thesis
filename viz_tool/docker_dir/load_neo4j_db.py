#When adding new data always load the nodes data first
import argparse
import os
import sys
from neo4j_graph_service import Neo4jGraphService

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(
            description="Load node and edge data into a Neo4j database.",
            epilog="""
    Examples:

        Load nodes from a CSV file (located inside the Neo4j import folder):
            python load_neo4j_db.py --function load_nodes_to_neo4j

        Load edges from a CSV file (located inside the Neo4j import folder):
            python load_neo4j_db.py --function load_edges_to_neo4j
        
        Remove edges from a Neo4j database:
            python load_neo4j_db.py --function remove_edges_from_neo4j

    Options can be used to specify file names and directory:

        Example with custom file names and location:
            python load_neo4j_db.py --function load_nodes_to_neo4j \\
                --data_dir /var/lib/neo4j/import/custom_folder \\
                --nodes_file_name custom_nodes.csv

        Neo4j must have access to the file path used (typically: /var/lib/neo4j/import).
        Ensure the Neo4j config allows CSV loading from the given folder.
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Add arguments to the parser to choose the function to run
    parser.add_argument('--function', required=True, choices=[
        'load_nodes_to_neo4j',
        'load_edges_to_neo4j',
        'remove_edges_from_neo4j',
    ], help='Function to run.')

    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    #nodes options
    parser.add_argument('--data_dir', type=str, default=f'{os.path.join(base_dir, 'data', 'data_targeted', 'latest')}',help='Path to dir where the data lays.')
    parser.add_argument('--nodes_file_name', type=str, default='nodes.csv')
    parser.add_argument('--edges_file_name', type=str, default='edges.csv')

    
    # Parse the arguments
    args = parser.parse_args()

    # Instantiate the service
    neo4j_service = Neo4jGraphService()

    # Check which function to run based on the arguments
    # The default name for the file
    if args.function == "load_nodes_to_neo4j":
        neo4j_service.load_nodes_to_neo4j(
            args.data_dir,
            args.nodes_file_name
        )
    elif args.function == "load_edges_to_neo4j":
        neo4j_service.load_edges_to_neo4j(
            args.data_dir,
            args.edges_file_name
        )
    elif args.function == "remove_edges_from_neo4j":
        neo4j_service.remove_edges_from_neo4j(
            args.data_dir,
            args.edges_file_name
        )
    else:
        print("Please specify an action: --function load_nodes_to_neo4j (options) or --function load_edges_to_neo4j (options)")

if __name__ == "__main__":
    main()