import argparse
import sys
import time
import pkg_resources

from ecmapper.create_database import download_database
from ecmapper.ecmapper import workflow

MESSAGE = """ecmapper -g FILE -d FOLDER -o FOLDER
ecmapper -d FOLDER"""

def cli():
    """Console script for ecmapper."""
    start_time = time.time()
    parser = argparse.ArgumentParser(
        "ecmapper",
        description=MESSAGE + "\nFor specific help on each subcommand use: ecmapper -help", formatter_class=argparse.RawTextHelpFormatter
    )
    # Parent parsers
    parent_parser_g = argparse.ArgumentParser(add_help=False)
    parent_parser_g.add_argument(
        "-g",
        "--gbk",
        help="genbank files",
        required=True,
    )
    parent_parser_d = argparse.ArgumentParser(add_help=False)
    parent_parser_d.add_argument(
        "-d",
        "--database",
        required=True,
        help="database path",
        metavar="OUPUT_DIR"
    )
    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        "-o",
        "--output",
        help="output folder",
        required=True,
        type=str
    )
   # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest="cmd")
    database_parser = subparsers.add_parser(
        "database",
        help="downlaod database use by ecmapper",
        parents=[
            parent_parser_d
        ],
        description=
        "downlaod database use by ecmapper"
    )
    mapper_parser = subparsers.add_parser(
        "map",
        help="map EC to bigg and modelseed",
        parents=[
            parent_parser_g, parent_parser_d, parent_parser_o
        ],
        description="map EC to bigg and modelseed"
    )
    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()

    if args.cmd == "database":
        database_folder = args.database
        download_database(database_folder)
    if args.cmd == "map":
        gbk_file = args.gbk
        database_folder = args.database
        output_folder = args.output
        workflow(gbk_file, database_folder, output_folder)
