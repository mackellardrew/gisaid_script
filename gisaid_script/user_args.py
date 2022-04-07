import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument(help="Your GISAID submitter ID", dest="submitter")
parser.add_argument(
    "-i",
    "--indir",
    help=(
        "Path to directory containing input tables "
        "(mutually exclusive with '--terra' and '--dashboard' options)"
    ),
    type=str,
    dest="indir",
    default=os.getcwd(),
)
parser.add_argument(
    "-t",
    "--terra",
    help=(
        "File path to the table of results from Terra; "
        "can be a list, with each item preceded by the "
        "'-t' or '--terra' flag."
    ),
    nargs="*",
    dest="terra_table",
    default=None,
)
parser.add_argument(
    "-d",
    "--dashboard",
    help="File path to the table of metadata from the WAPHL dashboard",
    nargs="*",
    dest="dashboard_table",
    default=None,
)
parser.add_argument(
    "-v",
    "--vocs",
    help=(
        "File path to a list of variants of concern/interest. "
        "File should have two columns: VOC and VOI, named "
        "on the first line, with Nextclade clade designations "
        "or Pango Lineages listed underneath for variants of "
        "concern and interest, respectively"
    ),
    type=str,
    dest="voc_list",
    default=None,
)
parser.add_argument(
    "-o",
    "--outdir",
    help=(
        "Path to directory to which outputs should be written; "
        "default is the '--indir', if provided, or else the "
        "current working dir"
    ),
    type=str,
    dest="outdir",
    default=os.getcwd(),
)
parser.add_argument(
    "-w",
    "--workflow",
    help=("Workflow used to generate Terra table: 'titan' or 'lang'"),
    type=str,
    dest="workflow",
    default="titan",
)
parser.add_argument(
    "-g",
    "--gsutil",
    help=("Absolute path to gsutil tool, in case not in active PATH"),
    type=str,
    dest="gsutil_path",
    default="gsutil",
)
parser.add_argument(
    "-s",
    "--skip",
    help=(
        "Skip downloading assemblies; assumes they're present in this dir "
        "in a subdir named 'assemblies'"
    ),
    type=bool,
    dest="skip_assembly_download",
    default=False,
)
parser.add_argument(
    "--no_auto_qc",
    help=("If TRUE, ignore genome QC criteria in generating outputs"),
    type=bool,
    dest="no_auto_qc",
    default=False,
)
parser.add_argument(
    "--author_list",
    help=(
        "Absolute path to file containing a list of Authors to be attributed "
        "in repositories.  Names should be on a single line, separated by semicolons."
    ),
    type=str,
    dest="author_list",
)

user_args = vars(parser.parse_args())