from gisaid_script.user_args import user_args
import sys
import os
import pandas as pd
from functools import partial
from glob import glob
from datetime import datetime
# from IPython.display import display
# from logging.handlers import RotatingFileHandler
# from typing import Tuple

INDIR = user_args.get("indir")

# For ease of use, look for the files in INDIR
for key in ("terra_table", "dashboard_table"):
    if user_args[key] is None:
        key_part = key.split("_")[0]
        keypattern = f"*{key_part}*"
        matches = glob(os.path.join(INDIR, keypattern))
        if len(matches) > 0:
            user_args[key] = matches
        else:
            missing_input_message = (
                f"Could not find required input {key_part}. "
                f"Please either provide with --{key_part} flag, "
                f"or ensure that file with '{key_part}' in file name "
                f"is present in {INDIR}."
            )
            print(missing_input_message)
            sys.exit()


SUBMITTER = user_args.get("submitter")
TERRA_TABLE = user_args.get("terra_table")
DASHBOARD_TABLE = user_args.get("dashboard_table")
OUTDIR = user_args.get("outdir")
VOC_LIST = user_args.get("voc_list")
GSUTIL_PATH = user_args.get("gsutil_path")
NO_AUTO_QC = user_args.get("no_auto_qc")
WORKFLOW = user_args.get("workflow").lower()
SKIP_ASSEMBLY_DOWNLOAD = user_args.get('skip_assembly_download')

ASSEMBLY_DIR = os.path.join(OUTDIR, "assemblies")
EXTENSION_HANDLERS = {
    ".csv": pd.read_csv,
    ".tsv": partial(pd.read_csv, sep="\t"),
    ".txt": partial(pd.read_csv, sep="\t"),
    ".xls": partial(pd.read_excel, engine="xlrd"),
    ".xlsx": partial(pd.read_excel, engine="openpyxl"),
}
TODAY = datetime.now().strftime("%Y-%m-%d")
DEFAULT_AUTHORS = (
            "Drew MacKellar, Philip Dykema, Denny Russell, "
            "Holly Halstead, "
            "Kathryn Sickles, Kristin Roche, Ardizon Valdez, "
            "Claire Howell, Alex Lathum, JohnAric Peterson, "
            "Avi Singh, Rebecca Cao"
) 
AUTHORS_PATH = user_args.get('author_list')
if AUTHORS_PATH:
    with open(AUTHORS_PATH, 'r') as f:
        AUTHORS = [name.strip() for name in f.readline().split(';')]
else:
    AUTHORS = DEFAULT_AUTHORS