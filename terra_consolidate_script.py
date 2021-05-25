import argparse
import datetime
import os
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument(
    help=("Path to dir containing all input Terra tables"), type=Path, dest="indir",
)
parser.add_argument(
    "-e",
    "--exclude_controls",
    help=("If 'FALSE', don't attempt to scrub controls from the data tables"),
    type=bool,
    dest="exclude_controls",
    default=True,
)
user_args = vars(parser.parse_args())
INDIR = user_args.get("indir")
EXCLUDE_CONTROLS = user_args.get("exclude_controls")


def consolidate_terra(INDIR, exclude_controls=EXCLUDE_CONTROLS):
    """Gathers all TSV files in `INDIR`, concatenates into a 
    single output table, and writes to the CWD"""
    # Gather the tables together into a list
    filenames = [filename for filename in os.listdir(INDIR) if filename[-4:] == ".tsv"]
    for filename in filenames:
        print(f"Found Data table: {filename}")
    paths = [os.path.join(INDIR, filename) for filename in filenames]
    dfs = [pd.read_csv(path, sep="\t") for path in paths]

    # Change the first column name so Terra will rename the table
    today = datetime.datetime.now().strftime("%Y%m%d")
    newname = f"entity:{today}_RUNNING_TOTAL_id"
    for df in dfs:
        if newname in df.columns:
            dfs.remove(newname)
        cols = df.columns.copy().tolist()
        cols[0] = newname
        df.columns = cols

    # Gather into a single DataFrame and sort
    df = pd.concat(dfs)

    # If desired, exclude any non-WA_no samples
    if exclude_controls:
        pattern = "WA[0-9]{7}"
        matches = df[newname].dropna().str.contains(pattern)
        df = df[matches]
        # df.sort_values(newname)

    # Write the new table out
    outfilename = "consolidated_terra.tsv"
    outpath = os.path.join(os.getcwd(), outfilename)
    df.set_index(newname, drop=True, inplace=True)
    df.to_csv(outpath, sep="\t")

    print(f"Wrote consolidated Data table to {outpath}")


if __name__ == "__main__":
    consolidate_terra(INDIR)
