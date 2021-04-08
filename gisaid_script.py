#! /usr/bin/python

import argparse
import datetime
import os
import pathlib
import re
import shlex
import subprocess
import sys
import pandas as pd
from Bio import SeqIO
from functools import partial
from glob import glob
from tqdm import tqdm


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
    type=argparse.FileType("r"),
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
    "-g",
    "--gsutil",
    help=("Absolute path to gsutil tool, in case not in active PATH"),
    type=str,
    dest="gsutil_path",
    default="gsutil",
)

user_args = vars(parser.parse_args())
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

# For additional ease of use, attempt to determine whether TSV, CSV, or Excel
extension_handlers = {
    ".csv": pd.read_csv,
    ".tsv": partial(pd.read_csv, sep="\t"),
    ".txt": partial(pd.read_csv, sep="\t"),
    ".xls": pd.read_excel,
    ".xlsx": pd.read_excel,
}


def infer_input_format(filepath):
    _, ext = os.path.splitext(filepath)
    df = extension_handlers.get(ext, pd.read_excel)(filepath)
    return df


if VOC_LIST:
    _, voc_ext = os.path.splitext(VOC_LIST)
    voc_df = extension_handlers.get(voc_ext)(VOC_LIST)
else:
    voc_df = None


def load_tables(table_list, terra_table=False):
    df_list = list()
    for table in table_list:
        df = infer_input_format(table)
        if terra_table:
            cols = df.columns.copy().tolist()
            cols[0] = "sample_name"
            df.columns = cols
        df_list.append(df)
    df = pd.concat(df_list)
    return df


dashboard_df = load_tables(DASHBOARD_TABLE)
terra_df = load_tables(TERRA_TABLE, terra_table=True)

for df, path in zip((dashboard_df, terra_df), (DASHBOARD_TABLE, TERRA_TABLE)):
    if df.shape[0] < 1:
        bad_input_message = (
            f"Could not interpret or could not find appropriate data "
            "in {path}; please check format of input file and try again."
        )
        sys.exit()

pattern = ".*(WA[0-9]{7}).*"
terra_df["wa_no"] = terra_df["sample_name"].str.extract(pattern)

dashboard_df.columns = ["_".join(col.lower().split()) for col in dashboard_df.columns]

merged_df = pd.merge(
    terra_df, dashboard_df, left_on="wa_no", right_on="folderno", how="left"
)

# Download Assemblies
assembly_dir = os.path.join(OUTDIR, "assemblies")
pathlib.Path(assembly_dir).mkdir(exist_ok=True)
pipes = {key: subprocess.PIPE for key in ("stdout", "stderr")}
stdouts, stderrs = dict(), dict()


def gsutil_download(row):
    wa_no, url = row[0], row[1]
    cmd = f"{GSUTIL_PATH} cp {url} {assembly_dir}"
    proc = subprocess.Popen(shlex.split(cmd), **pipes)
    stdout, stderr = proc.communicate()
    stdouts.update({wa_no: stdout.decode("utf-8")})
    stderrs.update({wa_no: stderr.decode("utf-8")})


def handle_counties(county):
    valid_wa_counties = "'adams; asotin; benton; chelan; clallam; clark; columbia; cowlitz; douglas; ferry; franklin; garfield; grant; grays harbor; island; jefferson; king; kitsap; kittitas; klickitat; lewis; lincoln; mason; okanogan; pacific; pend oreille; pierce; san juan; skagit; skamania; snohomish; spokane; stevens; thurston; wahkiakum; walla walla; whatcom; whitman; yakima'".split(
        "; "
    )

    if county.lower() in valid_wa_counties:
        return_str = f"North America / USA / Washington / {county.capitalize()}"
    else:
        return_str = f"North America / USA / Washington"
    return return_str


def prep_metadata(df):
    # Configure output metadata spreadsheet
    new_fields = {
        ("submitter", "Submitter"): SUBMITTER,
        ("fn", "FASTA filename"): "all_sequences.fasta",
        ("covv_virus_name", "Virus name"): df["seq_id"],
        ("covv_type", "Type"): "betacoronavirus",
        ("covv_passage", "Passage details/history"): "Original",
        ("covv_collection_date", "Collection date"): df["collected_date"],
        ("covv_location", "Location"): df["county"].fillna("").apply(handle_counties),
        ("covv_add_location", "Additional location information"): None,
        ("covv_host", "Host"): "Human",
        ("covv_add_host_info", "Additional host information"): None,
        ("covv_sampling_strategy", "Sampling Strategy"): None,
        ("covv_gender", "Gender"): "unknown",
        ("covv_patient_age", "Patient age"): "unknown",
        ("covv_patient_status", "Patient status"): "unknown",
        ("covv_specimen", "Specimen source"): None,
        ("covv_outbreak", "Outbreak"): None,
        ("covv_last_vaccinated", "Last vaccinated"): None,
        ("covv_treatment", "Treatment"): None,
        ("covv_seq_technology", "Sequencing technology"): "Illumina MiSeq",
        ("covv_assembly_method", "Assembly method"): df["ivar_version_consensus"],
        ("covv_coverage", "Coverage"): df["depth_trim"],
        (
            "covv_orig_lab",
            "Originating lab",
        ): "Washington State Department of Health Public Health Laboratories",
        ("covv_orig_lab_addr", "Address"): "1610 NE 150th St., Shoreline, WA 98155",
        ("covv_provider_sample_id", "Sample ID given by originating laboratory"): None,
        (
            "covv_subm_lab",
            "Submitting lab",
        ): "Washington State Department of Health Public Health Laboratories",
        ("covv_subm_lab_addr", "Address"): "1610 NE 150th St., Shoreline, WA 98155",
        ("covv_subm_sample_id", "Sample ID given by the submitting laboratory"): None,
        ("covv_authors", "Authors",): (
            "Drew MacKellar, Philip Dykema, Denny Russell, "
            "Joenice Gonzalez, Hannah Gray, Geoff Melly, "
            "Vanessa De Los Santos, Darren Lucas, JohnAric Peterson, "
            "Avi Singh, Rebecca Cao"
        ),
        ("covv_comment", "Comment"): None,
        ("comment_type", "Comment Icon"): None,
    }

    new_output_df = (
        pd.DataFrame(new_fields)
        .sort_values(("covv_virus_name", "Virus name"))
        .reset_index(drop=True)
    )
    new_output_df.dropna(
        axis=0, subset=[("covv_virus_name", "Virus name")], inplace=True
    )

    return new_output_df


# Gather assemblies and output with new header lines
file_df = (
    merged_df[["wa_no", "seq_id", "consensus_seq"]]
    .copy()
    .sort_values("seq_id")
    # .dropna(axis=0, subset=["seq_id", "consensus_seq"])
)

new_pattern = r".*/(?:call-consensus|cacheCopy)/(.*)"
file_df["consensus_file"] = (
    assembly_dir + os.sep + file_df["consensus_seq"].str.extract(new_pattern)
)


def gather_seqs(row, out_buffer):
    # wa_no, seq_id, seq_path = row['wa_no'], row["seq_id"], row["consensus_file"]
    try:
        with open(row["consensus_file"], "r") as seq_buffer:
            rec = next(SeqIO.parse(seq_buffer, "fasta"))
            rec.id = row["seq_id"]
            rec.description = ""
            # print(rec.description)
            SeqIO.write(rec, out_buffer, "fasta")
    except (TypeError, FileNotFoundError) as err:
        stderrs[row["wa_no"]] = err


def stderr_handling():
    failures = [
        sample
        for sample, msg in stderrs.items()
        if ("Exception" in msg) or ("Error" in msg)
    ]
    if len(failures) > 0:
        print(
            "NOTE: the following sequences could not "
            "be downloaded successfully and will be omitted "
            "from the outputs:",
            "\n".join(sample for sample in failures),
            sep="\n",
        )
    stderrs_msgs = " ".join(stderrs.values())
    if "AccessDeniedException" in stderrs_msgs:
        access_msg = (
            "NOTE: Cannot access the stored object(s); "
            "please run 'gsutil config' before re-running "
            "script."
        )
        print(access_msg)
    elif ("CommandException" in stderrs_msgs) or ("Error" in stderrs_msgs):
        url_msg = (
            "Please check formatting of 'consensus_seq' "
            "column in input Terra tables for samples "
            "listed above."
        )
        print(url_msg)

    return failures


def main():
    # Note: the assembly download execution step is the most
    # Costly & time-intensive; comment out line below if they're already present
    print(
        "Downloading consensus genome assemblies "
        "from the cloud; this may take some time..."
    )
    tqdm.pandas()
    _ = terra_df[["wa_no", "consensus_seq"]].progress_apply(gsutil_download, axis=1)
    not_downloaded = stderr_handling()
    new_df = merged_df[~merged_df["wa_no"].isin(not_downloaded)]
    new_output_df = prep_metadata(new_df)

    outpath = os.path.join(OUTDIR, "gisaid_metadata.csv")
    new_out_df = new_output_df.droplevel(1, axis=1)
    new_out_df.set_index("submitter").to_csv(outpath)
    print(f"GISAID metadata file written to {outpath}")

    all_seq_path = os.path.join(OUTDIR, "all_sequences.fa")
    with open(all_seq_path, "w") as out_buffer:
        _ = file_df.apply(gather_seqs, axis=1, out_buffer=out_buffer)
    print(f"Consolidated genome assemblies written to {all_seq_path}")

    print("Done")


if __name__ == "__main__":
    main()

