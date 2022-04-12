import argparse
import logging
import os
import pathlib
import re
import shlex
import subprocess
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from functools import partial
from glob import glob
from datetime import datetime
from IPython.display import display
from logging.handlers import RotatingFileHandler
from typing import Tuple

import gisaid_script.config as config

def setup_logger(output_dir: str) -> logging.Logger:
    """Returns a logger object writing to 'gisaid_script_logs.txt'."""
    log_filename = "gisaid_script_logs.txt"
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logger = logging.getLogger("gisaid_script_logger")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_filename)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.addHandler(logging.StreamHandler())

    return logger


def load_tables(table_list: list, terra_table=False) -> pd.DataFrame:
    """Load input tables and consolidate into pandas DataFrames"""

    def load_single_table(filepath):
        """Attempt to determine whether a given file is 
        TSV, CSV, or Excel, and return the given table as 
        a pandas DataFrame"""
        _, ext = os.path.splitext(filepath)
        df = config.EXTENSION_HANDLERS.get(ext, pd.read_excel)(filepath)
        return df

    df_list = list()
    for table in table_list:
        single_df = load_single_table(table)
        if terra_table:
            cols = single_df.columns.copy().tolist()
            cols[0] = "sample_name"
            single_df.columns = cols
        df_list.append(single_df)
    df = pd.concat(df_list)
    return df


def merge_tables(
    terra_df: pd.DataFrame, dashboard_df: pd.DataFrame, logger: logging.Logger
) -> pd.DataFrame:
    """Merges the Dashboard and Terra tables, and reformats slightly"""
    pattern = ".*(WA[0-9]{7}).*"
    terra_df["wa_no"] = terra_df["sample_name"].str.extract(pattern)
    missing_wa_nos = terra_df[terra_df["wa_no"].isna()]["sample_name"].tolist()
    if len(missing_wa_nos) > 0:
        missing_wa_nos_msg = "\n".join(
            [
                ("No WA number could be determined for the following " "samples: "),
                *missing_wa_nos,
                "and they will be omitted from the outputs",
            ]
        )
        logger.warning(missing_wa_nos_msg)
        print()
    terra_df.dropna(subset=["wa_no"], inplace=True)

    dashboard_df.columns = [
        "_".join(col.lower().split()) for col in dashboard_df.columns
    ]

    merged_df = pd.merge(
        terra_df, dashboard_df, left_on="wa_no", right_on="specimenid", how="left"
    )
    # "sample_id" field required by PHA4GE standard, later
    merged_df['sample_id'] = merged_df['seq_id'].fillna('').str.split('/', expand=True).iloc[:, 2]

    return merged_df


def get_column_map(workflow: str) -> dict:
    cols_needed = (
        "sequence",
        "coverage",
        "ivar_version",
        "nextclade_clade",
        "pangolin_lineage",
    )
    titan_col_names = (
        "assembly_fasta",
        "percent_reference_coverage",
        "ivar_version_consensus",
        "nextclade_clade",
        "pango_lineage",
    )
    lang_col_names = (
        "consensus_seq",
        "coverage_trim",
        "ivar_version_consensus",
        "nextclade_clade",
        "pangolin_lineage",
    )
    workflow_cols = dict()
    for key, col_names in zip(("titan", "lang"), (titan_col_names, lang_col_names)):
        workflow_cols[key] = {
            col_needed: col_name for col_needed, col_name in zip(cols_needed, col_names)
        }
    return workflow_cols[workflow]


def auto_qc(
    merged_df: pd.DataFrame, col_names, logger: logging.Logger
) -> list:
    conditions = merged_df[col_names.get("coverage")] < 60
    key_metrics = [col_names.get("coverage")]
    bad_samples = merged_df[conditions]["wa_no"].dropna().tolist()
    if len(bad_samples) > 0:
        auto_qc_msg = "\n".join(
            [
                (
                    "The following samples failed to meet the minimum QC metrics "
                    f"set for genome quality (currently {key_metrics})"
                ),
                *bad_samples,
                "and will be omitted from the outputs",
            ]
        )
        logger.warning(auto_qc_msg)
        print()
    return bad_samples


def download_assemblies(
    merged_df: pd.DataFrame, col_names
) -> Tuple[dict, dict]:
    """For each sample represented in the Terra results, attempts to 
    download the corresponding genome assembly"""

    pathlib.Path(config.ASSEMBLY_DIR).mkdir(exist_ok=True, parents=True)
    pipes = {key: subprocess.PIPE for key in ("stdout", "stderr")}
    download_stdouts, download_stderrs = dict(), dict()

    def gsutil_download(row):
        wa_no, url = row[0], row[1]
        cmd = f"{config.GSUTIL_PATH} cp {url} {config.ASSEMBLY_DIR}"
        proc = subprocess.Popen(shlex.split(cmd), **pipes)
        stdout, stderr = proc.communicate()
        download_stdouts.update({wa_no: stdout.decode("utf-8")})
        download_stderrs.update({wa_no: stderr.decode("utf-8")})

    _ = merged_df[["wa_no", col_names.get("sequence")]].apply(gsutil_download, axis=1)

    return download_stdouts, download_stderrs


def handle_counties(county: str):
    """Ensure that any fields reported for the County in which sample 
    was collected are validly named WA counties; else return just state 
    for localization field in metadata"""
    wa_counties_lower = (
        "adams; asotin; benton; chelan; clallam; clark; columbia; cowlitz; "
        "douglas; ferry; franklin; garfield; grant; grays harbor; island; "
        "jefferson; king; kitsap; kittitas; klickitat; lewis; lincoln; mason; "
        "okanogan; pacific; pend oreille; pierce; san juan; skagit; skamania; "
        "snohomish; spokane; stevens; thurston; wahkiakum; walla walla; "
        "whatcom; whitman; yakima"
    ).split("; ")
    no_county = "North America / USA / Washington"
    if county.lower() in wa_counties_lower:
        words = [word.capitalize() for word in county.split()]
        new_county = (" ".join(words) if len(words) > 1 else words[0])
        return_str = f"{no_county} / {new_county} County"
    else:
        return_str = no_county
    return return_str


def get_platform(sample_index: str) -> str:
    instrument_type = "Illumina NextSeq"
    miseqs = ["M4796", "M5130", "M5916"]
    for miseq in miseqs:
        if miseq in sample_index:
            instrument_type = "Illumina MiSeq"
    return instrument_type


def get_collecting_lab_address(lab_name: str) -> str:
    """Fetch the address of the submitting lab, if available"""
    lab_addr_map = {
        'atlas genomics': (
            '2296 W. Commodore Way, Suite 220, Seattle, WA 98199, USA'
        ),
        'confluence': (
            '1201 South Miller St Wenatchee, WA 98801, USA'
        ),
        'incyte diagnostics spokane': (
            '15912 East Marietta Avenue, Suite 200, '
            'Spokane Valley, WA 99216, USA'
        ),
        'interpath laboratory': (
            '8660 Emerald St # 102, Boise, ID 83704, USA'
        ),
        'northwest laboratory': (
            '3548 Meridian St, Suite 101, Bellingham, WA 98225, USA'
        ),
        'wa phl': (
            '1610 NE 150th St., Shoreline, WA 98155'
        ),
    }
    pha4ge_collecting_lab_addr = lab_addr_map.get(
        lab_name.lower(), 'WA, USA'
    )
    return pha4ge_collecting_lab_addr


def handle_reasons(reason: str) -> str:
    """Convert WAPHL categories of reasons for sequencing to those that fit PHA4GE standard"""
    reasons_map = {
        "phl diagnostic": "Baseline surveillance (random sampling)",
        "suspected vaccine breakthrough": "Vaccine escape surveillance",
        "suspected reinfection": "Re-infection surveillance",
        "outbreak investigation": "Cluster/Outbreak investigation",
        "travel associated": "Travel-associated surveillance",
        "other": "Not Provided",
        "s-dropout": "Screening for Variants of Concern (VOC)",
        "sentinel surveillance": "Baseline surveillance (random sampling)",
        "pt": "DROP",
    }
    pha4ge_reason = reasons_map.get(
        reason.lower(), "NO_REASON_FOUND"
    )
    return pha4ge_reason


def get_pha4ge_metadata(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """Generates the proper field names & formats for metadata
    compatible with the PHA4GE standard."""
    pha4ge_map = {
        'specimen_collector_sample_id': df["sample_id"],
        'bioproject_umbrella_accession': "PRJNA615625",
        'bioproject_accession': "PRJNA749781",
        'biosample_accession': None,
        'genbank_ena_ddbj_accession': None,
        'gisaid_accession': None,
        'gisaid_virus_name': df["seq_id"],
        'sample_collected_by': df["submittinglab"],
        'sample_collector_contact_address': (
            df['submittinglab'].fillna("").apply(get_collecting_lab_address)
        ),
        'sequence_submitted_by': "Washington State Department of Health Public Health Laboratories",
        'sequence_submitter_contact_address': "1610 NE 150th St., Shoreline, WA 98155",
        'sample_collection_date': df["collected_date"],
        'geo_loc_name_country': "USA",
        'geo_loc_name_state_province_territory': "Washington",
        'geo_loc_name_county_region': df["county"].fillna("").apply(handle_counties),
        'organism': "Severe acute respiratory syndrome coronavirus 2",
        'isolate': df[["sample_id", "collected_date"]].apply(
            lambda x: f"SARS-CoV-2/Homo sapiens/USA/{x[0]}/{x[1]}",
            axis=1
        ),
        'host_scientific_name': "Homo sapiens",
        'host_disease': "COVID-19",
        'purpose_of_sequencing': df["reason"].apply(handle_reasons),
        'sequencing_instrument': df["sample_name"].apply(get_platform),
        'sequencing_protocol_name': "Illumina COVIDSeq",
        'raw_sequence_data_processing_method': "FASTQC, Trimmomatic: quality and adapter trimming",
        'dehosting_method': "NCBI SRA human scrubber",
        'consensus_sequence_software_name': "iVar",
        'consensus_sequence_software_version': (
            df["ivar_version_consensus"].fillna(" . . Unknown").str.split(expand=True).iloc[:, -1]
        ),
        'breadth_of_coverage_value': df["percent_reference_coverage"].fillna("0").apply(lambda x: f"{int(x)}%"),
        'depth_of_coverage_value': df["assembly_mean_coverage"].fillna("0").apply(lambda x: f"{int(x)}x"),
        'bioinformatics_protocol': (
            "https://github.com/theiagen/public_health_viral_genomics/blob/main/workflows/wf_titan_illumina_pe.wdl"
        ),
    }

    pha4ge_metadata_df = pd.DataFrame(pha4ge_map)

    pha4ge_metadata_df.dropna(
        axis=0, subset=['specimen_collector_sample_id'], inplace=True
    )

    return pha4ge_metadata_df


def get_gisaid_metadata(pha4ge_metadata_df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """Configure output metadata spreadsheet for upload to GISAID repository 
    from input PHA4GE standard"""
    gisaid_fields = {
        ("submitter", "Submitter"): config.SUBMITTER,
        ("fn", "FASTA filename"): "all_sequences.fasta",
        ("covv_virus_name", "Virus name"): pha4ge_metadata_df["gisaid_virus_name"],
        ("covv_type", "Type"): "betacoronavirus",
        ("covv_passage", "Passage details/history"): "Original",
        ("covv_collection_date", "Collection date"): pha4ge_metadata_df["sample_collection_date"],
        ("covv_location", "Location"): pha4ge_metadata_df["geo_loc_name_county_region"],
        ("covv_add_location", "Additional location information"): None,
        ("covv_host", "Host"): "Human",
        ("covv_add_host_info", "Additional host information"): None,
        ("covv_sampling_strategy", "Sampling Strategy"): pha4ge_metadata_df['purpose_of_sequencing'],
        ("covv_gender", "Gender"): "unknown",
        ("covv_patient_age", "Patient age"): "unknown",
        ("covv_patient_status", "Patient status"): "unknown",
        ("covv_specimen", "Specimen source"): None,
        ("covv_outbreak", "Outbreak"): None,
        ("covv_last_vaccinated", "Last vaccinated"): None,
        ("covv_treatment", "Treatment"): None,
        ("covv_seq_technology", "Sequencing technology"): pha4ge_metadata_df["sequencing_instrument"],
        ("covv_assembly_method", "Assembly method"): (
            pha4ge_metadata_df['consensus_sequence_software_version'].apply(lambda x: f"iVar v{x}")
        ),
        ("covv_coverage", "Coverage"): pha4ge_metadata_df['depth_of_coverage_value'],
        (
            "covv_orig_lab",
            "Originating lab",
        ): pha4ge_metadata_df['sample_collected_by'],
        ("covv_orig_lab_addr", "Address"): (
            pha4ge_metadata_df['sample_collector_contact_address']
        ),
        (
            "covv_provider_sample_id", 
            "Sample ID given by originating laboratory"
        ): None,
        (
            "covv_subm_lab",
            "Submitting lab",
        ): "Washington State Department of Health Public Health Laboratories",
        ("covv_subm_lab_addr", "Address"): "1610 NE 150th St., Shoreline, WA 98155",
        ("covv_subm_sample_id", "Sample ID given by the submitting laboratory"): None,
        ("covv_authors", "Authors",): config.AUTHORS,
        ("covv_comment", "Comment"): None,
        ("comment_type", "Comment Icon"): None,
    }

    gisaid_metadata_df = (
        pd.DataFrame(gisaid_fields)
        .sort_values(("covv_virus_name", "Virus name"))
        .reset_index(drop=True)
    )
    # new_output_df.dropna(
    #     axis=0, subset=[("covv_virus_name", "Virus name")], inplace=True
    # )

    return gisaid_metadata_df


def format_collection_date(date_col: pd.Series) -> pd.Series:
    """Reformat collection dates from MM/DD/YYYY in LIMS Dashboard output to 
    YYYY-MM-DD format required by NCBI."""
    formatted_date_col = pd.to_datetime(date_col, format='%Y-%m-%d')
    return formatted_date_col


def format_ictv_isolate(specimen_col: pd.Series) -> pd.Series:
    """Reformat pha4ge specimen_collector_sample_id into 
    ISTC format isolate name"""
    formatted_specimen_col = specimen_col.apply(
        lambda x: f"SARS-CoV-2/Human/USA/{x}/{config.TODAY[:4]}"
    )
    return formatted_specimen_col


def get_biosample_metadata(pha4ge_metadata_df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """Configure output metadata spreadsheet for upload to NCBI BioSample 
    database from input PHA4GE standard"""
    biosample_fields = {
        'sample_name': pha4ge_metadata_df['specimen_collector_sample_id'],
        'bioproject_accession': pha4ge_metadata_df['bioproject_accession'],
        # 'umbrella_bioproject_accession': pha4ge_metadata_df['bioproject_umbrella_accession'],
        'organism': pha4ge_metadata_df['organism'],
        'collected_by': pha4ge_metadata_df['sample_collected_by'],
        'collection_date': format_collection_date(pha4ge_metadata_df['sample_collection_date']),
        'geo_loc_name': (
            pha4ge_metadata_df['geo_loc_name_county_region'].str.split(' / ', expand=True).iloc[:, 1:]
            .fillna('').agg(': '.join, axis=1).str.strip(': ')
        ),
        'host': pha4ge_metadata_df['host_scientific_name'],
        'host_disease': pha4ge_metadata_df['host_disease'],
        'isolate': format_ictv_isolate(
            pha4ge_metadata_df['specimen_collector_sample_id']
        ),
        'isolation_source': 'Clinical/Nasal Swab',
        # Going to have to leave out the GISAID accession field for now;
        # can't figure out how to easily accommodate that within this script
        # 'gisaid_accession': pha4ge_metadata_df['gisaid_accession'],
        'gisaid_virus_name': pha4ge_metadata_df['gisaid_virus_name'],
        'purpose_of_sequencing': pha4ge_metadata_df['purpose_of_sequencing'],
        'sequenced_by': pha4ge_metadata_df['sequence_submitted_by'],
    }
    biosample_metadata_df = pd.DataFrame(biosample_fields)
    biosample_metadata_df.rename(columns={
        'umbrella_bioproject_accession': 'bioproject_accession',
    },
    inplace=True
    )
    return biosample_metadata_df


def get_genbank_metadata(pha4ge_metadata_df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """Configure output metadata spreadsheet for upload to GenBank repository 
    from input PHA4GE standard"""
    genbank_fields = {
        'sequence_ID': pha4ge_metadata_df['specimen_collector_sample_id'],
        'isolate': format_ictv_isolate(
            pha4ge_metadata_df['specimen_collector_sample_id']
        ),
        'collection-date': format_collection_date(
            pha4ge_metadata_df['sample_collection_date']
        ),
        'host': 'Homo sapiens',
        'country': 'USA: Washington',
        'isolation-source': 'Nasal swab',
        'BioProject Accession': pha4ge_metadata_df['bioproject_accession'],
        'BioSample Accession': pha4ge_metadata_df['biosample_accession'],
    }
    genbank_metadata_df = pd.DataFrame(genbank_fields)
    return genbank_metadata_df



def generate_fasta(df: pd.DataFrame, col_names: dict, logger: logging.Logger):
    """Gather assemblies and output with new header lines"""
    file_df = (
        df[["wa_no", "seq_id", col_names.get("sequence")]]
        .copy()
        .sort_values("seq_id")
    )

    seq_id_pattern = r".*/(?:call-consensus|cacheCopy)/(.*)"
    file_df["consensus_file"] = (
        config.ASSEMBLY_DIR
        + os.sep
        + file_df[col_names.get("sequence")].fillna("").str.extract(seq_id_pattern)
    )
    fasta_generation_errs = dict()

    def gather_seqs(row, out_buffer):
        try:
            with open(row["consensus_file"], "r") as seq_buffer:
                rec = next(SeqIO.parse(seq_buffer, "fasta"))
                rec.id = row["seq_id"]
                rec.description = ""
                SeqIO.write(rec, out_buffer, "fasta")
        except (TypeError, AttributeError, FileNotFoundError) as err:
            fasta_generation_errs[row["wa_no"]] = err
            pass

    all_seq_path = os.path.join(config.OUTDIR, "all_sequences.fa")
    with open(all_seq_path, "w") as out_buffer:
        _ = file_df.dropna(subset=["wa_no"]).apply(
            gather_seqs, axis=1, out_buffer=out_buffer
        )
    logger.info(f"Consolidated genome assemblies written to {all_seq_path}")
    if len(fasta_generation_errs) > 0:
        fasta_generation_msg = "\n".join(
            [
                (
                    "There were problems with gathering/renaming the "
                    "genome sequences for the following samples, and they "
                    "were omitted from the outputs:"
                ),
                *fasta_generation_errs.keys(),
            ]
        )
        logger.warning(fasta_generation_msg)
        print()
    return fasta_generation_errs


def handle_missing_data(
    df: pd.DataFrame, req_fields: list, logger: logging.Logger
):
    samples_missing_data = list()
    working_df = df.copy().set_index("wa_no")[req_fields]  # .astype(str)
    for req_field in req_fields:
        missing_mask = working_df[req_field].replace("", None).isna()
        samples_missing_data.extend(working_df[missing_mask].index.astype(str).tolist())
    if len(samples_missing_data) > 0:
        missing_data_msg = "\n".join(
            [
                (
                    "The following samples are missing data in the required "
                    f"fields {req_fields}:"
                ),
                *samples_missing_data,
                ("and will be omitted from the outputs."),
            ]
        )
        logger.warning(missing_data_msg)
        print()
    return samples_missing_data


def handle_missing_genomes(
    df: pd.DataFrame, col_names: dict, 
    download_stderrs: dict, logger: logging.Logger
):
    """Detects errors that occur when downloading genomes, 
    and prints/logs messages warning the user they will be 
    omitted from the outputs."""

    download_failures = list()
    download_failure_msgs = list()

    for sample, msg in download_stderrs.items():
        if ("AccessDeniedException" in msg) or ("CommandException" in msg):
            download_failures.append(sample)
            download_failure_msgs.append(msg)

    if len(download_failures) > 0:
        download_failures_msg = "\n".join(
            [
                "NOTE: the following sequences could not "
                "be downloaded successfully and will be omitted "
                "from the outputs: ",
                "\n".join(sample for sample in download_failures),
            ]
        )
        logger.warning(download_failures_msg)
        print()
    failure_msgs = " ".join(download_failure_msgs)
    if "AccessDeniedException" in failure_msgs:
        access_msg = (
            "Cannot access the stored object(s); "
            "please run 'gsutil config' before re-running "
            "script."
        )
        logger.warning(access_msg)
        print()
    elif ("CommandException" in failure_msgs) or ("Error" in failure_msgs):
        url_msg = (
            f"Please check formatting of '{col_names.get('sequence')}' "
            "column in input Terra tables for samples "
            "listed above."
        )
        logger.warning(url_msg)
        print()

    return download_failures


def get_vocs():
    if config.VOC_LIST:
        _, voc_ext = os.path.splitext(config.VOC_LIST)
        voc_df = config.EXTENSION_HANDLERS.get(voc_ext)(config.VOC_LIST)
        voc_values = voc_df.T.values
        vocs, vois = voc_values[0], voc_values[1]
    else:
        vocs, vois = [], []
    return vocs, vois


def handle_vocs(
    vocs: list, vois: list, terra_df: pd.DataFrame, 
    col_names: dict, logger: logging.Logger
):
    if len(vocs) == 0 and len(vois) == 0:
        return [], []
    clades = (
        terra_df[
            [
                "wa_no",
                col_names.get("nextclade_clade"),
                col_names.get("pangolin_lineage"),
            ]
        ]
        .dropna()
        .set_index("wa_no")[
            [col_names.get("nextclade_clade"), col_names.get("pangolin_lineage")]
        ]
    )
    voc_samples = clades[
        (clades[col_names.get("nextclade_clade")].isin(vocs))
        | (clades[col_names.get("pangolin_lineage")].isin(vocs))
    ]
    voi_samples = clades[
        (clades[col_names.get("nextclade_clade")].isin(vois))
        | (clades[col_names.get("pangolin_lineage")].isin(vois))
    ]

    vocs_msg, vois_msg = None, None
    for sample_df, label, list_, msg in zip(
        (voc_samples, voi_samples), ("VOC", "VOI"), (vocs, vois), (vocs_msg, vois_msg)
    ):
        if sample_df.shape[0] > 0:
            samples = sample_df.index.values.tolist()
            msg = "\n".join(
                [
                    (
                        "The following samples were found to be in the designated "
                        f"{label} list ({list_}): "
                    ),
                    *samples,
                    (
                        "Please notify the Epidemiologist group at "
                        "'wgs-epi@doh.wa.gov' prior to upload to GISAID"
                    ),
                ]
            )
            logger.info(msg)
            print()
    # Kind of janky solution to show which sample is which lineage:
    # display in terminal, then write to separate TSV file
    vocs_vois_df = pd.concat([voc_samples, voi_samples])
    if vocs_vois_df.shape[0] > 0:
        outpath = os.path.join(config.OUTDIR, "vocs_vois_table.tsv")
        vocs_vois_df.to_csv(outpath, sep="\t")
        vocs_vois_out_msg = (
            "The following table of samples "
            "and Pango Linage/NextClade Clade "
            f"was written to {outpath}."
        )
        logger.info(vocs_vois_out_msg)
        display(vocs_vois_df)
        print()

    return voc_samples, voi_samples