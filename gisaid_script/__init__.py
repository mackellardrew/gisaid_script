import os
import pathlib
import sys

import gisaid_script.config as config
from gisaid_script.core import (auto_qc, download_assemblies,
                                format_collection_date, format_ictv_isolate,
                                generate_fasta,
                                get_biosample_metadata,
                                get_collecting_lab_address, get_column_map,
                                get_genbank_metadata, get_gisaid_metadata,
                                get_pha4ge_metadata, get_platform, get_vocs,
                                handle_counties,
                                handle_missing_data, handle_missing_genomes,
                                handle_reasons, handle_vocs,
                                load_tables, merge_tables, setup_logger)


def main():
    """Run the functions of this script in order, to process data in 
    preparation for uploading to GISAID"""
    print()
    pathlib.Path(config.OUTDIR).mkdir(exist_ok=True)
    logger = setup_logger(config.OUTDIR)
    col_names = get_column_map(config.WORKFLOW)
    dashboard_df = load_tables(config.DASHBOARD_TABLE)
    terra_df = load_tables(config.TERRA_TABLE, terra_table=True)

    for df, path in zip(
        (dashboard_df, terra_df), (config.DASHBOARD_TABLE, config.TERRA_TABLE)
    ):
        if df.shape[0] < 1:
            bad_input_message = (
                "Could not interpret or could not find appropriate data "
                f"in {path}; please check format of input file and try again."
            )
            logger.critical(bad_input_message)
            sys.exit()

    # DCM 20210805 note: the following block is just to get
    # a local copy of the `merged_df` table to look at
    merged_df = merge_tables(terra_df, dashboard_df, logger=logger)
    
    req_fields = ["collected_date"]
    samples_missing_data = handle_missing_data(merged_df, req_fields, logger)
    vocs, vois = get_vocs()
    _, _ = handle_vocs(vocs, vois, merged_df, col_names, logger)
    if not config.NO_AUTO_QC:
        bad_samples = auto_qc(merged_df, col_names, logger)
    else:
        bad_samples = []

    # Note: the assembly download execution step is the most
    # Costly & time-intensive; comment out lines below if they're already present
    print(
        "Downloading consensus genome assemblies "
        "from the cloud; this may take some time...",
        end="\n",
    )
    failed_samples = samples_missing_data + bad_samples

    if config.SKIP_ASSEMBLY_DOWNLOAD:
        download_stderrs = dict()
    else:
        _, download_stderrs = download_assemblies(
            merged_df[~merged_df["wa_no"].isin(failed_samples)]
        )

    missing_genomes = handle_missing_genomes(
        merged_df, download_stderrs, col_names, logger
    )
    failed_samples.extend(missing_genomes)
    fasta_generation_errs = generate_fasta(
        merged_df[~merged_df["wa_no"].isin(failed_samples)], col_names, logger
    )
    failed_samples.extend(list(fasta_generation_errs.keys()))
    pha4ge_metadata_df = get_pha4ge_metadata(
        merged_df[~merged_df["wa_no"].isin(failed_samples)],
        logger
    )#.set_index('specimen_collector_sample_id', drop=True)
    drop_reasons = pha4ge_metadata_df[
        pha4ge_metadata_df['purpose_of_sequencing'] == 'DROP'
    ]
    gisaid_metadata_df = (
        get_gisaid_metadata(pha4ge_metadata_df, logger)
        .droplevel(1, axis=1)
        .set_index("submitter")
    )
    # Next step is to generate genbank_metadata_df
    genbank_metadata_df = (
        get_genbank_metadata(pha4ge_metadata_df, logger)
    )

    biosample_metadata_df = get_biosample_metadata(pha4ge_metadata_df, logger)
    genbank_metadata_df = get_genbank_metadata(pha4ge_metadata_df, logger)

    pha4ge_outpath = os.path.join(config.OUTDIR, "pha4ge_metadata.csv")
    pha4ge_outpath, gisaid_outpath, biosample_outpath, genbank_outpath = (
        os.path.join(config.OUTDIR, filename) for filename in 
        ('pha4ge_metadata.csv', 'gisaid_metadata.csv', 'biosample_metadata.csv', 'genbank_metadata.csv')
    )

    # A brief block to provide any final formatting to metadata files
    for metadata_df, new_index, metadata_table_name, outpath in zip(
        (pha4ge_metadata_df, gisaid_metadata_df, biosample_metadata_df, genbank_metadata_df),
        ('specimen_collector_sample_id', 'submitter', 'sample_name', 'sequence_ID'),
        ("PHA4GE", "GISAID", "BioSample", "GenBank"),
        (pha4ge_outpath, gisaid_outpath, biosample_outpath, genbank_outpath),
    ):
        if metadata_table_name != "GISAID":
            metadata_df.set_index(new_index, drop=True, inplace=True)
        metadata_df.to_csv(outpath)
        logger.info(f"{metadata_table_name} metadata file written to {outpath}")

    print("Done", end="\n\n")

if __name__ == "__main__":
    main()
