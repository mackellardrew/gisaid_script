# gisaid_script

This is a script to aid submission of SARS-CoV-2 genome sequences prepared by the [WA DOH PHL](https://www.doh.wa.gov/forpublichealthandhealthcareproviders/publichealthlaboratories) to the [GISAID](https://www.gisaid.org/) sequence repository, for use in public health surveillance work during the COVID-19 pandemic.

The script requires a single positional arg as input: the GISAID username of the person submitting the sequences.  In addition, it requires two input files to function:

1. A TSV/Excel dump of the latest sample inputs from the DOH "Dashboard" tracking sample statuses
    * This file is needed primarily for the updated virus name field to identify the samples once published, but also for metadata required by GISAID
2. A TSV/Excel table of results from the [Terra](https://app.terra.bio/) bioinformatics platform
    * This table contains relevant QC metrics for assessing the suitability of genome publication, as well as links to the genome sequences within a Google Cloud Storage container
  
The absolute paths to these files can be provided using the `--terra` and `--dashboard` flags.  Alternatively, the simplest way to use it is to put both files in a common dir, as the only files containing "terra" and "dashboard" in their respective filenames, and either providing the path to this dir with `--indir`, or simply running the script from within that dir.  More than one Terra table can be passed at a time, if desired.

The script will download all genome sequences to a new "assemblies" subdir of the current working dir.  Major outputs are two files: a `gisaid_metadata.csv` file and an `all_sequences.fa` FASTA file containing the sequence data.
