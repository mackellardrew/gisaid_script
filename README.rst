=============
gisaid_script
=============


.. image:: https://img.shields.io/pypi/v/gisaid_script.svg
        :target: https://pypi.python.org/pypi/gisaid_script

.. image:: https://img.shields.io/travis/mackellardrew/gisaid_script.svg
        :target: https://travis-ci.com/mackellardrew/gisaid_script

.. image:: https://readthedocs.org/projects/gisaid-script/badge/?version=latest
        :target: https://gisaid-script.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




A package for formatting metadata for SARS-CoV-2 metadata


* Free software: GNU General Public License v3
* Documentation: https://gisaid-script.readthedocs.io.


Features
--------

The main focus of this repository is a "**gisaid_script**" to aid
submission of SARS-CoV-2 genome sequences prepared by the `WA DOH
PHL <https://www.doh.wa.gov/forpublichealthandhealthcareproviders/publichealthlaboratories>`__
to the `GISAID <https://www.gisaid.org/>`__ sequence repository, for use
in public health surveillance work during the COVID-19 pandemic.

The script requires a single positional arg as input: the GISAID
username of the person submitting the sequences. In addition, it
requires two input files to function:

1. A TSV/Excel dump of the latest sample inputs from the DOH "Dashboard"
   tracking sample statuses

   -  This file is needed primarily for the updated virus name field to
      identify the samples once published, but also for metadata
      required by GISAID

2. A TSV/Excel table of results from the
   `Terra <https://app.terra.bio/>`__ bioinformatics platform

   -  This table contains relevant QC metrics for assessing the
      suitability of genome publication, as well as links to the genome
      sequences within a Google Cloud Storage container

The absolute paths to these files can be provided using the ``--terra``
and ``--dashboard`` flags. Alternatively, the simplest way to use it is
to put both files in a common dir, as the only files containing "terra"
and "dashboard" in their respective filenames, and either providing the
path to this dir with ``--indir``, or simply running the script from
within that dir. More than one Terra table can be passed at a time, if
desired.

The script will download all genome sequences to a new "assemblies"
subdir of the current working dir. Major outputs are two files: a
``gisaid_metadata.csv`` file and an ``all_sequences.fa`` FASTA file
containing the sequence data.

Also included in this repository is a second, much simpler
"**terra_consolidate_script**", which is meant to aid combining
periodically the data tables produced by Terra workflows for individual
runs into a single, larger data table, to reduce clutter in the WA DOH
PHL Terra workspaces.

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
