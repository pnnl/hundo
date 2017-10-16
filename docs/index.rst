.. hundo


.. figure:: _static/logo.png
   :alt: logo


Snakemake-based amplicon processing protocol for 16S and ITS sequences:

-  Performs quality control based on quality, can trim adapters, and
   remove sequences matching a contaminant database;
-  Handles paired-end read merging;
-  Integrates *de novo* and reference-based chimera filtering;
-  Clusters sequences and annotates using databases that are downloaded
   as needed;
-  Generates standard outputs for these data like a newick tree, a
   tabular OTU table with taxonomy, and .biom.

This workflow is built using
`Snakemake <https://snakemake.readthedocs.io/en/stable/>`__ and makes
use of `Bioconda <https://bioconda.github.io/>`__ to install its
dependencies.

|DOI|

.. |DOI| image:: https://zenodo.org/badge/83449413.svg
   :target: https://zenodo.org/badge/latestdoi/83449413


Contents:
---------

.. toctree::
   :maxdepth: 2

   install
   usage
   parameters
   output
