|DOI| |Documentation Status| |bioconda-badge|

-  Performs quality control based on quality, can trim adapters, and
   remove sequences matching a contaminant database
-  Handles paired-end read merging
-  Integrates *de novo* and reference-based chimera filtering
-  Clusters sequences and annotates using databases that are downloaded
   as needed
-  Generates standard outputs for these data like a newick tree, a
   tabular OTU table with taxonomy, and .biom.

This workflow is built using
`Snakemake <https://snakemake.readthedocs.io/en/stable/>`__ and makes
use of `Bioconda <https://bioconda.github.io/>`__ to install its
dependencies.

Documentation
=============

For complete documentation and install instructions, see:

https://hundo.readthedocs.io

Install
=======

This protocol leverages the work of Bioconda and depends on ``conda``.
For complete setup of these, please see:

https://bioconda.github.io/#using-bioconda

Really, you just need to make sure ``conda`` is executable and you've
set up your channels (numbers 1 and 2). Then:

::

    conda install python>=3.6 click \
        pyyaml snakemake>=5.1.4 biopython
    pip install hundo


Usage
=====

Running samples through annotation requires that input FASTQs be
paired-end, named in a semi-conventional style starting sample ID,
contain "\_R1" (or "\_r1") and "\_R2" (or "\_r2") index identifiers, and
have an extension ".fastq" or ".fq". The files may be gzipped and end
with ".gz". By default, both R1 and R2 need to be larger than 10K in
size. This cutoff is arbitrary and can be set using
``--prefilter-file-size``.

Using the example data of the mothur SOP located in our tests directory,
we can annotate across SILVA using:

::

    cd example
    hundo annotate \
        --filter-adapters qc_references/adapters.fa.gz \
        --filter-contaminants qc_references/phix174.fa.gz \
        --out-dir mothur_sop_silva \
        --database-dir annotation_references \
        --reference-database silva \
        mothur_sop_data

Dependencies are installed by default in the results directory defined
on the command line as ``--out-dir``. If you want to re-use dependencies
across many analyses and not have to re-install each time you update the
output directory, use Snakemake's ``--conda-prefix``:

::

    hundo annotate \
        --out-dir mothur_sop_silva \
        --database-dir annotation_references \
        --reference-database silva \
        mothur_sop_data \
        --conda-prefix /Users/brow015/devel/hundo/example/conda

Output
======

**OTU.biom**

Biom table with raw counts per sample and their associated taxonomic
assignment formatted to be compatible with downstream tools like
phyloseq.

**OTU.fasta**

Representative DNA sequences of each OTU.

**OTU.tree**

Newick tree representation of aligned OTU sequences.

**OTU.txt**

Tab-delimited text table with columns OTU ID, a column for each sample,
and taxonomy assignment in the final column as a comma delimited list.

**OTU\_aligned.fasta**

OTU sequences after alignment using Clustal Omega.

**all-sequences.fasta**

Quality-controlled, dereplicated DNA sequences of all samples. The
header of each record identifies the sample of origin and the count
resulting from dereplication.

**blast-hits.txt**

The BLAST assignments per OTU sequence.

**summary.html**

Captures and summarizes data of the experimental dataset. Things like
sequence quality, counts per sample at varying stages of pre-processing, and
summarized taxonomic composition per sample across phylum, class, and
order.

Citing
======

Cite this as:

    Brown J, Zavoshy N, Brislawn CJ, McCue LA. (2018) Hundo: a Snakemake workflow for microbial community sequence data. PeerJ Preprints 6:e27272v1 https://doi.org/10.7287/peerj.preprints.27272v1


.. |DOI| image:: https://zenodo.org/badge/83449413.svg
   :target: https://zenodo.org/badge/latestdoi/83449413
.. |Documentation Status| image:: https://readthedocs.org/projects/hundo/badge/?version=latest
   :target: http://hundo.readthedocs.io/en/latest/?badge=latest
.. |bioconda-badge| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
   :target: http://bioconda.github.io
