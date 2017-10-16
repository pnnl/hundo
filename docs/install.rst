Installation
============

This protocol leverages the work of Bioconda and depends on ``conda``.
For complete setup of these, please see:

https://bioconda.github.io/#using-bioconda

Really, you just need to make sure ``conda`` is executable and you've
set up your channels (numbers 1 and 2). Then:

::

    conda install python=3.5 pyyaml snakemake biopython biom-format=2.1.5
    pip install hundo
