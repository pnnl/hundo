Installation
============

This protocol leverages the work of Bioconda and depends on ``conda``.
For complete setup of these, please see:

https://bioconda.github.io/#using-bioconda

Really, you just need to make sure ``conda`` is executable and you've
set up your channels (steps 1 and 2). Then:

::

    conda install python=3.6 \
        pyyaml snakemake>=5.1.4 biopython \
        biom-format=2.1.6 numpy pandas \
        plotly=2.7.0
    pip install hundo

To update to the newest version of Hundo, run

::

    pip install --upgrade hundo

Alternatively, if you do not want any new executables in your environment
you can install into a new conda environment, e.g. hundo_env::

    conda create --name hundo_env python=3.6 \
        pyyaml snakemake>=5.1.4 biopython \
        biom-format=2.1.6 numpy pandas \
        plotly=2.7.0
    source activate hundo_env
    pip install hundo

To leave the environment::

    source deactivate
