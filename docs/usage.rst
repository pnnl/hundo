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
