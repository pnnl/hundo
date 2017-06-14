Configuration Settings
======================

Defaults are displayed in the example.

Author
------

Prints the author in the footer of the report.

::

    author: hundo (https://github.com/pnnl/hundo)

Threads
-------

Some codes within ``hundo`` can use multiple threads. This sets the maximum
number of threads to be used for any given sub-command of the protocol. The
default is the CPU count on the machine invoking ``hundo make-config``.

::

    threads: 1

Optionally, ``threads_small`` can be specified to further split quality
trimming jobs.

::

    threads_small: 1

Database Directory
------------------

Databases will be downloaded the first time the protocol is executed.
Subsequent runs will utilize already downloaded files. This parameter sets
where to store or where to look for the database reference files.

::

    database_dir: ./databases

Quality Control
---------------

Optionally, when performing quality trimming, the user can supply a FASTA
file containing adapter sequences to be trimmed from reads.

::

    adapters: ""

During quality filtering, the user may also supply a contaminant database,
e.g. PhiX, from which reads will be aligned and filtered.

::

    filter_contaminants: ""

Kmer mismatches allowed during the adapter trimming process is set with::

    allowable_kmer_mismatches: 1

The length of kmer to search against sequences::

    reference_kmer_match_length: 31

Look for shorter kmers at read tips down to this length (0 disables)::

    reduced_kmer_min: 8

Set the passing single-end read length, prior to merging with::

    minimum_passing_read_length: 100

Trim ends using a minimum base quality::

    minimum_base_quality: 10

Being more lenient on quality here is advised as qualities may improve after
merging R1 and R2. Poor reads should then be filtered out using
``maximum_expected_error``.

Merging
-------

Set the minimum merge length and the allowable error rate::

    minimum_merge_length: 150
    maximum_expected_error: 1

Chimeric Sequences
------------------

Optionally perform *de novo* chimera filtering::

    denovo_chimera_filter: true

To skip reference-based chimera filtering, the option can be omitted or set to
false. To perform reference-based filtering, one may supply either the path
to the desired database or set to true to use the preconfigured BLAST database::

    reference_chimera_filter: true

Clustering
----------

Set the minimum allowable sequence abundance. Set greater than 1 to dereplicate::

    minimum_sequence_abundance: 2


Set the maximum difference between an OTU member sequence and the
representative sequence of that OTU::

    percent_of_allowable_difference: 3

Annotation
----------

Perform error correction on OTU seed sequences prior to aligning reads for
counts::

    perform_error_correction: true

Select a reference database for annotation. For 16S, use either 'greengenes'
or 'silva', and for ITS use 'unite'::

    reference_database: silva

Annotation alignments are performed using BLAST followed by LCA. Hits may be
filtered by minimum bitscore::

    blast_minimum_bitscore: 100

Or, hits can be filtered by their bitscore relative to the best hit::

    blast_top_fraction: 0.95

When mapping reads back to OTUs for quantification, set their mapping
identity using::

    read_identity_requirement: 0.97
