![logo](resources/logo.png)

[![DOI](https://zenodo.org/badge/83449413.svg)](https://zenodo.org/badge/latestdoi/83449413)
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io)


# Install

The simplest way to do this is to grab the Python 3 version of [Anaconda](https://www.continuum.io/downloads), then add the prerequisite channels:

```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
```

And install dependencies

```
conda install python=3.5 snakemake biom-format fasttree bbmap blast clustalo
```

Unfortunately, [USEARCH](http://www.drive5.com/usearch/download.html) is currently a dependency of this protocol and it is not easily distributed.

Due to the size of this repo, as reference data is included, I'm hesitant to place this within Bioconda for now.

Finally, clone this code:

```
git clone git@github.com:pnnl/hundo.git
```

# Usage

## Experimental Data

Place demultiplexed, uncompressed reads into `results/<EXPERIMENT NAME>/demux`
with a `.fastq` file extension.

For `test-experiment` will look like:

```
cd hundo
tree results
└── test-experiment
    └── demux
        ├── sample-1_R1.fastq
        ├── sample-1_R2.fastq
        ├── sample-2_R1.fastq
        └── sample-2_R2.fastq
```

It's important that sample names end with some form of underscore then read identifier, e.g. _R1.fastq or _r1.fastq.

Sample names should not include underscores.

## Preparing Databases

Uncompress the BLAST database and taxonomy file:

```
cd ref/silva_123
sh extract.sh
gunzip SLV_123_SSU.tax.gz
```

## Prepare Other Executables

Build `lca`:

```
cd resources/lca_src
make
```

## Executing the Workflow

To run the workflow across 24 cores:

```
cd hundo
snakemake --configfile resources/16s.config.yaml --config eid=test-experiment
```

`eid` is our experiment ID and represents how we structure our results directory. It can also be specified in the configuration file by adding:

```
eid: test-experiment
```

Running the same command with `eid` defined in the configuration file would look like:

```
snakemake --configfile our-new-config.yaml
```

# Results

Using the above example, our results will be written to:

```
results/test-experiment
```

Methods and summary data can be found in:

```
results/test-experiment/97/blast/README.html
```

The summary portion of an example `README.html`:

![readme](resources/readme_summary.png)

# Configuration

See sample configuration files in resources for defaults. Reference database information has been pre-filled to an extent and encompasses databases included in this repository.

Notable options are:

```
annotation_method: blast
```

This performs alignment using BLAST. The alternative, if you prefer to use USEARCH, is:

```
annotation_method: utax
```

Reference-based chimera filtering is optional and can be disabled with:

```
chimera_filter_seed_sequences: false
```

Using both `utax` and `blast` on the same experiment will not overwrite any existing data. Output will be available for both methods.

Altering `percent_of_allowable_difference` also write an entirely new output directory under which analyses are performed.

This looks like:

```
test-experiment/
    ├── 97                                      # clustering pairwise identity threshold
    │   ├── blast
    │   │   ├── blast_hits.txt                  # raw blast hits per OTU seed seq
    │   │   ├── lca_assignments.txt             # raw lca results TSV from blast hits
    │   │   ├── OTU.biom                        # tax annotated biom (no metadata, no normalization)
    │   │   ├── OTU_tax.fasta                   # otu seqs with tax in FASTA header
    │   │   ├── OTU.txt                         # tab delimited otu table with taxonomy
    │   │   └── README.html                     # results report when annotation method is 'blast'
    │   ├── logs
    │   │   ├── cluster_sequences.log
    │   │   ├── fasttree.log
    │   │   └── uniques.log
    │   ├── OTU_aligned.fasta                   # multiple alignment file of otu seed seqs
    │   ├── OTU.fasta                           # otu seqs without taxonomy
    │   ├── OTU.tree                            # newick tree of multiple alignment
    │   └── utax
    │       ├── logs
    │       │   └── utax.log
    │       ├── OTU.biom                        # tax annotated biom (no metadata, no normalization)
    │       ├── OTU_tax.fasta                   # otu seqs with tax in FASTA header
    │       ├── OTU.txt                         # tab delimited otu table with taxonomy
    │       ├── README.html                     # results report when annotation method is 'utax'
    │       └── utax_hits.txt                   # raw UTAX hits per OTU seed seq
    ├── demux
    │   ├── *.fastq.count
    │   └── *.fastq
    ├── logs
    │   ├── quality_filtering_stats.txt
    │   └── *.count
    ├── merged_?.fasta                          # error corrected FASTA prior to clustering into OTU seqs
    ├── merged.fastq                            # all sample reads merged into single file with updated headers
    └── quality_filter
        └── *.fastq                             # files that should have been cleaned up!
```
