# Install


Install dependencies using Bioconda repository:

```
conda install -c bioconda snakemake biom-format fasttree bbmap blast clustalo
```

Unfortunately, [USEARCH](http://www.drive5.com/usearch/download.html) is currently a dependency of this protocol.

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
