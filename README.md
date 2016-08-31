<img src="resources/dag.png" height="500"/>

# Install

+ snakemake - `pip install snakemake`
+ biom-format - `pip install biom-format`
+ usearch - `module load usearch`
+ clustalo - `module use /people/brow015/modulefiles && module load cbb/clustalo/1.2.0`
+ fasttree - `conda install -c bioconda fasttree`
+ bbtools - `conda install -c bioconda bbmap`
+ blast


# Usage

Place demultiplexed, uncompressed reads into `results/<EXPERIMENT NAME>/demux`
with a `.fastq` file extension.

Uncompress the BLAST database and taxonomy file:

```
cd ref/silva_123
sh extract.sh
gunzip SLV_123_SSU.tax.gz
```

Build `lca`:

```
cd resources/lca_src
make
```

To run the workflow across 24 cores:

```
snakemake -j 24 --configfile resources/16s.config.yaml --config eid=<EXPERIMENT NAME>
```

Results are written to `results/<EXPERIMENT NAME>/<PERCENT ID>/`.

Methods and summary data can be found in `results/<id>/<percent_identity>/README.html`

![readme](resources/readme_example.png)
