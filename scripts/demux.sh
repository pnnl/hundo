#!/bin/sh
#SBATCH --account=mint
#SBATCH --partition=shared
#SBATCH --time=90
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="demux"
#SBATCH --output="demux_%A.out"
#SBATCH --error="demux_%A.err"


set -o nounset

function usage () {
    echo "Demultiplex a run with a barcodes TSV file"
    echo "Usage: $(basename $0) <RUN ID> <EXPERIMENT ID>"
    exit 1
}

if [ $# -lt 2 ]; then
    usage
fi

runid=$1
eid=$2

r1=/pic/projects/mint/hofmockel/data/${runid}/Undetermined_S0_L001_R1_001.fastq.gz
r2=/pic/projects/mint/hofmockel/data/${runid}/Undetermined_S0_L001_R2_001.fastq.gz
i1=/pic/projects/mint/hofmockel/data/${runid}/Undetermined_S0_L001_I1_001.fastq.gz

barcodes=/pic/projects/mint/otu-16s/results/$eid/demux/${runid}_barcodes.txt
stats=/pic/projects/mint/otu-16s/results/$eid/demux/${runid}_demultiplexing_stats.txt

if [[ ! -d /pic/projects/mint/otu-16s/results/$eid/demux ]]; then
    mkdir /pic/projects/mint/otu-16s/results/$eid/demux
fi

fastq-multx -m 0 -l $barcodes $i1 $r1 $r2 -o n/a \
    -o /pic/projects/mint/otu-16s/results/$eid/demux/%_R1.fastq \
    -o /pic/projects/mint/otu-16s/results/$eid/demux/%_R2.fastq > $stats

rm /pic/projects/mint/otu-16s/results/$eid/demux/unmatched_R1.fastq
rm /pic/projects/mint/otu-16s/results/$eid/demux/unmatched_R2.fastq
