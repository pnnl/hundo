#!/usr/bin/env python
# coding=utf-8
"""
DNA 16S Sequence ID --> DNA 16S Barcode --> DNA 16S Job ID

sample\tbarcode\tgroup

from the sequence and sample id, generate a list of output files that we should expect
generate a map file and run demultiplexing on any run that is missing a particular output file
"""

import argparse
import logging
import numpy as np
import pandas as pd
import os
import sys
from subprocess import check_call


EXPECTED_COLS = ['DNA 16S Sequence ID', 'DNA 16S Barcode', 'DNA 16S Job ID',
                 'DNA ITS Sequence ID', 'DNA ITS Barcode', 'DNA ITS Job ID',
                 'cDNA 16S Sequence ID', 'cDNA 16S Barcode', 'cDNA 16S Job ID',
                 'cDNA ITS Sequence ID', 'cDNA ITS Barcode', 'cDNA ITS Job ID']
DEMUXSH = "/people/brow015/mint/otu-16s/scripts/demux.sh"

def valid_columns(cols, df):
    logging.debug("All columns: %s" % df.columns)
    for col in cols:
        if col not in df.columns:
            logging.critical("Missing: %s" % col)
            return False
    return True


def transform_cols(df, cols):
    def _transform(x):
        if isinstance(x, int):
            return "{}".format(x)
        elif isinstance(x, float):
            return "{:.0f}".format(x)
        elif isinstance(x, str):
            return x.strip()
        else:
            return x
    df = df.applymap(_transform)
    df = df.replace("nan", "")
    return df


def map_dicts(df, amplicon):
    runinfo = dict()

    # runinfo['150114-16s']['SP0514-FP0469-D1'] = "ACTG"
    for job_id in df.columns:
        if not "Job ID" in job_id: continue
        for key, groupdf in df.groupby(job_id):
            # not sequenced yet
            if not key: continue

            if "16s" in job_id.lower() and "16s" in amplicon.lower():
                jid = "%s-16s" % key
                runinfo[jid] = dict()
                for i, r in groupdf.iterrows():
                    if r['DNA 16S Sequence ID'] and r['DNA 16S Barcode']:
                        runinfo[jid][r['DNA 16S Sequence ID']] = r['DNA 16S Barcode']
                    if r['cDNA 16S Sequence ID'] and r['cDNA 16S Barcode']:
                        runinfo[jid][r['cDNA 16S Sequence ID']] = r['cDNA 16S Barcode']

            elif "its" in job_id.lower() and "its" in amplicon.lower():
                jid = "%s-ITS" % key
                runinfo[jid] = dict()
                for i, r in groupdf.iterrows():
                    if r['DNA ITS Sequence ID'] and r['DNA ITS Barcode']:
                        runinfo[jid][r['DNA ITS Sequence ID']] = r['DNA ITS Barcode']
                    if r['cDNA ITS Sequence ID'] and r['cDNA ITS Barcode']:
                        runinfo[jid][r['cDNA ITS Sequence ID']] = r['cDNA ITS Barcode']

            else:
                pass

    return runinfo


def files_from_dict(dct, eid, output_dir):
    translation = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    reverse_complement = lambda s: str.translate(s[::-1], translation)

    mapping_files = {}
    for runid, sample_meta in dct.items():
        # barcodes=results/$eid/demux/${runid}_barcodes.txt
        output_file = os.path.join(output_dir, eid, "demux", "%s_barcodes.txt" % runid)
        # make the full directory tree at this point
        os.makedirs("%s/%s/%s" % (output_dir, eid, "demux"), exist_ok=True)

        all_exist = True
        for sample_name, sample_barcode in sample_meta.items():
            if not os.path.exists(os.path.join(output_dir, eid, "%s_R1.fastq" % sample_name)) or not os.path.exists(os.path.join(output_dir, eid, "%s_R2.fastq" % sample_name)):
                all_exist = False

        if not all_exist:
            # write a new mapping file
            with open(output_file, 'w') as fh:
                for sample_name, sample_barcode in sample_meta.items():
                    # writes original and reverse-complement groups for fastq-multx
                    print(sample_name, sample_barcode, "original", sep="\t", file=fh)
                    print(sample_name, reverse_complement(sample_barcode), "reverse-complement", sep="\t", file=fh)

            # only include mapping files for runs that need demuxing
            mapping_files[runid] = output_file
            logging.info("Wrote new barcodes file for run %s containing %d samples for experiment %s" % (runid, len(sample_meta), eid))
        else:
            logging.info("Skipping run %s since all FASTQ files are already present")

    return mapping_files


def launch_jobs(files, eid):
    for runid, map_file in files.items():
        logging.info("Submitting demultiplexing job for run %s" % runid)
        cmd = "sbatch {script} {runid} {eid}".format(script=DEMUXSH, runid=runid, eid=eid)
        check_call(cmd, shell=True)
448

def main(metadata, eid, amplicon, barcode_length, output_dir):
    metadf = pd.read_excel(metadata)
    metadf.columns = [i.strip() for i in metadf.columns]
    if not valid_columns(EXPECTED_COLS, metadf):
        logging.critical("The expected columns for DNA, cDNA, and ITS are not all present in %s" % metadata)
        sys.exit(1)
    metadf = transform_cols(metadf, EXPECTED_COLS)
    run_mapping = map_dicts(metadf, amplicon)
    map_files = files_from_dict(run_mapping, eid, output_dir)
    launch_jobs(map_files, eid)


if __name__ == "__main__":

    def _file_exists(parser, arg):
        if not os.path.exists(arg):
            parser.error("The file %s does not exist" % arg)
        if not os.path.isfile(arg):
            parser.error("Expected file, not folder (%s)" % arg)
        return arg

    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("metadata", type=lambda x: _file_exists(p, x), help="Sample metadata in XLSX format with headers defined in the description.")
    p.add_argument("eid", metavar="ID", help="Experiment ID. This is prepended to the map file and is used when processing this run.")
    p.add_argument("amplicon", metavar="AMPLICON", choices=['16s', 'ITS'], help="Amplicon type")
    p.add_argument("--barcode-length", default=12, action="append", help="Expected barcode length. Can be specified more than once.")
    p.add_argument("--output-dir", default="results", help="Top level output directory, e.g 'results' becomes 'results/<run id>/<experiment id>_map.txt'.")
    p.add_argument("--overwrite", action="store_true", help="Overwrite any existing mapping file and re-process sequencing run.")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s")
    main(args.metadata, args.eid, args.amplicon, args.barcode_length, args.output_dir)
