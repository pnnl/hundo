#!/usr/bin/env python
# coding=utf-8
"""
copy source directory files to destination directory. where a duplicate exists, if the source has
greater than n reads, concatenate the files. if exclude, files in source are not copied to
destination directory.
"""

import argparse
import logging
import os
import shutil
import subprocess


def read_count(fastq):
    """Very fast read counter (assuming the FASTQ is formatted properly).

    Args:
        fastq (str): fastq file path

    Returns:
        int
    """
    return int(subprocess.check_output("awk '{n++}END{print n/4}' " + fastq, shell=True).decode())


def main(src, dst, exclude=[], min_count=5000):
    dst_files = os.listdir(dst)
    for src_file in os.listdir(src):
        if src_file in exclude:
            logging.info("%s has been excluded" % src_file)
            continue
        if read_count(os.path.join(src, src_file)) < min_count:
            logging.info("excluding %s due to read count" % src_file)
            continue
        if src_file in dst_files:
            logging.info("Joining source and destination for %s" % src_file)
            subprocess.check_call("cat %s >> %s" % (os.path.join(src, src_file), os.path.join(dst, src_file)), shell=True)
        else:
            logging.info("Copying over source of %s" % src_file)
            shutil.copy(os.path.join(src, src_file), os.path.join(dst, src_file))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("src")
    p.add_argument("dst")
    p.add_argument("--exclude", action="append", default=[], help="files to explicitly exclude")
    p.add_argument("--min-count", type=int, default=5000, help="ignore if source file has fewer than INT reads")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s")
    main(os.path.abspath(os.path.expanduser(args.src)), os.path.abspath(os.path.expanduser(args.dst)), args.exclude, args.min_count)
