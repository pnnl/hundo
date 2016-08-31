#!/usr/bin/env python
# coding=utf-8
"""
Convert utax fasta to blast fasta and taxonomy file.
"""
import argparse
import gzip
import os


gzopen = lambda f: gzip.open(f) if f.endswith(".gz") else open(f)


def convert_header(fasta_header):
    """
    >>> convert_header(">FJ362054|SH204600.07FU;tax=d:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Marasmiaceae,g:Moniliophthora;")
    ['FJ362054|SH204600.07FU', 'k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Agaricales; f__Marasmiaceae; g__Moniliophthora; s__?']
    >>> convert_header(">DQ663396|SH221479.07FU;tax=d:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Cortinariaceae,g:Cortinarius,s:Cortinarius_eufulmineus_SH221479.07FU;")
    ['DQ663396|SH221479.07FU', 'k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Agaricales; f__Cortinariaceae; g__Cortinarius; s__Cortinarius_eufulmineus_SH221479.07FU']
    """
    assert ";tax=d:" in fasta_header
    header, _, taxonomy = fasta_header.partition(";tax=")
    taxonomy = taxonomy.strip(";").replace("d:", "k:").replace('"', '').split(",")
    tax_dict = dict(x.split(":") for x in taxonomy)
    return [header[1:], "; ".join(["%s__%s" % (idx, tax_dict.get(idx, "?")) for idx in "kpcofgs"])]


def main(utax, blast, tax):
    with gzopen(utax) as fasta_file, open(blast, "w") as blast_file, open(tax, "w") as tax_file:
        headers = set()
        for line in fasta_file:
            line = line.strip()
            # sequence line
            if not line.startswith(">"):
                print(line, file=blast_file)
            # header line
            else:
                header, taxonomy = convert_header(line)
                assert header not in headers, header
                headers.add(header)
                print(">%s" % header, file=blast_file)
                print(header, taxonomy, sep="\t", file=tax_file)


if __name__ == "__main__":

    def _file_exists(parser, arg):
        if not os.path.exists(arg):
            parser.error("The file %s does not exist" % arg)
        if not os.path.isfile(arg):
            parser.error("Expected file, not folder (%s)" % arg)
        return arg

    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("utax", type=lambda x: _file_exists(p, x), help="UTAX formatted FASTA file")
    p.add_argument("blast", help="FASTA with blast formatted headers")
    p.add_argument("tax", help="tab delimited header<tab>taxonomy file")
    args = p.parse_args()
    main(args.utax, args.blast, args.tax)
