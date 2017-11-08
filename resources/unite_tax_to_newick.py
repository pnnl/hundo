#!/usr/bin/env python
# coding=utf-8
"""
"""
import click
import os
from collections import defaultdict
from hundo.fasta import read_fasta, format_fasta_record


@click.group()
@click.pass_context
def cli(obj):
    """
    """


@cli.command("tax-to-newick")
@click.argument("tax", type=click.File("r"))
@click.argument("fasta", type=click.File("r"))
@click.argument("outfasta", type=click.File("w"))
@click.argument("outmap", type=click.File("w"))
@click.argument("outtre", type=click.File("w"))
def tax_to_newick(tax, fasta, outfasta, outmap, outtre):
    """
    Tax and FASTA input files represent clusters at 99%% identity via:
    https://unite.ut.ee/sh_files/sh_mothur_release_10.10.2017.zip
    """
    def tree():
        return defaultdict(tree)

    def tree_add(t, path):
      for node in path:
        t = t[node]

    def tree_to_newick(root):
        items = []
        for k in root.keys():
            s = ''
            if len(root[k].keys()) > 0:
                sub_tree = tree_to_newick(root[k])
                if sub_tree != '':
                    s += '(' + sub_tree + ')'
            s += k
            items.append(s)
        return ','.join(items)

    t = tree()
    saved = set()
    for line in tax:
        toks = line.strip().split("\t")
        taxonomies = toks[1].strip(";").split(";")
        if not taxonomies[0] == "k__Fungi": continue
        assert(len(taxonomies) == 7)
        tree_add(t, taxonomies)
        print(toks[0], taxonomies[6], sep="\t", file=outmap)
        saved.add(toks[0])

    tree_str = tree_to_newick(t)
    print(tree_str, file=outtre)

    for name, seq in read_fasta(fasta):
        if name in saved:
            print(format_fasta_record(name, seq), file=outfasta)


if __name__ == '__main__':
    cli()
