import functools
from itertools import zip_longest

from Bio import Phylo


def memoize(func):
    cache = func.cache = {}

    @functools.wraps(func)
    def memoized_func(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = func(*args, **kwargs)
        return cache[key]

    return memoized_func


class Tree(object):
    """"""

    def __init__(self, newick_tree):
        self.tree = Phylo.read(
            newick_tree, "newick", values_are_confidence=True, rooted=True
        )

    @memoize
    def find_clades(self, tx):
        """
        Save clades for all observed species assignments to speed up
        annotations.
        """
        return [i for i in self.tree.find_clades(tx)][0]

    @memoize
    def get_path(self, clade):
        cpath = self.tree.get_path(clade)
        cpath.insert(0, self.tree.root)
        return cpath

    def get_clade(self, tx, percent_id):
        if not tx:
            return None

        if percent_id < 0.80:
            cutoff = "k"
        elif percent_id >= 0.80 and percent_id < 0.85:
            cutoff = "p"
        elif percent_id >= 0.85 and percent_id < 0.90:
            cutoff = "c"
        elif percent_id >= 0.90 and percent_id < 0.95:
            cutoff = "o"
        elif percent_id >= 0.95 and percent_id < 0.97:
            cutoff = "f"
        elif percent_id >= 0.97 and percent_id < 0.99:
            cutoff = "g"
        else:
            cutoff = "s"

        c = self.find_clades(tx)
        cpath = self.get_path(c)
        for clade in cpath:
            if clade.name.startswith(cutoff):
                # cutoff due to alignment percentage
                return clade
        # target was less specific than cutoff
        return c

    @staticmethod
    def tax_str(tx):
        if isinstance(tx, str):
            tx = [tx]
        for i, t in zip_longest("kpcofgs", tx):
            if not t:
                tx.append("%s__?" % i)
            else:
                assert t.startswith("%s__" % i)
        return tx

    def lca(self, tx, percent_id):
        if isinstance(tx, str):
            tx = [tx]
        c = set([self.get_clade(i, percent_id) for i in tx])
        c = self.tree.common_ancestor(c)
        if c == self.tree.root:
            return self.tax_str([self.tree.root.name])
        else:
            l = self.get_path(c)
            return self.tax_str([c.name for c in l])
