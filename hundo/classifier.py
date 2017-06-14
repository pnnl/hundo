import logging
import os
import sys
from collections import defaultdict, OrderedDict
from hundo.blast import BlastHits
from hundo.fasta import format_fasta_record, read_fasta


RANKS = {8:"s", 7:"g", 6:"f", 5:"o", 4:"c", 3:"p", 2:"k"}


class Node:
    def __init__(self, name=None, parent=None, node_id=None, depth=None):

        self.name = name
        self.parent = parent
        self.children = []
        self.node_id = node_id
        self.depth = depth
        if parent:
            self.parent.add_child(child_node=self)

    def add_child(self, child_node):
        self.children += ([child_node])

    def is_root(self):
        return not self.parent

    def get_phylogeny(self, limit=False):
        """Return all nodes until root in a list, in order of decreasing depth
        (ACC,Species,Genus..Root).
        """
        phyl = [self]
        while not phyl[-1].is_root():
            phyl.append(phyl[-1].parent)
        if limit:
            while phyl[0].depth > limit:
                phyl = phyl[1:]
        return phyl

    def to_tax_string(self, root=True, limit=False, sep=","):
        tax_defs = {}
        phyl = self.get_phylogeny(limit)

        if not root:
            phyl = phyl[:-1]

        for p in phyl:
            try:
                depth_str = RANKS[p.depth]
            except KeyError:
                continue
            if depth_str in tax_defs:
                continue
            tax_defs[depth_str] = p.name

        complete_taxonomy = []
        for idx in "kpcofgs":
            complete_taxonomy.append("%s__%s" % (idx, tax_defs.get(idx, "?")))
        return sep.join(complete_taxonomy)

    def highest_rank(self):
        """Return the highest rank to which this node is classified, assuming
        intermediate parents without rank.
        """
        max_rank = 0
        for node in self.get_phylogeny():
            if node.depth and node.depth > max_rank:
                max_rank = node.depth
        return max_rank

    def _get_children_by_rank(self, rank_children, rank):
        """Recursively add to rank_children if the depth is right, too high
        or if it has no more kids and the depth is. Otherwise iterate with
        all all kids.
        """
        # If the rank is right, or no more children and the rank remains
        # unknown and not found (once)

        if ((self.depth and self.depth >= rank) or
            (((not self.depth) or
              self.depth == Tree.NORANK) and
             (not self.children) and
             self.highest_rank() == rank - 1)):
                if self.parent:
                    p = self.parent
                    self.parentPrintName = p.name
                    while p.parent and (not p.depth or p.depth == Tree.NORANK):
                        p = p.parent
                        self.parentPrintName = (p.name + " : " +
                                                self.parentPrintName)
                else:
                    self.parentPrintName = ""
                rank_children.append(self)
        else:
            # not matching
            for child in self.children:
                rank_children = child._get_children_by_rank(rank_children, rank)
        return rank_children

    def _prune_unassigned(self, to_remove):
        #Exhaustive recursion child by child stopping when node removed
        if not self.is_root():
            to_remove.append(self)
        else:
            for child in self.children:
                child._prune_unassigned(to_remove)

    def __repr__(self):
        if self.is_root():
            return "Root Node: %s" % (self.name)
        else:
            return "Node: %s, parent: %s" % (self.name, self.parent.name)


class Tree:

    NORANK = 0
    META = 1
    DOMAIN = 2
    PHYLUM = 3
    CLASS = 4
    ORDER = 5
    FAMILY = 6
    GENUS = 7
    SPECIES = 8
    SUBSPECIES = 9
    SUB2 = 10

    depths = {META: "base", DOMAIN: "domain", PHYLUM: "phylum",
              CLASS: "class", ORDER: "order", FAMILY: "family",
              GENUS: "genus", SPECIES: "species"}

    # domain 10000 to account for map coming from SILVA
    map_codes = {NORANK: 0, META: 0, DOMAIN: [0, 10000], PHYLUM: 2, CLASS: 3,
                 ORDER: 4, FAMILY: 5, GENUS: 98, SPECIES: 100, SUBSPECIES: 101}

    def __init__(self, root_node=None, name=None):
        if not root_node:
            root_node = Node(name="root", node_id=1)
        self.root = root_node
        self.name = name
        self.nodes = {self.root.name: self.root}
        self.nodesNoBrackets = {}
        self.idCount = 10000000

    def add_node(self, node):
        self.nodes[node.name] = node

    def get_node(self, name):
        return self.nodes[name] if name in self.nodes else None

    def new_id(self):
        self.idCount += 1
        return self.idCount

    def move_node(self, node, new_parent):
        if isinstance(node, str):
            node = self.get_node(node)
        i = node.parent.children.index(node)
        del node.parent.children[i]
        if isinstance(new_parent, str):
            new_parent = self.get_node(new_parent)
        node.parent = new_parent
        new_parent.add_child(node)

    def delete_node(self, node):
        """Removes node from tree structure. Still retains the disconnected
        node in self.nodes dict.
        """
        if isinstance(node, str):
            nodename = node
            node = self.get_node(nodename)
            if not node:
                logging.warning("Node '%s' not found" % nodename)
                return
        i = node.parent.children.index(node)
        del node.parent.children[i]
        del self.nodes[node.name]

    def prune_unassigned(self):
        to_remove = []
        self.root._prune_unassigned(to_remove)
        for node in to_remove:
            self.delete_node(node)

    def read_tree(self, map_file, tre_file):
        with open(map_file) as tax_map:
            code_map = {}

            for k, v in Tree.map_codes.items():
                if isinstance(v, list):
                    for i in v:
                        code_map[i] = k
                else:
                    code_map[v] = k

            self.node_ids = {}
            self.node_ids[self.root.node_id] = self.root.name
            self.ranks = {}
            for line in tax_map:
                parts = line.split("\t")
                node_id = int(parts[0])
                self.node_ids[node_id] = parts[1]
                try:
                    self.ranks[node_id] = code_map[int(parts[-1])]
                except:
                    sys.stderr.write("Error: Depth key not found: %s\n" %
                                     int(parts[-1]))
                    self.ranks[node_id] = None

        tree = ""
        with open(tre_file) as tre:
            for line in tre:
                tree += line.strip()

        parent = self.root
        idx = ""
        for i in range(1, len(tree) + 1):
            char = line[-i]
            if char == ";":
                pass
            elif char == "," or char == "(" or char == ")":
                if not idx == "":
                    node_id = int(idx)
                    idx = ""
                    if node_id == 1:
                        n = self.root
                    else:
                        d = self.ranks[node_id]
                        # Fix for MEGAN compability forcing all ranks < phylum
                        # to be 2
                        if d == 0:
                            pl = len(parent.get_phylogeny())
                            if pl < Tree.PHYLUM:
                                d = pl
                        #--Fix--
                        n = Node(name=self.node_ids[node_id], depth=d,
                                 node_id=node_id, parent=parent)
                        self.add_node(n)

                if char == "(":
                    parent = parent.parent
                elif char == ")":
                    parent = n
            else:
                idx = char + idx

    def get_children_by_rank(self, rank):
        return self.root._get_children_by_rank([], rank)


class Classifier(Tree):

    def __init__(self, fasta_file):
        Tree.__init__(self, Node(name="root", node_id=1), name="ref")
        self.noHits = Node("No hits", parent=self.root)
        self.add_node(self.noHits)
        self.seqs = OrderedDict()
        self.read_node_assignments = {}

        self.taxonomy_limits = {Tree.SPECIES: .99,
                                Tree.GENUS: .97,
                                Tree.FAMILY: .95,
                                Tree.ORDER: .90,
                                Tree.CLASS: .85,
                                Tree.PHYLUM: .80}

        for name, seq in read_fasta(fasta_file):
            self.seqs[name] = seq

    def parse_blast(self, blast_hits, top_fraction=0.95, min_bitscore=100):
        hsps = defaultdict(lambda: BlastHits(top_fraction=top_fraction))
        with open(blast_hits) as blast_hits_fh:
            for hsp in blast_hits_fh:
                toks = dict(zip(BLAST6, hsp.strip().split("\t")))
                if float(toks["bitscore"]) < min_bitscore: continue
                hsps[toks["qseqid"]].add(toks["sseqid"],
                                         float(toks["pident"]),
                                         toks["bitscore"])
        for otu_name, hits in hsps.items():
            target = hits.best_hit()
            sequence = self.seqs[otu_name]
            node = self.get_node(target)
            if not node:
                logging.critical("%s was not found in your MAP" % target)
                sys.exit(1)
            parents = node.get_phylogeny()[1:]
            for hit_name in hits.names:
                # hit_name = hit[0]
                tree_node = self.get_node(hit_name)
                if not tree_node:
                    logging.warning("Node %s was not found" % hit[0])
                else:
                    tree_parent = tree_node.parent
                    while tree_parent not in parents:
                        tree_parent = tree_parent.parent
                    parents = parents[parents.index(tree_parent):]
            hsp_sim = hits.percent_ids[-1]

            max_rank_limit = Tree.SPECIES
            max_rank = max_rank_limit
            d = max_rank_limit
            ranks = list(self.taxonomy_limits.keys())
            ranks.sort(reverse=True)
            for rank in ranks:
                if hsp_sim < self.taxonomy_limits[rank]:
                    max_rank = rank - 1
                else:
                    break

            while (max_rank < Tree.SPECIES and max_rank < parents[0].highest_rank()):
                d = min(parents[0].highest_rank(), max_rank_limit)
                parents = parents[1:]

            if d < max_rank_limit:
                novelName = ("Unknown %s %s" % (parents[0].name, Tree.depths[d]))
                nn = self.get_node(novelName)
                if nn:
                    novelNode = nn
                else:
                    depth = parents[0].highest_rank() + 1
                    novelNode = Node(novelName, parent=parents[0], depth=depth)
                    self.add_node(novelNode)
                parents = [novelNode] + parents

            self.read_node_assignments[otu_name] = parents[0]

        for name, seq in self.seqs.items():
            if name in self.read_node_assignments: continue
            self.read_node_assignments[name] = self.noHits

    def print_sequences(self, out_file):
        with open(out_file, "w") as ofh:
            for name, seq in self.seqs.items():
                node = self.read_node_assignments[name]
                full_name = "{name};tax={taxonomy}".format(name=name, taxonomy=node.to_tax_string(root=False))
                print(format_fasta_record(full_name, seq), file=ofh)

    def print_table(self, out_file):
        with open(out_file, "w") as ofh:
            for name in self.seqs:
                node = self.read_node_assignments[name]
                taxonomy = node.to_tax_string(root=False)
                print(name, taxonomy, sep="\t", file=ofh)
