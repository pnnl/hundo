import bisect
from collections import defaultdict, deque

BLAST6 = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


class BlastHits(object):
    def __init__(self, names=None, max_hits=10, top_fraction=None):
        """Class that represents BLAST hits for a single target sequence. Hits are added to queues
        for bitscore and ID and ordered by increasing bitscore.
        Args:
            names (Optional[list]): when initiated with a name list; :func:`best_hit` and
                :func:`add` will no longer operate as intended
            max_hits (int): maximum number of hits to consider for this :class:`BlastHits` group
            top_fraction (float): fraction cutoff from best bitscore, e.g. 0.3 will filter out 699 when best bitscore is 1000
        Notes:
            max_hits and top_fraction work in conjunction of one another
        """
        if names is None:
            # increasing bitscore sorted
            self.names = deque()
            self.percent_ids = deque()
            self.bitscores = deque()
        else:
            self.names = names
        self.max_hits = max_hits
        self.top_fraction = top_fraction

    def __repr__(self):
        return "{cls}[{tax}]".format(cls=self.__class__.__name__, tax=self.names)

    def add(self, name, percent_id, bitscore):
        """Add entry to this :class:`BlastHits` group.
        Args:
            name (str): hit identifier
            bitscore (str): bitscore for hit
        """
        percent_id = float(percent_id)
        bitscore = float(bitscore)
        if self.top_fraction and self.bitscores:
            # the filter
            if bitscore < (self.bitscores[-1] * self.top_fraction):
                bitscore = None
            # new best
            elif bitscore > self.bitscores[-1]:
                score = self.bitscores[0]
                while score < bitscore * self.top_fraction:
                    self.names.popleft()
                    self.percent_ids.popleft()
                    self.bitscores.popleft()
                    score = self.bitscores[0]
        if bitscore:
            # insert into sorted list
            idx = bisect.bisect_left(self.bitscores, bitscore)
            self.bitscores.insert(idx, bitscore)
            self.percent_ids.insert(idx, percent_id)
            self.names.insert(idx, name)

            if len(self.names) > self.max_hits:
                # remove lowest bitscore
                self.names.popleft()
                self.percent_ids.popleft()
                self.bitscores.popleft()

    def best_hit(self):
        """Returns the hit ID of the best scoring alignment."""
        return self.names[-1]

    def majority(self):
        """Returns the hit ID of the best scoring hit ID that is repeated or the best hit when
        no items are repeated.
        """
        # no repeated names
        if len(self.names) == len(set(self.names)):
            return self.best_hit()
        else:
            # count each taxonomy, grab top taxonomy
            most_common = Counter(self.names).most_common(1)[0][0]
            # need to flip to grab best bitscore
            names_reversed = self.names.copy()
            names_reversed.reverse()
            # left most index match
            idx = names_reversed.index(most_common)
            return names_reversed[idx]


def parse_blasthits(blasthits, min_score=155, top_fraction=0.98):
    hsps = defaultdict(lambda: BlastHits(top_fraction=top_fraction))
    for hsp in blasthits:
        toks = dict(zip(BLAST6, hsp.strip().split("\t")))
        if float(toks["bitscore"]) < min_score:
            continue
        hsps[toks["qseqid"]].add(
            toks["sseqid"], float(toks["pident"]) / 100, toks["bitscore"]
        )
    return hsps


class VsearchHits(object):
    def __init__(self, names=None, max_hits=10, top_fraction=None):
        """Class that represents vsearch hits for a single target sequence. Hits are added to queues
        for precent id and seq ID and ordered by increasing percent ID.
        Args:
            names (Optional[list]): when initiated with a name list; :func:`best_hit` and
                :func:`add` will no longer operate as intended
            max_hits (int): maximum number of hits to consider for this :class:`VsearchHits` group
             top_fraction (float): fraction cutoff for percent identity. Anything below this will be dropped
        Notes:
            max_hits and top_fraction work in conjunction of one another
        """
        if names is None:
            self.names = deque()
            self.percent_ids = deque()
        else:
            self.names = names
        self.max_hits = max_hits
        self.top_fraction = top_fraction

    def __repr__(self):
        return "{cls}[{tax}]".format(cls=self.__class__.__name__, tax=self.names)

    def add(self, name, percent_id):
        """Add entry to this :class:`VsearchHits` group.
        Args:
            name (str): hit identifier
            percent id (str): percent id for hit
        """
        percent_id = float(percent_id)

        # filter instead based only on percent ID rather than bitscore
        if self.top_fraction and self.percent_ids:
            if percent_id < (self.percent_ids[-1] * self.top_fraction):
                percent_id = None
            elif percent_id >= self.percent_ids[-1]:
                score = self.percent_ids[0]
                while score < percent_id * self.top_fraction:
                    self.names.popleft()
                    self.percent_ids.popleft()
                    score = self.percent_ids[0]
        if percent_id:
            # insert into sorted list
            idx = bisect.bisect_left(self.percent_ids, percent_id)
            self.percent_ids.insert(idx, percent_id)
            self.names.insert(idx, name)
            if len(self.names) > self.max_hits:
                self.names.popleft()
                self.percent_ids.popleft()
                # self.bitscores.popleft()

    def best_hit(self):
        """Returns the hit ID of the best scoring alignment."""
        # the last one in the list is the best
        return self.names[-1]

    def majority(self):
        """Returns the hit ID of the best scoring hit ID that is repeated or the best hit when
        no items are repeated.
        """
        # no repeated names
        if len(self.names) == len(set(self.names)):
            # if the length of names is equal to the number of unique names in
            # the list, then return the best hit
            return self.best_hit()
        else:
            # count each taxonomy, grab top taxonomy
            most_common = Counter(self.names).most_common(1)[0][0]
            names_reversed = self.names.copy()
            names_reversed.reverse()
            # left most index match
            idx = names_reversed.index(most_common)
            return names_reversed[idx]


def parse_vsearchhits(file, min_pid=.85, top_fraction=0.98):
    hsps = defaultdict(lambda: VsearchHits(top_fraction=top_fraction))
    for hsp in file:
        toks = dict(zip(BLAST6, hsp.strip().split("\t")))
        if (float(toks["pident"]) / 100) < min_pid:
            continue
        hsps[toks["qseqid"]].add(toks["sseqid"], float(toks["pident"]) / 100)
    return hsps
