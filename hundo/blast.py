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

                # print('bitscore',bitscore)
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
            print("most common")
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


##vsearch.py


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
            # increasing bitscore sorted
            # below commands just take the arguements passed to the classifier and make them empty deque objects.
            self.names = deque()
            self.percent_ids = deque()
            # self.bitscores = deque()
        else:
            self.names = names
        self.max_hits = max_hits  # always 10
        self.top_fraction = top_fraction  # defined by user

    def __repr__(self):
        return "{cls}[{tax}]".format(
            cls=self.__class__.__name__, tax=self.names
        )  # this makes it human readable if you just call teh class it will tell you what it is

    def add(self, name, percent_id):
        """Add entry to this :class:`VsearchHits` group.
        Args:
            name (str): hit identifier
            percent id (str): percent id for hit
        """
        percent_id = float(percent_id)
        # bitscore = float(bitscore)

        if self.top_fraction and self.percent_ids:
            # the filter

            if percent_id < (
                self.percent_ids[-1] * self.top_fraction
            ):  # if the last item in the bitscores * the fraction decided is less than
                # the bitscore provided, then apply a None value

                percent_id = None

            # new best

            elif percent_id >= self.percent_ids[-1]:
                # if bitscore provided is larger than the last element in the list, then
                # set score to be the first item of the bitscores list

                score = self.percent_ids[0]

                while (
                    score < percent_id * self.top_fraction
                ):  # as long as the score is not larger than the bitscore * the fraction determined
                    # then the following commands will remove the first entry in that desque object
                    self.names.popleft()
                    self.percent_ids.popleft()
                    # self.bitscores.popleft()
                    score = self.percent_ids[
                        0
                    ]  # once you have removed that first entry, set score to be equal to the "new"
                    # first entry(basically the second)

        if percent_id:

            # insert into sorted list
            idx = bisect.bisect_left(self.percent_ids, percent_id)
            # insert x(the bitscore) in a (the deque of bitscores)
            # this will save the position of where that bitscore should go to maintain the sorted order. clever.

            # the below inserts the wanted information into the specific desque objects. The first arguement is the position
            # provided above, along with the information to be included.
            # self.bitscores.insert(idx, bitscore)
            self.percent_ids.insert(idx, percent_id)
            self.names.insert(idx, name)

            if (
                len(self.names) > self.max_hits
            ):  # if there are more names than there are allowed maximum hits, remove the
                # remove lowest bitscore object which would be the first one in the list
                self.names.popleft()
                self.percent_ids.popleft()
                # self.bitscores.popleft()

    def best_hit(self):
        """Returns the hit ID of the best scoring alignment."""
        return self.names[-1]  # the last one in the list is the best

    def majority(self):
        """Returns the hit ID of the best scoring hit ID that is repeated or the best hit when
        no items are repeated.
        """
        # no repeated names
        if len(self.names) == len(
            set(self.names)
        ):  # if the length of names is equal to the number of unique names in
            # the list, then return the best hit
            return self.best_hit()
        else:

            # count each taxonomy, grab top taxonomy
            most_common = Counter(self.names).most_common(1)[0][
                0
            ]  # so of the names, take the most common,
            # this is returned as a list of tuples, so grab the first item of the list, and then the first item of the tuple
            # need to flip to grab best bitscore
            names_reversed = self.names.copy()  # copy the names
            names_reversed.reverse()  # reverse their order
            # left most index match
            idx = names_reversed.index(
                most_common
            )  # grab the position of the most common when the order is reversed
            return names_reversed[idx]  # return the position of the reversed order


def parse_vsearchhits(file, min_pid=.85, top_fraction=0.98):
    hsps = defaultdict(lambda: VsearchHits(top_fraction=top_fraction))
    # so this makes a dictionary that will generate
    # the default value to be the top blast hit
    # print(hits)
    for hsp in file:
        # print(hsp)

        toks = dict(zip(BLAST6, hsp.strip().split("\t")))
        # print(toks)#so for each value in the blasthits object,
        # create a dictionary called toks where the keys are the string of what the blast6 headers are and then the values
        # are the different items in the blasthits object
        # if toks["qseqid"] == "OTU_126":
        #     print("first:", min_pid, toks)
        if (float(toks["pident"]) / 100) < min_pid:
            continue  # if the bitscore is lower than the min socore we set of 155, then
        # take the value in toks that has the qseqid and add the ssequid, the percent identity and the bitscore
        hsps[toks["qseqid"]].add(toks["sseqid"], float(toks["pident"]) / 100)

    return hsps
