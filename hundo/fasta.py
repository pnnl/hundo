from itertools import groupby


def read_fasta(fh):
    """Fasta iterator.

    Accepts file handle of .fasta and yields name and sequence.

    Args:
        fh (file): Open file handle of .fasta file

    Yields:
        tuple: name, sequence

    >>> import os
    >>> from itertools import groupby
    >>> f = open("test.fasta", 'w')
    >>> f.write("@seq1\nACTG")
    >>> f.close()
    >>> f = open("test.fastq")
    >>> for name, seq in read_fastq(f):
            assert name == "seq1"
            assert seq == "ACTG"
    >>> f.close()
    >>> os.remove("test.fasta")
    """
    for header, group in groupby(fh, lambda line: line[0] == ">"):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = "".join(line.strip() for line in group)
            yield name, seq


def format_fasta_record(name, seq, wrap=80):
    """Fasta __str__ method.

    Convert fasta name and sequence into wrapped fasta format.

    Args:
        name (str): name of the record
        seq (str): sequence of the record
        wrap (int): length of sequence per line

    Yields:
        tuple: name, sequence

    >>> format_fasta_record("seq1", "ACTG")
    ">seq1\nACTG"
    """
    record = ">" + name + "\n"
    if wrap:
        for i in range(0, len(seq), wrap):
            record += seq[i : i + wrap] + "\n"
    else:
        record += seq + "\n"
    return record.strip()
