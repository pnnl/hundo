#!/usr/bin/env python

# import libraries
import sys # for argument vector
import optparse # to parse arguments
from src.LCAStar import *
import re

# describe what the script does
what_i_do = "A simple script for printing files"
usage = "myscript.py -i <input_file> -o output_file"

# initialize the parser
parser = optparse.OptionParser(usage = usage, epilog = what_i_do)
parser.add_option("-i", "--input_file", dest="input_file", default=None,
                   help='file to print out to the screen [Required]')
parser.add_option("-t", "--taxa_map", dest="taxa_map", default=None,
                   help='file to print out to the screen [Required]')
parser.add_option("-m", "--megan_map", dest="megan_map", default=None,
                   help='MEGAN file with preferred names [Required]')
parser.add_option("-l", "--level", dest="level", default=None,
                   help='file to print out to the screen [Required]')
parser.add_option("-o", "--output_file", dest="output_file", default=None,
                   help='output file to write to')

taxa_pattern = re.compile("\(.*\)")

def checkArguments(opts):
    # checks arguents and returns True if there is a problem
    opts_dict = opts.__dict__
    for key, value in opts_dict.items():
        if not value:
            print "Error:", key, "not defined"
            print usage
            exit()

def translate_to_prefered_name(id, ncbi_megan_map, lcastar):
    # This maps an NCBI Taxonomy Database ID to the prefered MEGAN Taxonomy name and
    # reports the default name on the NCBI Taxonomy Database otherwise.
    # id: id to be translated
    # ncbi_megan_map: dictionary of NCBI ID to MEGAN name
    # lcastar: instance of the LCAStar class
    id_str = str(id)
    if id_str in ncbi_megan_map:
        return ncbi_megan_map[id_str] + " (" + id_str + ")"
    else:
        res = lcastar.translateIdToName(id_str)
        if res:
            return res + " (" + id_str + ")"
        else:
            return "Unknown (" + id_str + ")"

def createTaxaLine(lin, ncbi_megan_map, lcastar):
    res = ";".join( [ translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in lin ] )
    return res

def getTaxonomyFromLine(line):
    fields = line.split("\t")

    if ((len(fields) > 7) and (fields[0] != "")):
        search_res = taxa_pattern.search(fields[8])
        if search_res:
            res = fields[8].split("(")
            res = res[1].split(")")
            res = res[0]
            return res
        else:
            return None
    else:
        return None

def getLevelCounts(func_taxa_table_file, lcastar, megan_map, level):

    level_counts = {}

    with open(func_taxa_table_file, "r") as fh:
        header = fh.readline()
        for line in fh:
            taxa_id = getTaxonomyFromLine(line)
            if taxa_id:
                lin = lcastar.get_lineage(taxa_id)
                tmp_lin = lin[::-1]
                if len(tmp_lin) > level:
                    tmp_lin = tmp_lin[0:level]
                dec_taxa_str = createTaxaLine(tmp_lin, megan_map, lcastar)
                if dec_taxa_str not in level_counts:
                    level_counts[dec_taxa_str] = 0
                level_counts[dec_taxa_str] += 1
            else:
                continue

    return level_counts

# the main function of the script
def main():
    # load and extract object
    (opts, args) = parser.parse_args()
    # checkArguments(opts)

    ## create preferred mapping
    ncbi_megan_map = {} # hash map from given taxonomy to preferred one used by megan
    with open(opts.megan_map, 'r') as meganfile:
        for line in meganfile:
             fields = line.split("\t")
             fields = map(str.strip, fields)
             ncbi_megan_map[fields[0]] = fields[1]

    lcastar = LCAStar(opts.taxa_map)
    level_counts = getLevelCounts(opts.input_file, lcastar, ncbi_megan_map, int(opts.level))

    for level, counts in level_counts.items():
        print level + "\t" + str(counts)

    # print createTaxaLine(lin[::-1], ncbi_megan_map, lcastar)










if __name__ == "__main__":
    main()
