import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--threshold", dest="th", help="Output file")


(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    """ import data either from a gzipped or or uncrompessed file or from STDIN"""
    import gzip

    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


MATRIX = d(list)
Names = []

for l in load_data(options.IN):
    if l.startswith(">"):
        Names.append(l.rstrip()[1:])
        C = 0
        continue
    for i in range(len(l.rstrip())):
        C += 1
        MATRIX[C].append(l.rstrip()[i])

for k, v in list(sorted(MATRIX.items())):
    # print(k, v.count("-") / len(v))
    if v.count("-") / len(v) > float(options.th):
        # print(len(MATRIX))
        del MATRIX[k]

Print = d(list)

for k, v in list(sorted(MATRIX.items())):
    for i in range(len(Names)):
        Print[Names[i]].append(v[i])

for k, v in sorted(Print.items()):
    print(">" + k + "\n" + "".join(v))
