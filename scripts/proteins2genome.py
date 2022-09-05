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
parser.add_option("--Tax", dest="TA", help="Output file")
parser.add_option("--TaxList", dest="TL", help="Output file")
parser.add_option("--FreqTH", dest="FT", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

TaxHash = d(str)


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


for l in open(options.TL, "rt"):
    a = l.rstrip().split(",")
    if len(a) < 3:
        continue
    if options.TA not in l:
        continue
    ID = a[2].replace(" ", "_")
    TaxHash[ID] = "_".join([a[0], ID])

Genes = d(int)
SeqHash = d(lambda: d(list))
for l in load_data(options.IN):
    if l.startswith(">"):
        ID = l.split("[")[-1][:-2].replace(" ", "_")
        a = l.rstrip().split()
        Gene = "_".join(a[1:-2])
        Genes[Gene] += 1
        continue
    SeqHash[ID][Gene].append(l.rstrip())

for Gene, Count in list(Genes.items()):
    if Count / len(SeqHash.keys()) < float(options.FT):
        del Genes[Gene]

for ID, v in sorted(SeqHash.items()):
    # print(ID)
    if ID not in TaxHash:
        continue
    Seq = []
    for Gene, Count in sorted(Genes.items()):
        Seq.extend(v[Gene])
    if Seq == []:
        continue
    print(">" + TaxHash[ID])
    print("\n".join(Seq))
