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
parser.add_option(
    "--logical", dest="log", help="logical parameter", action="store_true"
)
parser.add_option("--param", dest="param",
                  help="numerical parameter", default=1)

(options, args) = parser.parse_args()
parser.add_option_group(group)

Tree = open(options.IN, "rt").readline()[:-2]

Taxa = [x.split(":")[0]
        for x in Tree.replace("(", "").replace(")", "").split(",")]

Rename = d(str)

for Taxon in Taxa:
    Name = "_".join(Taxon.split("_")[2:])
    Tree = Tree.replace(Taxon, Name)

print(Tree + ";")
# print(Taxa)
