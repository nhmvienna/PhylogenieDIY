import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--tree", dest="tree", help="Input file")
parser.add_option("--taxa", dest="taxa", help="Output file")
parser.add_option("--list", dest="list", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

Tree = open(options.tree, "rt").readline()[:-2]
Taxa = [x.split(":")[0]
        for x in Tree.replace("(", "").replace(")", "").split(",")]

ORD = options.list.split(",")
List = []
for l in open(options.taxa, "rt"):
    a = l.rstrip().split(",")
    if len(a) < 4:
        # print(l)
        continue
    Ord = a[5]
    ID = "_".join(a[2].split(" ")[:2])
    if Ord not in ORD or ID not in Taxa:
        continue
    List.append(ID)
print(",".join(["'" + x + "'" for x in List]))
