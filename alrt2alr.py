#
# clean_trees.py
# 
# Victor Hanson-Smith
# victorhansonsmith@gmail.com
#
# This script will remove branch support values from a Newick-formatted tree.
# This is useful if you're using PAML, which seems to be unable to deal with
# branch labels.
#
# USAGE: 
#     python clean_trees.py FILEIN > FILEOUT
#
# . . .where FILEIN is a text file containing one or more Newick-formatted trees,
# with one tree on each line.  FILEOUT is the desired output filename.  
# After the script finished, FILEOUT will contain the support-less versions 
# of each tree in FILEIN. If FILEOUT does not exist, then it will be created de novo.
#
import math,re,sys,os
fin = open(sys.argv[1], "r")
newline = ""
for line in fin.readlines():
    #print line
    tokens = line.split(":")
    for t in tokens:
        if t.__contains__(")") and False == t.__contains__(";"):
            #print t
            ts = t.split(")")
            newline += ts[0] + ")"
            alrt = float(ts[1])
            if alrt > 1418:
                alrt = 1418.0
            alr = math.exp(alrt/2.0)
            #print alr
            #newline += alr.__str__()
            newline += "%.4e"%alr
            newline += ":"
        else:
            newline += t
            newline += ":"
    print newline
fin.close()
