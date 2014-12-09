import os, sys, re
from tools import *

fout = open("summary.ancestors.txt", "w")
for d in DIRS:
    for m in get_phylo_modelnames(con):
        for i in ap.params["ingroup"]:
            [start,stop] = get_boundary_sites(  get_phylippath(d), ap.params["seedtaxa"][i]  )
            print d, start, stop
            fin = open(d + "/asr." + get_runid(d,m) + "/" + i + ".txt", "r")
            lines = fin.readlines()
            seq = lines[ lines.__len__()-2 ]
            seq = seq[start-1:stop]
            seq = re.sub("-", "", seq)
            fout.write(d + "\t" + m + "\t" + i + "\t" + seq + "\n")
fout.close()
