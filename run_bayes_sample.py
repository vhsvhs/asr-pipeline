#
#
#
import random,sys,os
from configuration import *
from tools import *


for d in get_alignment_method_names(con):
    for m in get_phylo_modelnames(con):
        print "\n. Sampling Bayesian Ancestors. . .", d, m
        runid = get_runid(d,m)
        for a in ap.params["ingroup"]:
            [start,stop] = get_boundary_sites( get_phylippath(d), ap.params["asrseedtaxa"][a] )
            print "\n\t .", a
            apath = d + "/asr." + runid + "/" + a + ".dat"
            command = "python SCRIPTS/asr_bayes.py --in " + apath 
            command += " --out " + d + "/asr." + runid + "/" + a + ".bayes.txt"
            command += " --id " + a
            command += " --n " + N_BAYES_SAMPLES.__str__()
            command += " --start " + start.__str__()
            command += " --stop " + stop.__str__()
            print command
            os.system(command)
