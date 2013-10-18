#
#
#
import random,sys,os

#from configuration import *
from tools import *

for d in DIRS:
    for m in models:
        print "\n. Sampling Bayesian Ancestors. . .", d, m
        runid = get_runid(d,m)
        for a in ingroups:
            [start,stop] = get_boundary_sites( get_fullphylippath(d), ingroup_seed[a] )
            print "\n\t .", a
            apath = d + "/asr." + runid + "/anc." + a + ".dat"
            command = "python SCRIPTS/asr_bayes.py --in " + apath 
            command += " --out " + d + "/asr." + runid + "/anc." + a + ".bayes.txt"
            command += " --id " + a
            command += " --n " + N_BAYES_SAMPLES.__str__()
            command += " --start " + start.__str__()
            command += " --stop " + stop.__str__()
            print command
            os.system(command)