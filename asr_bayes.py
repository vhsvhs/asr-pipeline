#
# Sample Bayesian posterior ancestral sequences from amino acid
# probability distributions, a la Williams PLoS Computational Biology
#
#

import random,sys,os
from argParser import *
from tools import *
ap = ArgParser(sys.argv)
from asrpipelinedb_api import *


def cdf(state_prob):
    for tuple in state_prob:
        if tuple[0] == "-":
            return ""
    
    x = random.random()
    
    total = 0.0
    for tuple in state_prob:
        total += tuple[1]
    
    #print "41:", total
    x *= total
    
    sum = 0.0

    for tuple in state_prob:
        sum += tuple[1]
        if sum >= x:
            return tuple[0]


def sample_data(data):
    """ Returns a sampled sequence from data.
    data[site][state] = pp   """
    seq = ""
    sites = data.keys()
    sites.sort()
    for s in sites:
        seq += cdf( data[s] )
    return seq


def make_test_data():
    d = {}
    d[1] = {'A':1.0}
    d[2] = {'N':0.5,'D':0.5}
    d[3] = {'-':100}
    d[4] = {'I':0.1,'L':0.9}
    return d


def run_test():
    seqs = []
    for ii in range(0, 100000):
        seqs.append( sample_data( make_test_data() ) )
    d = {}
    d[1] = {}
    d[2] = {}
    d[3] = {}
    d[4] = {}
    for s in seqs:
        #print s
        for ii in range(0, s.__len__()):
            state = s[ii]
            site = ii+1
            if state not in d[site]:
                d[site][state] = 1
            else:
                d[site][state] += 1
    print d


def run_asr_bayes(con, ap):
    #print "77: run_asr_bayes"
    cur = con.cursor()
    sql = "select id, name from AlignmentMethods"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        """Get the seed sequence, and then find the start and end sites."""
        alid = ii[0]
        alname = ii[1]
        
        sql = "select modelid, name from PhyloModels"
        cur.execute(sql)
        y = cur.fetchall()
        for jj in y:
            modelid = jj[0]
            modelname = jj[1]
        
            print "\n. Sampling Bayesian Ancestors. . .", alname, modelname
            runid = get_runid(alname,modelname)
                        
            for a in ap.params["ingroup"]:                    
                ancdir = alname + "/asr." + modelname + "/tree1"
                if False == os.path.exists( ancdir ):
                    continue
                if False == os.path.exists( ancdir + "/BAYES_SAMPLES" ):
                    os.system("mkdir " + ancdir + "/BAYES_SAMPLES" )
                for f in os.listdir(ancdir):
                    if f.__contains__(".dat"):

                        print alname + "/asr." + modelname + "/tree1/" + f            
                        data = get_pp_distro(alname + "/asr." + modelname + "/tree1/" + f)
                        if data.keys().__len__() < 1:
                            print "I found no data in the ancestral file ", f
                        n = 100                        
                        fout = open(alname + "/asr." + modelname + "/tree1/BAYES_SAMPLES/bayes." + f, "w")
                        mls = get_ml_sequence_from_file(alname + "/asr." + modelname + "/tree1/" + f)
                        shortname = f.split(".")[0]
                        fout.write(">" + shortname + "ML_ancestral_sequence\n")
                        fout.write(mls + "\n")
                        for ii in range(0, n):
                            fout.write(">" + shortname + ".posterior.sample." + alname.__str__() + "\n")
                            l = sample_data( data ) 
                            fout.write(l + "\n")



