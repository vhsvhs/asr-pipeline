#
# Sample Bayesian posterior ancestral sequences from amino acid
# probability distributions, a la Williams PLoS Computational Biology
#
#

import random,sys,os
from argParser import *
from tools import *
ap = ArgParser(sys.argv)

inpath = ap.getArg("--in")
outpath = ap.getArg("--out")
startsite = ap.getOptionalArg("--start")
if startsite == False:
    startsite = 0
else:
    startsite = int(startsite)
stopsite = ap.getOptionalArg("--stop")
if stopsite == False:
    stopsite = -1
else:
    stopsite = int(stopsite)

def get_data(path):
    fin = open( inpath , "r")
    site_state_pp = {}
    for l in fin.xreadlines():
        if l.__len__() > 2:
            tokens = l.split()
            site = int(tokens[0])
            if site < startsite:
                continue # skip this site
            if site > stopsite and stopsite > 0:
                continue # skip this site
            if site not in site_state_pp:
                site_state_pp[site] = {}
            for ii in range(1,tokens.__len__()):
                if ii%2 == 1:
                    state = tokens[ii]
                    #print state
                    #print tokens, ii
                    prob = float(tokens[ii+1])
                    site_state_pp[ site ][state] = prob
                    #print site, state, prob
    return site_state_pp



def cdf(state_prob):
    #print "33:", state_prob
    if state_prob.keys().__contains__('-'):
        return ""
    
    x = random.random()
    
    total = 0.0
    for c in state_prob.keys():
        total += state_prob[c]
    
    #print "41:", total
    x *= total
    
    sum = 0.0

    for c in state_prob.keys():
        sum += state_prob[c]
        if sum >= x:
            return c
    

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


#
# main:
#

#run_test()
#exit()

data = get_data(inpath)
if data.keys().__len__() < 1:
    print "I found no data! ?"
    
n = 1
x = ap.getOptionalArg( "--n") 
if x != False:
    n = int(x)
    
id = ap.getOptionalArg("--id")
if id == False:
    id = ""

fout = open(outpath, "w")
mls = get_ml_sequence(data, start=startsite, stop=stopsite)
fout.write(">ML_ancestral_sequence\n")
fout.write(mls + "\n")
for ii in range(0, n):
    fout.write("> posterior.sample." + ii.__str__() + "." + id  + "\n")
    l = sample_data( data ) 
    fout.write(l + "\n")
