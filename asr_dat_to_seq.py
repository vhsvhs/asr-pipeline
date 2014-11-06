#
# USAGE:
# python asr_dat_to_seq.py <ancestral dat file path>
#

import sys, os

dat_to_seed = {} # key = site number in dat sequence (includes indels), value = site number in seed sequence

def getprobs(inpath):
    fin = open(inpath, "r")
    lines = fin.readlines()
    fin.close()

    site_states_probs = {}
    for l in lines:
        tokens = l.split()
        site = int(tokens[0])
        site_states_probs[ site ] = {}
        i = 1
        while i < tokens.__len__():
            s = tokens[i]
            foundgap = False
            if s == "-":
                p = 0.0
                foundgap = True
            else:
                p = float(tokens[i+1])
            if p > 1.0:
                p = 0.0
                foundgap = True
            site_states_probs[site][s] = p
            i += 2
            if foundgap:
                i = tokens.__len__() # stop early
    return site_states_probs

def get_ml_sequence(site_states_probs):
    mlseq = ""
    sites = site_states_probs.keys()
    sites.sort()
    for site in sites:
        maxp = 0.0
        maxc = ""
        for c in site_states_probs[site]:
            #print site_states_probs[site][c]
            if site_states_probs[site][c] > maxp:
                maxp = site_states_probs[site][c]
                maxc = c
        if maxc != "-":
            mlseq += maxc
    return mlseq
    #return mlseq
        
