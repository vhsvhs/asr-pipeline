from configuration import *
from dendropy import Tree
import math, re, subprocess, os
from Bio import Phylo # Note, this must be version 1.63 or newer.

from argParser import *
ap = ArgParser(sys.argv)

def run_script(path):
    exe = None
    if ap.params["usempi"]:
        exe = ap.params["mpirun_exe"] + " " + path
    else:
        exe = ap.params["run_exe"] + " " + path
    
    # to-do: why does version 1 fail on the Johnson lab server?
    
    #
    # version 1
    #
    #args = exe.split()
    #proc = subprocess.Popen( args, preexec_fn=os.setsid ) # see http://pymotw.com/2/subprocess/
    #proc.wait()
    
    #
    # version 2
    #
    print "\n. ASR pipeline: Running the system command:"
    print exe
    os.system(exe)

def run_subprocess(command):
    args = command.split()
    proc = subprocess.Popen( args, preexec_fn=os.setsid)
    proc.wait()
    return proc

def get_mean(values):
    sum = 0.0
    for v in values:
        sum += float(v)
    return sum / float(values.__len__())

def get_sd(values):
    mean = get_mean(values)
    sumofsquares = 0.0
    for v in values:
        sumofsquares += (v - mean)**2
    return math.sqrt( sumofsquares / float(values.__len__()) )


def get_runid(dir, model):
    nick = DIR_nick[dir]
    runid = nick + "." + model
    return runid

def get_phylippath(DIR):
    nick = DIR_nick[DIR]
    return DIR + "/" + ap.params["geneid"] + SEP + nick + SEP + "phylip"

def get_full_phylippath(DIR):
    nick = DIR_nick[DIR]
    return DIR + "/" + ap.params["geneid"] + SEP + nick + SEP + "full" + SEP + "phylip"

def get_fastapath(DIR):
    nick = DIR_nick[DIR]
    return DIR + "/" + ap.params["geneid"] + SEP + nick + SEP + "fasta"

def get_full_fastapath(DIR):
    nick = DIR_nick[DIR]
    return DIR + "/" + ap.params["geneid"] + SEP + nick + SEP + "full" + SEP + "fasta"

def get_asr_fastapath(DIR):
    return get_fastapath(DIR)

def get_asr_phylippath(DIR):
    return get_phylippath(DIR)

def get_phylipstats(path):
    """Input: a path to a phylip-formatted alignment. Returns tuple (ntaxa, nsites)"""
    fin = open(path, "r")
    header = fin.readline()
    fin.close()
    tokens = header.split()
    ntaxa = int(tokens[0])
    nsites = int(tokens[1])
    return (ntaxa, nsites)

def get_raxml_infopath(DIR, model):
    runid = get_runid(DIR,model)
    return DIR + "/RAxML_info." + runid

def get_raxml_logpath(DIR, model):
    runid = get_runid(DIR,model)
    return DIR + "/RAxML_log." + runid

#
# The path to the RAxML ML tree
#
def get_raxml_treepath(DIR, runid):
    return DIR + "/RAxML_result." + runid

#
# The path to the ML tree with ALR branch support.  These ALR values are
# calculated from ALRT values, generated by PhyML
#
def get_alrt_treepath(DIR, model):
    phylippath = get_phylippath(DIR)
    return phylippath + "_phyml_tree_" + model + ".alrt.txt"

def get_alr_treepath(DIR, model):
    phylippath = get_phylippath(DIR)
    return phylippath + "_phyml_tree_" + model + ".alr.tre"

def get_tree_length(path):
    """Input: path to newick tree. Returns the sum of branches on the tree."""
    t = Tree()
    t.read_from_path(path, "newick")
    return t.length()

def get_tree_newick(path):
    """Input: path to newick tree. Returns the tree's Newick string"""
    t = Tree()
    t.read_from_path(path, "newick")
    edges = t.get_edge_set()
    for e in edges:
        e.length = e.length = float( "%.4f"%e.length )
    return t.__str__()

def reroot_tree(tstr):
    """Input: a tree path to a Newick tree.  Output: a re-rooted version of the tree, based on the outgroup defined in configuration.py"""
    t = Tree()
    #print tstr
    t.read_from_string(tstr, "newick")
    #print "t:", t.__str__()
    og = ap.params["outgroup"]
    og = re.sub("\[", "", og)
    og = re.sub("\]", "", og)
    og = re.sub("\"", "", og)
    ogs = og.split(",")
    mrca = t.mrca(taxon_labels=ogs)
    t.reroot_at_edge(mrca.edge, update_splits=False)
    ret = t.as_string("newick")
    ret = re.sub("\[\&\R\] ", "", ret)
    ret = ret.strip()
    return ret


def get_cladogram_path(d, model):
    tpath = d + "/asr." + model + "/tree1/tree1.txt"
    fin = open(tpath, "r")
    cline = fin.readlines()[3]
    cline = cline.strip()
    cstr = re.sub(" ", "", cline)
    cstr = re.sub(";", ";", cstr)
    cstr = reroot_tree( cstr )
    
    cladopath = d + "/asr." + model + "/cladogram.tre"
    fout = open(cladopath, "w")
    fout.write( cstr + "\n")
    fout.close()
    return cladopath

def get_sequence(msapath, taxa):
    """msapath must be a phylip file.  Returns the seed sequence."""
    fin = open(msapath, "r")
    for l in fin.readlines():
        if l.startswith(taxa):
            tokens = l.split()
            return tokens[1]

def get_ml_sequence(site_states_probs, start=0, stop=-1):
    mlseq = ""
    sites = site_states_probs.keys()
    sites.sort()
    for site in sites:
        if site < start:
            continue
        if site > stop and stop > 0:
            continue
        maxp = 0.0
        maxc = ""
        for tup in site_states_probs[site]:
            #print site_states_probs[site][c]
            state = tup[0]
            p = tup[1]
            if  p > maxp:
                maxp = p
                maxc = state
        if maxc != "-" and maxc != '-':
            mlseq += maxc.upper()
        #print site_states_probs, site, maxc
    return mlseq

def get_ml_sequence_from_file(path, getindels=False):
    fin = open(path, "r")
    mlseq = ""
    for l in fin.xreadlines():
        if l.__len__() > 3:
            tokens = l.split()
            state = tokens[1]
            if state != "-":
                mlseq += state.upper()
            elif getindels:
                mlseq += "-"
    return mlseq

def get_pp_distro(path):
    fin = open( path , "r")
    site_state_pp = {}
    for l in fin.xreadlines():
        if l.__len__() > 2:
            tokens = l.split()
            site = int(tokens[0])
            if site not in site_state_pp:
                site_state_pp[site] = []
            for ii in range(1,tokens.__len__()):
                if ii%2 == 1:
                    state = tokens[ii].upper()
                    #print state
                    #print tokens, ii
                    prob = float(tokens[ii+1])
                    site_state_pp[ site ].append( [state,prob] )
    return site_state_pp

def get_model_path(model, ap):
    modelstr = "~/Applications/paml44/dat/lg.dat"
    if model.__contains__("JTT"):
        modelstr = ap.params["mmfolder"] + "/jones.dat"
    elif model.__contains__("WAG"):
        modelstr = ap.params["mmfolder"] + "/wag.dat"
    elif model.__contains__("LG"):
        modelstr = ap.params["mmfolder"] + "/lg.dat"
    return modelstr

def get_pp_distro_stats(data):
    """Input: the output from get_pp_distro.  Output: mean and s.d. PP."""
    pps = []
    for site in data:
        pps.append(data[site][1])
    sum = 0.0

# returns a bin number for this P value
def binForProb(p):
    return int(p / 0.05)

# return the P value of the floor of this bin
def probForBin(b):
    x = float(b*5) / float(100)
    if x == 1.00:
        return x
    return x + 0.025
    
def get_boundary_sites(msapath, taxa):
    """Input: a sequence alignment, a seed taxa, and start/end motifs determined in the config. file.
    This function then finds those motifs in the seed sequence (which includes gaps),
    and then returns the translate start/end sites for this alignment.
    NOTE: msapath must point to a PHYLIP-formatted alignment."""
    start_motif = ap.params["start_motif"] #"YQLI"
    end_motif = ap.params["end_motif"] #"MPFF"

    #print "248:", start_motif
    #print "249:", end_motif

    seq = get_sequence(msapath, taxa)
    #print seq
    startsite = 1
    endsite = seq.__len__()
                                                                  
    if start_motif != None:
        if start_motif.__len__() > 0:
            for i in range(0, seq.__len__()):
                #print "258:", i, seq[i], start_motif[0]
                if seq[i] == start_motif[0]:
                    here = ""
                    j = i
                    while here.__len__() < start_motif.__len__() and j < seq.__len__():
                        #print "262:", j, here
                        if seq[j] != "-":
                            here += seq[j]
                        j += 1
                    
                    if here  == start_motif:
                        startsite = i + 1         
                        break
        
    if end_motif != None:
        if end_motif.__len__() > 0:
            for i in range(i, seq.__len__()):
                if seq[i] == end_motif[0]:
                    here = ""
                    j = i
                    while here.__len__() < end_motif.__len__() and j < seq.__len__():
                        if seq[j] != "-":
                            here += seq[j]
                        j += 1          
                    if here  == end_motif:
                        endsite = j
                        break
    return [startsite, endsite]
