import re,sys,os
from tools import *

def read_config_file(ap):
    cpath = ap.getArg("--configpath")
    fin = open(cpath, "r")

    # Default values:
    #
    ap.params["end_motif"] = None
    ap.params["start_motif"] = None

    for l in fin.xreadlines():
        l = l.strip()
        if l.startswith("#"):
            continue
        if l.__len__() < 2:
            continue
        tokens = l.split("=")
        if tokens.__len__() < 1: 
            continue
        
        if tokens[0].startswith("GENE_ID"):
            ap.params["geneid"] = re.sub(" ", "", tokens[1])
        
        if tokens[0].startswith("PROJECT_TITLE"):
            ap.params["project_title"] = re.sub(" ", "", tokens[1])

        if tokens[0].startswith("SEQUENCES"):
            ap.params["ergseqpath"] = re.sub(" ", "", tokens[1])
        
        elif tokens[0].startswith("RAXML"):
            ap.params["raxml_exe"] = tokens[1].strip()
        
        elif tokens[0].startswith("PHYML"):
            ap.params["phyml_exe"] = tokens[1].strip()

        elif tokens[0].startswith("LAZARUS"):
            ap.params["lazarus_exe"] = tokens[1].strip()

        elif tokens[0].startswith("MARKOV_MODEL_FOLDER"):
            ap.params["mmfolder"] = tokens[1].strip()
        
        elif tokens[0].startswith("MPIRUN"):
            ap.params["mpirun_exe"] = tokens[1]
        
        elif tokens[0].startswith("RUN"):
            ap.params["run_exe"] = tokens[1]
        
        elif tokens[0].startswith("MSAPROBS"):
            ap.params["msaprobs_exe"] = tokens[1].strip()

        elif tokens[0].startswith("MUSCLE"):
            ap.params["muscle_exe"] = tokens[1].strip()

        elif tokens[0].startswith("PRANK"):
            ap.params["prank_exe"] = tokens[1].strip()
        
        elif tokens[0].startswith("ALIGNMENT_ALGORITHMS"):
            x = tokens[1].split()
            ap.params["msa_algorithms"] = []
            for i in x:
                ap.params["msa_algorithms"].append( i )
        
        elif tokens[0].startswith("MODELS_RAXML"):
            x = tokens[1].split()
            ap.params["raxml_models"] = []
            for i in x:
                ap.params["raxml_models"].append( i )
        
        elif tokens[0].startswith("START_MOTIF"):
            ap.params["start_motif"] = re.sub(" ", "", tokens[1])
        
        elif tokens[0].startswith("END_MOTIF"):
            ap.params["end_motif"] = re.sub(" ", "", tokens[1])
        
        elif tokens[0].startswith("SEED_MOTIF_TAXA"):
            ap.params["seed_motif_seq"] = re.sub(" ", "", tokens[1])
        
        elif tokens[0].startswith("N_BAYES_SAMPLES"):
            ap.params["n_bayes_samples"] = int(tokens[1])
        
        elif tokens[0].startswith("OUTGROUP"):
            ap.params["outgroup"] = re.sub(" ", "", tokens[1])
        
        elif tokens[0].startswith("ANCESTORS"):
            x = tokens[1].split()
            ap.params["ancestors"] = []
            for i in x:
                ap.params["ancestors"].append( i )
        
        elif tokens[0].startswith("INGROUP"):
            if "ingroup" not in ap.params:
                ap.params["ingroup"] = {}
            anc = tokens[0].split()[1]
            ingroup = tokens[0].split()[2]
            ap.params["ingroup"][ anc ] = ingroup
        
        elif tokens[0].startswith("SEED"):
            if "seedtaxa" not in ap.params:
                ap.params["seedtaxa"] = {}
            anc = tokens[0].split()[1]
            seed = tokens[0].split()[2]
            ap.params["seedtaxa"][ anc ] = re.sub(" ", "", seed)

        elif tokens[0].startswith("COMPARE"):
            anc1 = tokens[0].split()[1]
            anc2 = tokens[0].split()[2]
            if "compareanc" not in ap.params:
                ap.params["compareanc"] = []
            ap.params["compareanc"].append( (anc1,anc2) )
            
        elif tokens[0].startswith("USE_MPI"):
            print tokens[0]
            answer = re.sub(" ", "", tokens[0].split()[1])
            if answer == "on" or answer == "True":
                ap.params["usempi"] = True
            else:
                ap.params["usempi"] = False
    fin.close()

def verify_config(ap):
    """Will return nothing if the configuration is ok.  Will error and quit if the
    configuration is flawed."""
    if "ancestors" in ap.params:
        for a in ap.params["ancestors"]:
            print a, ap.params["seedtaxa"]
            if a not in ap.params["seedtaxa"]:
                print "\n. ERROR: You did not specify a SEED for the ancestor", a
                exit()
        for c in ap.params["compareanc"]:
            a1 = c[0]
            a2 = c[1]
            if a1 not in ap.params["ancestors"]:
                print "\n. ERROR: you specified a comparison between ancestors", a1, "and", a2, "but", a1,"was not defined in the ANCESTORS line."
                exit() 
            if a1 not in ap.params["ancestors"]:
                print "\n. ERROR: you specified a comparison between ancestors", a1, "and", a2, "but", a2,"was not defined in the ANCESTORS line."
                exit()
                
    if False == os.path.exists(ap.params["ergseqpath"]):
        print "\n. I could not find your sequences at", ap.params["ergseqpath"]
        exit()

    if ap.params["start_motif"] == None:
        ap.params["start_motif"] = ""
    if ap.params["end_motif"] == None:
        ap.params["end_motif"] = ""
        
    for msa in ap.params["msa_algorithms"]:
        if msa == "MUSCLE" and "muscle_exe" not in ap.params:
            print "\n. Something is wrong. Your config file doesn't have an executable path for MUSCLE."
        if msa == "MSAPROBS" and "msaprobs_exe" not in ap.params:
            print "\n. Something is wrong. Your config file doesn't have an executable path for MSAPROBS."
        if msa == "PRANK" and "prank_exe" not in ap.params:
            print "\n. Something is wrong. Your config file doesn't have an executable path for PRANK."
        
    return ap

def print_config(ap):
    for p in ap.params:
        print p, ":", ap.params[p]


def setup_workspace(ap):
    for msa in ap.params["msa_algorithms"]:
        if False == os.path.exists(msa):
            os.system("mkdir " + msa)
    if False == os.path.exists("SCRIPTS"):
        os.system("mkdir SCRIPTS")
            