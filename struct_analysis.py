import os, sys, time
from configuration import *
from tools import *
from phyloxml_helper import *
from asr_dat_to_seq import *
from asrpipelinedb_api import *

PDBDIR = "pdb"

def write_ancseq_fasta(con, ap):
    """Writes a FASTA file into the PDBDIR, containing
    the sequences for all the ancestors listed with an ingroup
    in the configuration file."""
    if os.path.exists(PDBDIR) == False:
        os.system("mkdir " + PDBDIR)
    
    fout = open(PDBDIR + "/ancseqs.fasta", "w")
            
    for model in get_phylo_modelnames(con):
        for msa in get_alignment_method_names(con):
            for anc in ap.params["ingroup"]:
                datpath = msa + "/asr." + model + "/" + anc + ".dat"
                probs = getprobs(datpath)
                mls = get_ml_sequence(probs)
                fout.write(">" + datpath + "\n")
                fout.write(mls + "\n")
    fout.close()