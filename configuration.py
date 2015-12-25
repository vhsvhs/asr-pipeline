from argParser import *

# What character will be used to separate pieces of filenames?
# (e.g., mygene.muscle.fasta uses periods).
SEP = "."

DIR_nick = {}
DIR_nick["muscle"] = "muscle"
DIR_nick["prank"] = "prank"
DIR_nick["msaprobs"] = "msaprobs"
DIR_nick["mafft"] = "mafft"

def get_msa_nickname(msaname):
    if msaname in DIR_nick:
        return DIR_nick[msaname]
    else:
        return msaname