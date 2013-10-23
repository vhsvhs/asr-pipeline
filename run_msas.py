#
# USAGE:
#
# $> python run_msas.py P N
#
# P is the path to a fasta-formatted file with sequences to be aligned.
# N is the short nickname of the analysis.
#
from argParser import *
from tools import *

""""
The original shelll script....

mkdir MUSCLE
mkdir PRANK
mkdir MSAPROBS
cp cgmc.erg.raw.fasta cgmc.txt
python SCRIPTS/clean_seqs.py cgmc.txt > cgmc.fasta
mpirun -np 4 SCRIPTS/run_msas.mpi.sh
"""

def write_msa_commands(ap):
    p = "SCRIPTS/msas.commands.sh"
    fout = open(p, "w")
    for msa in ap.params["msa_algorithms"]:
        if msa == "muscle":
            fout.write(ap.params["muscle_exe"] + " -in " + ap.params["ergseqpath"] + " -out " + get_fastafull_path(msa, ap) + "\n")
        elif msa == "prank":
            fout.write(ap.params["prank_exe"] + " -d=" + ap.params["ergseqpath"] + " -o=" + get_fastafull_path(msa, ap) + "\n")
        elif msa == "msaprobs": 
            fout.write(ap.params["msaprobs_exe"] + " " + ap.params["ergseqpath"] + " > " + get_fastafull_path(msa, ap) + "\n")
    fout.close()
    return p
    #os.system("mpirun -np 4 --machinefile hosts.txt /common/bin/mpi_dispatch SCRIPTS/msas.commands.sh")

def convert_all_fasta_to_phylip(ap):
    for msa in ap.params["msa_algorithms"]:
        f = get_fastafull_path(msa, ap)
        fasta_to_phylip(f)

def trim_alignments(ap):
    """Trims the alignment to match the start and stop motifs.
    The results are written as both PHYLIP and FASTA formats."""
    convert_all_fasta_to_phylip(ap)
    for msa in ap.params["msa_algorithms"]:
        fin = open(get_phylipfull_path(msa, ap), "r")
        lines = fin.readlines()
        ntaxa = lines[0].split()[0]
        seqlen = int( lines[0].split()[1] )
        start = 1
        stop = seqlen
        
        if "start_motif" in ap.params and "end_motif" in ap.params:
            [start,stop] = get_boundary_sites( get_phylipfull_path(msa, ap), ap)
        
        pfout = open(get_phylip_path(msa,ap), "w")
        pfout.write(ntaxa + "  " + (stop-start+1).__str__() + "\n")
        
        ffout = open(get_fasta_path(msa,ap), "w")
        
        for ii in range(1, lines.__len__()):
            tokens = lines[ii].split()
            taxa = tokens[0]
            seq = tokens[1]
            trimmed_seq = seq[start-1:stop]
            pfout.write(taxa + "  " + trimmed_seq + "\n")
            ffout.write(">" + taxa + "\n" + trimmed_seq + "\n")
        pfout.close()
        ffout.close()