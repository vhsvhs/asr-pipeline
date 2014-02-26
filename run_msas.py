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

def clean_erg_seqs(ap):
    fin = open(ap.params["ergseqpath"], "r")
    lines = fin.readlines()
    fin.close()
    
    outlines = []
    taxanames = []
    for l in lines:
        l = l.strip()
        if l.__len__() <= 1:
            pass
        elif l.startswith(">"):
            taxaname = re.sub(">", "", l.split()[0] )
            if taxaname in taxanames:
                print "Something is wrong. I found the sequence name", taxaname, "twice in your sequences."
                exit()
            outlines.append(">" + taxaname)
        else:
            outlines.append(l)
    
    cleanpath = ap.params["ergseqpath"]
    fout = open(cleanpath, "w")
    for l in outlines:
        fout.write(l + "\n")
    fout.close()
    
def write_msa_commands(ap):
    p = "SCRIPTS/msas.commands.sh"
    fout = open(p, "w")
    for msa in ap.params["msa_algorithms"]:
        if msa == "muscle":
            fout.write(ap.params["muscle_exe"] + " -in " + ap.params["ergseqpath"] + " -out " + get_fastapath(msa) + "\n")
        elif msa == "prank":
            fout.write(ap.params["prank_exe"] + " -d=" + ap.params["ergseqpath"] + " -o=" + get_fastapath(msa) + "\n")
        elif msa == "msaprobs": 
            fout.write(ap.params["msaprobs_exe"] + " " + ap.params["ergseqpath"] + " > " + get_fastapath(msa) + "\n")
    fout.close()
    return p
    #os.system("mpirun -np 4 --machinefile hosts.txt /common/bin/mpi_dispatch SCRIPTS/msas.commands.sh")

def fasta_to_phylip(inpath, outpath):
    fin = open(inpath, "r")
    last_taxa = None
    taxa_seq = {}
    for l in fin.xreadlines():
        if l.startswith(">"):
            taxa = re.sub(">", "", l)
            taxa = taxa.split()[0] # trim any lingering GenBank labels. As of 03/2014, this line may be unecessary - split()[0] part of the parsing occurs when the RAW sequences are read into a FASTA file
            taxa = taxa.strip()
            last_taxa = taxa
            taxa_seq[ taxa ] = ""
        elif l.__len__() > 0:
            l = l.strip()
            taxa_seq[ last_taxa ] += l
    fin.close()
    
    fout = open(outpath, "w")
    fout.write("  " + taxa_seq.__len__().__str__() + "  " + taxa_seq[ taxa_seq.keys()[0] ].__len__().__str__() + "\n")
    for taxa in taxa_seq:
        fout.write(taxa + "   " + taxa_seq[taxa] + "\n")
    fout.close()

def convert_all_fasta_to_phylip(ap):
    for msa in ap.params["msa_algorithms"]:
        f = get_fastapath(msa)
        p = get_phylippath(msa)
        fasta_to_phylip(f, p)

#
# depricated method:
#
def trim_alignments(ap):
    """Trims the alignment to match the start and stop motifs.
    The results are written as both PHYLIP and FASTA formats."""
    for msa in ap.params["msa_algorithms"]:
        fin = open(get_phylippath(msa), "r")
        lines = fin.readlines()
        ntaxa = lines[0].split()[0]
        seqlen = int( lines[0].split()[1] )
        start = 1
        stop = seqlen
        
        if (ap.params["end_motif"] != None) and (ap.params["start_motif"] != None):
            [start,stop] = get_boundary_sites( get_phylippath(msa), ap)
        

        """Write PHYLIP to pfout and FASTA to ffout."""
        poutl = ""
        ffout = open(get_fasta_path(msa,ap), "w")
        count_good_taxa = 0
        for ii in range(1, lines.__len__()):
            tokens = lines[ii].split()
            taxa = tokens[0]
            seq = tokens[1]
            trimmed_seq = seq[start-1:stop]
            hasdata = False
            for c in trimmed_seq:
                if c != "-":
                    hasdata = True
            if hasdata == True:
                count_good_taxa += 1
                poutl += taxa + "  " + trimmed_seq + "\n"
                ffout.write(">" + taxa + "\n" + trimmed_seq + "\n")
        ffout.close()
        pfout = open(get_phylippath(msa,ap), "w")
        pfout.write(count_good_taxa.__str__() + "  " + (stop-start+1).__str__() + "\n")
        pfout.write(poutl)
        pfout.close()
