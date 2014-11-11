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
from asrpipelinedb import *
from log import *

""""
The original shelll script....

mkdir MUSCLE
mkdir PRANK
mkdir MSAPROBS
cp cgmc.erg.raw.fasta cgmc.txt
python SCRIPTS/clean_seqs.py cgmc.txt > cgmc.fasta
mpirun -np 4 SCRIPTS/run_msas.mpi.sh
"""
from asrpipelinedb_api import *

def read_fasta(fpath):
    """Ignores redundant taxa."""
    fin = open(fpath, "r")    
    taxanames = []
    taxa_seq = {}
    last_taxa = None
    last_seq = ""
    okay = True
    for l in fin.xreadlines():
        l = l.strip()
        if l.__len__() <= 1:
            pass
        elif l.startswith(">"):
            okay = True
            taxaname = re.sub(">", "", l.split()[0] )
            if taxaname in taxa_seq:
                okay = False
            else:
                if last_taxa != None:
                    taxa_seq[last_taxa] = last_seq
                last_taxa = taxaname
                last_seq = ""
                
        elif okay == True:
            last_seq += l
    taxa_seq[last_taxa] = last_seq
    fin.close()
    return taxa_seq

def import_and_clean_erg_seqs(con, ap):    
    taxa_seq = read_fasta(ap.params["ergseqpath"])
    
    for taxon in taxa_seq:
        import_original_seq(con, taxon, taxa_seq[taxon] )
    
    cleanpath = ap.params["ergseqpath"]
    fout = open(cleanpath, "w")
    for taxon in taxa_seq:
        fout.write(">" + taxon + "\n")
        fout.write(taxa_seq[taxon] + "\n")
    fout.close()
    
    cur = con.cursor()
    sql = "select count(*) from Taxa"
    cur.execute(sql)
    c = cur.fetchone()[0]
    write_log(con, "There are " + c.__str__() + " taxa in the database.")
    
def write_msa_commands(con, ap):
    p = "SCRIPTS/msas.commands.sh"
    fout = open(p, "w")
    for msa in ap.params["msa_algorithms"]:
        if msa == "muscle":
            fout.write(ap.params["muscle_exe"] + " -in " + ap.params["ergseqpath"] + " -out " + get_fastapath(msa) + "\n")
            import_alignment_method(con, msa, ap.params["muscle_exe"])
        elif msa == "prank":
            fout.write(ap.params["prank_exe"] + " -d=" + ap.params["ergseqpath"] + " -o=" + get_fastapath(msa) + "\n")
            import_alignment_method(con, msa, ap.params["prank_exe"])
        elif msa == "msaprobs": 
            fout.write(ap.params["msaprobs_exe"] + " " + ap.params["ergseqpath"] + " > " + get_fastapath(msa) + "\n")
            import_alignment_method(con, msa, ap.params["msaprobs_exe"])
        elif msa == "mafft": 
            fout.write(ap.params["mafft_exe"] + " --auto " + ap.params["ergseqpath"] + " > " + get_fastapath(msa) + "\n")
            import_alignment_method(con, msa, ap.params["mafft_exe"])
    fout.close()
    return p
    #os.system("mpirun -np 4 --machinefile hosts.txt /common/bin/mpi_dispatch SCRIPTS/msas.commands.sh")

def check_aligned_sequences(con, ap):
    cur = con.cursor()
    for msa in ap.params["msa_algorithms"]:
        expected_fasta = get_fastapath(msa)
        if False == os.path.exists( expected_fasta ):
            write_error(con, "I can't find the expected FASTA output from " + msa + " at " + expected_fasta, code=None)
        taxa_seq = read_fasta(expected_fasta)
    
        sql = "select id from AlignmentMethods where name='" + msa + "'"
        cur.execute(sql)
        x = cur.fetchall()
        if x.__len__() == 0:
            write_error(con, "The alignment method " + msa + " does not exist in the database.")
            exit()
        almethodid = x[0][0] 

        for taxa in taxa_seq:
            taxonid = get_taxonid(con, taxa)
            import_aligned_seq(con, taxonid, almethodid, taxa_seq[taxa])
                 
        cur = con.cursor()
        sql = "select count(*) from AlignedSequences where almethod=" + almethodid.__str__()
        cur.execute(sql)
        c = cur.fetchone()[0]
        write_log(con, "There are " + c.__str__() + " sequences in the " + msa + " alignment.")
     

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
        f = get_trimmed_fastapath(msa)
        p = get_trimmed_phylippath(msa)
        fasta_to_phylip(f, p)


def build_seed_sitesets(con, ap):
    cur = con.cursor()
    sql = "insert or replace into SiteSets(name) VALUES(\"seed\")"
    cur.execute(sql)
    con.commit()
    

def clear_sitesets(con, ap):
    cur = con.cursor()
    sql = "delete from SiteSets"
    cur.execute(sql)
    sql = "delete from SiteSetsAlignment"
    cur.execute(sql)
    con.commit()

def trim_alignments(con, ap):
    """Trims the alignment to match the start and stop motifs.
    The results are written as both PHYLIP and FASTA formats."""
    cur = con.cursor()
    
    """Create a SiteSet for the sites corresponding to the seed sequence(s)."""
    sql = "insert or replace into SiteSets(setname) VALUES('seed')"
    cur.execute(sql)
    con.commit()
    sql = "select id from SiteSets where setname='seed'"
    cur.execute(sql)
    sitesetid = cur.fetchone()[0]
    
    """Now iterate through each alignment, and find the start/stop sites
    corresponding to the seed sequence(s)."""
    sql = "select id, name from AlignmentMethods"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        """Get the seed sequence, and then find the start and end sites."""
        alid = ii[0]
        alname = ii[1]        
        seedid = get_taxonid(con, ap.params["seed_motif_seq"] )
        seedseq = get_aligned_seq(con, seedid, alid)
        start_motif = ap.params["start_motif"]
        end_motif = ap.params["end_motif"]
        [start, stop] = get_boundary_sites(seedseq, start_motif, end_motif)
        
        """Remember these start and end sites."""
        sql = "insert or replace into SiteSetsAlignment(setid, almethod, fromsite, tosite) VALUES("
        sql += sitesetid.__str__() + "," + alid.__str__() + "," + start.__str__() + "," + stop.__str__() + ")"
        cur.execute(sql)
        con.commit()
        
        """Write an alignment, trimmed to the seed."""
        ffout = open( get_trimmed_fastapath(alname), "w")
        sql = "select taxonid, alsequence from AlignedSequences where almethod=" + alid.__str__()
        cur.execute(sql)
        y = cur.fetchall()
        for jj in y: # for each sequence
            taxonid = jj[0]
            taxa = get_taxon_name(con, taxonid)
            seq = jj[1]
            trimmed_seq = seq[start-1:stop]
            hasdata = False 
            for c in trimmed_seq:
                if c != "-":
                    hasdata = True
            if hasdata == True: # Does this sequence contain at least one non-indel character?
                ffout.write(">" + taxa + "\n" + trimmed_seq + "\n")
            else:
                write_log(con, "After trimming the alignment " + alid.__str__() + " to the seed boundaries, the taxa " + taxa + " contains no sequence content. This taxa will be obmitted from downstream analysis.")
        ffout.close()
    
    #trim_by_zorro(ap)
    
def build_restriction_alignments(con, ap):
    #
    # continue here
    #
    
    # 1. run ZORRO
    # 2. parse the results
    # 3. optimize the set of sites we draw from ZORRO, and the stability of the output tree
    # 4. write the final alignment
    pass
    

def trim_by_zorro(ap):
    z_commands = []
    for msa in ap.params["msa_algorithms"]:
        c = "zorro -sample " + get_fastapath(msa) + " > " + get_fastapath(msa) + ".zorro"
        z_commands.append( c )
    
    fout = open("SCRIPTS/zorro_commands.sh", "w")
    for c in z_commands:
        fout.write( c + "\n" )
    fout.close()
    run_script("SCRIPTS/zorro_commands.sh")

    for msa in ap.params["msa_algorithms"]:
        fin = open(get_fastapath(msa) + ".zorro", "r")
        keep_sites = []
        site = 0
        for ll in fin.xreadlines():
            if ll.__len__() > 2:
                score = float( ll.strip() )
                if score > 2.0: # to-do, fix this:
                    keep_sites.append( site )
            site += 1
        fin.close()
        
        fout = open( get_trimmed_fastapath(msa), "w")
        fin = open( get_fastapath(msa), "r")
        site_count = 0
        for l in fin.xreadlines():
            if l.startswith(">"):
                fout.write("\n" + l)
                site_count = 0
            elif l.__len__() > 5:
                l = l.strip()
                for c in l:
                    if site_count in keep_sites:
                        fout.write(c)
                    site_count += 1
        fin.close()
        fout.close()
        
    
    
