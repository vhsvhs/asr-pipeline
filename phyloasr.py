#################################################################################
#
# run_phyloasr.py
#
# Version 1.0 - September 20, 2013
# VIctor Hanson-Smith
# victor.hanson-smith@ucsf.edu
#
# This script reads multiple sequence alignments, infers ML phylogenies
# using RAxML and PhyML, and then reconstructs ancestral sequences on
# those phylogenies, using Lazarus.
# Optionally, this program can compare ancestors using the 'hb' metric.
#
# USAGE:
# python run_phyloasr.py ID
#
# . . . where ID is the unique ID name of the gene family.
# This script assumes that other scripts, previously, created multiple
# sequence alignments with names such as ID.muscle (for the MUSCLE software)
# and ID.prank (for the PRANK software).
#
# If you choose to modify this script, please edit the header definitions,
# indicated below, to point to correct paths and directories.
#
# If you have questions, please contact Victor @ victor.hanson-smith@ucsf.edu
#

from tools import *



#
# Create RAxML commands for each alignment and each model.
#
def write_raxml_commands(ap):
    commands = []
    for msa in ap.params["msa_algorithms"]:
        phypath = get_phylip_path(msa, ap)
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model) 
            command = ap.params["raxml_exe"]
            command += "  -s " + phypath
            command += " -n " + runid
            command += " -e 0.001"
            command += " -m " + model
            command += " -p 12345"
            command += " > OUT." + msa + "/catch." + runid + ".txt" 
            commands.append(command)
    p = "SCRIPTS/raxml.commands.sh"
    fout = open(p, "w")
    for c in commands:
        fout.write(c + "\n")
    fout.close()
    return p
"""

#
# Run RAxML. . .
#
fout = open("raxml_commands.txt", "w")
for c in commands:
    fout.write( c + "\n" )
fout.close() 
#os.system(MPIRUN + " raxml_commands.txt")

for DIR in DIRS:
    os.system("mv ./*" + DIR_nick[DIR] + "* ./" + DIR + "/")

#
# Post-RAxML cleanup:
#
for DIR in DIRS:
    keyword = DIR_nick[DIR]
    for runid in DIR_runids[DIR]:
        # Renmae the RAxML tree with ".tre" suffix, so that future programs
        # will know it's a newick tree.
        os.system("mv " + get_raxml_treepath(DIR, runid) + "  " + get_ml_treepath(DIR, runid) )
        # Removed the .reduced alignment that's written by RAxML.
        # It will be redundant with our trimmed version of the alignment.
        os.system("rm " + get_phylippath(DIR) + ".reduced")


#
# Fetch ML results,
# Calculate min, max likelihoods for each model, 
# and find the best-fitting model.
#
DIR_runid_lnl = {} 
DIR_suml = {}
DIR_runid_pp = {}
DIR_runid_alpha = {}
for DIR in DIR_runids:
    DIR_runid_lnl[DIR] = {}
    DIR_runid_pp[DIR] = {}
    DIR_runid_alpha[DIR] = {}
    maxl = None
    minl = None
    for runid in DIR_runids[DIR]:
        # Get the output stats about this run:
        statspath = DIR_runid_statspath[DIR][runid]
        if False == os.path.exists(statspath):
            print "I can't find:", statspath
            continue
        fin = open(statspath, "r")
        lines = fin.readlines()
        for l in lines:
            if l.startswith("Final GAMMA-based Score of best tree"):
                tokens = l.split()
                lnl = float(tokens[6])
                DIR_runid_lnl[DIR][runid] = lnl
                print runid, lnl
                if maxl == None:
                    maxl = lnl
                if minl == None:
                    minl = lnl
                if lnl > maxl:
                    maxl = lnl
                if lnl < minl:
                    minl = lnl
            elif l.startswith("alpha[0]:"):
                l = l.strip()
                tokens = l.split()
                this_alpha = tokens[1]
                DIR_runid_alpha[DIR][runid] = this_alpha
        fin.close()

    DIR_suml[DIR] = 0.0
    for runid in DIR_runid_lnl[DIR]:
        DIR_suml[DIR] += (DIR_runid_lnl[DIR][runid] - minl)
    
    fout = open("lnl_log." + DIR_nick[DIR] + ".txt", "w")
    for runid in DIR_runid_lnl[DIR]:
        lnl = DIR_runid_lnl[DIR][runid]
        pp = (lnl-minl)/(DIR_suml[DIR])
        DIR_runid_pp[DIR][runid] = pp
        special = ""
        if lnl == maxl:
            special = " (ML) "
        line = DIR + "\t" + runid + "\t" + lnl.__str__() + "\t" + "%.4f"%pp + "\t" + special
        fout.write(line + "\n")
    fout.close()


#
# Use PhyML to calculate ALRTs on each RAxML ML tree.
#
alrt_commands = []
for DIR in DIR_runids:
    for runid in DIR_runid_lnl[DIR]:
        mltreepath = get_ml_treepath(DIR, runid)
        phylippath = get_phylippath(DIR)
        model = ""
        if runid.__contains__("JTT"):
            model = "JTT"
        elif runid.__contains__("WAG"):
            model = "WAG"
        elif runid.__contains__("LG"):
            model = "LG"
        gamma = ""
        if runid.__contains__("GAMMA"):
            this_alpha = DIR_runid_alpha[DIR][runid]
            gamma = " --nclasses 8 --alpha " + this_alpha.__str__()
        command = PHYML
        command += " --input " + phylippath
        command += " -u " + mltreepath
        command += " --model " + model
        command += gamma
        command += " -f m"  # use fixed pi values
        command += " -o n"  # don't optimize anything
        command += " -b -1" # calculate ALRTs
        command += " --run_id " + runid + ".alrt"
        command += " --datatype aa"
        command += " > " + DIR + "/catch.phyml." + runid + ".txt"
        alrt_commands.append( command )
fout = open("alrt_commands.sh", "w")
for c in alrt_commands:
    fout.write(c + "\n")
fout.close()
#os.system(MPIRUN + " alrt_commands.sh")        

#
# Convert ALRTs to ALRs
#
alr_commands = []
for DIR in DIR_runids:
    for runid in DIR_runid_lnl[DIR]:
        alrt_treepath = get_alrt_treepath(DIR, runid)
        alr_treepath = get_alr_treepath(DIR, runid)
        c = ALRT2ALR + " " + alrt_treepath + " > " + alr_treepath
        alr_commands.append(c)
#for c in alr_commands:
#    print c
#    os.system(c)

#exit()

##############################################
#
# How similar are the ML trees?
#
treepath_id = {} # key = treepath, value = integer ID for brevity.
ta_tb_distance = {} # key = treepath, value = hash, key = treepath, value = symmetric distance between the trees
ii = 1
for DIR in DIR_runids:
    for runid in DIR_runid_lnl[DIR]:
        mltreepath = get_ml_treepath(DIR, runid)
        treepath_id[ mltreepath ] = ii
        ii += 1
from dendropy import Tree
for ta in treepath_id:
    ta_tb_distance[ta] = {}
    treea = Tree()
    treea.read_from_path( ta, "newick" )
    for tb in treepath_id:
        treeb = Tree()
        treeb.read_from_path( tb, "newick" )
        d = treea.symmetric_difference( treeb )
        d = treeb.symmetric_difference( treea )
        ta_tb_distance[ta][tb] = d
print ta_tb_distance
#exit()

#################################################
#
# Run ASR
#
asr_commands = []
for DIR in DIR_runid_lnl:
    fin = open(get_fastapath(DIR), "r")
    lines = fin.readlines()
    fin.close()
    innames = [] # array of taxa in the trimmed alignment:
    for i in range(0, lines.__len__()):
        if lines[i].startswith(">"):
            l = lines[i].strip()
            name = re.sub(">", "", l)
            if name in innames:
                print ". Your alignment" + get_fastapath(DIR) + " has two copies of", name
                exit()
            print "182:", name
            innames.append(name)

    fin = open(get_fastapath(DIR), "r")
    outlines = []
    linecache = ""
    lines = fin.readlines()
    toggle = False
    for i in range(0, lines.__len__()):
        if lines[i].startswith(">"):
            l = lines[i].strip()
            name = re.sub(">", "", l)
            if name in innames:
                if linecache.__len__() > 2:
                    outlines.append(linecache)
                    linecache = ""
                toggle = True
                #print "adding", l
                outlines.append(l + "\n")
            else:
                #print "skipping", name
                toggle = False
        elif toggle == True:
            linecache += lines[i]
    if linecache.__len__() > 2:
        outlines.append(linecache)

    fout = open(DIR + "/" + gene + "." + DIR_nick[DIR] + ".asr.fasta", "w")
    for l in outlines:
        fout.write(l)
    fout.close()

    os.system("perl ~/Applications/seqConverter.pl -d" + DIR + "/" + gene + "." + DIR_nick[DIR] + ".asr.fasta -ope")

    for runid in DIR_runid_lnl[DIR]:
        os.system("mkdir " + DIR + "/asr." + runid)
        model = "~/Applications/paml44/dat/lg.dat"
        if runid.__contains__("JTT"):
            model = "~/Applications/paml44/dat/jones.dat"
        elif runid.__contains__("WAG"):
            model = "~/Applications/paml44/dat/wag.dat"
        elif runid.__contains__("LG"):
            model = "~/Applications/paml44/dat/lg.dat"
        asrtreepath = get_asr_treepath(DIR, runid)
        asr_commands.append("python ~/Documents/SourceCode/Lazarus/lazarus.py --alignment " + DIR + "/" + gene + "." + DIR_nick[DIR] + ".asr.fasta --tree " + asrtreepath + " --model " + model + " --outputdir " + DIR + "/asr." + runid + " --branch_lengths fixed --asrv 8 --codeml --gapcorrect True --outgroup " + outgroup + " --cleanup True")

#fout = open("asr_commands.sh", "w")
#for a in asr_commands:
#    fout.write(a + "\n")
#fout.close()
#os.system(MPIRUN + " asr_commands.sh")
#exit()

#
# Get ancestors
#
getanc_commands = []
for DIR in DIR_runid_lnl:
    for runid in DIR_runid_lnl[DIR]:
        here = os.getcwd()
        asrmsa = DIR + "/" + gene + "." + DIR_nick[DIR] + ".asr.fasta"
        asrtree = get_asr_treepath(DIR, runid)
        model = "~/Applications/paml44/dat/"
        if runid.__contains__("JTT"):
            model += "jones.dat"
        elif runid.__contains__("WAG"):
            model += "wag.dat"
        elif runid.__contains__("LG"):
            model += "lg.dat"
        for ing in ingroups:
            ingroup = ingroups[ing]
            getanc_commands.append("python ~/Documents/SourceCode/Lazarus/lazarus.py --alignment " + asrmsa + " --tree " + asrtree + " --model " + model + " --outputdir " + here + "/" + DIR + "/asr." + runid + " --outgroup " + outgroup + " --ingroup " + ingroup + " --getanc True")
            getanc_commands.append("mv " + DIR + "/asr." + runid + "/ancestor-ml.dat " + DIR + "/asr." + runid + "/anc." + ing + ".dat")
            getanc_commands.append("mv " + DIR + "/asr." + runid + "/ancestor.out.txt " + DIR + "/asr." + runid + "/anc." + ing + ".txt")

#fout = open("getanc_commands.txt", "w")
#for a in getanc_commands:
#    os.system(a)
#    fout.write(a + "\n")
#fout.close()

#exit()

def get_struct_sites():
    msapath = "MSAPROBS/" + gene + ".msaprobs.asr.phylip"
    seed = "Saccharomyces.cerevisiae.IME2"

    start_motif = START_MOTIF #"YQLI"
    end_motif = END_MOTIF #"MPFF"
    startsite = None
    endsite = None
    fin = open(msapath, "r")
    for l in fin.readlines():
        if l.startswith(seed):
            seq = l.split()[1]
            # find the starting motif                                                                                               
            for i in range(0, seq.__len__()):
                if seq[i] == start_motif[0]:
                    here = ""
                    j = i
                    while here.__len__() < start_motif.__len__() and j < seq.__len__():
                        if seq[j] != "-":
                            here += seq[j]
                        j += 1
                    if here  == start_motif:
                        startsite = i + 1
                        #print here, seq[startsite-1:(startsite-1)+4]                                                               
                        break
            for i in range(i, seq.__len__()):
                if seq[i] == end_motif[0]:
                    here = ""
                    j = i
                    while here.__len__() < end_motif.__len__() and j < seq.__len__():
                        if seq[j] != "-":
                            here += seq[j]
                        j += 1
                    #print here, end_motif                                                                                          
                    if here  == end_motif:
                        endsite = j
                        #print seq[endsite-1:(endsite-1)+3]
                        break
    fin.close()
    return [startsite, endsite]



#
# Compare ancestors
#
compare_commands = []
for pair in comparisons:
    #print pair
    outlines = []
    outlines.append("seed " + ingroup_seed[ pair[0] ])

    msapaths = "msapaths "
    for DIR in DIR_runid_lnl:
        for runid in DIR_runid_lnl[DIR]:
            if 0.01 < DIR_runid_pp[DIR][runid]:
                msapaths += DIR + "/asr." + runid + "/reformatted_alignment.phy "
    outlines.append(msapaths)

    for DIR in DIR_runid_lnl:
        for runid in DIR_runid_lnl[DIR]:
            if 0.01 < DIR_runid_pp[DIR][runid]:
                msapath = DIR + "/asr." + runid + "/reformatted_alignment.phy "
                outlines.append("msaname " + msapath + " " + DIR + "-" + runid)

    for DIR in DIR_runid_lnl:
        for runid in DIR_runid_lnl[DIR]:
            if 0.01 < DIR_runid_pp[DIR][runid]:            
                msapath = DIR + "/asr." + runid + "/reformatted_alignment.phy "
                outlines.append("compare " + DIR + "/asr." + runid + "/anc." + pair[0]+ ".dat " + DIR + "/asr." + runid + "/anc." + pair[1] + ".dat " + DIR + "-" + runid)

    for DIR in DIR_runid_lnl:
        for runid in DIR_runid_lnl[DIR]:
            if 0.01 < DIR_runid_pp[DIR][runid]:
                msapath = DIR + "/asr." + runid + "/reformatted_alignment.phy "
                outlines.append("msaweight " + DIR + "-" + runid + " " + DIR_runid_pp[DIR][runid].__str__())

    specpath = "compare_ancs." + pair[0] + "-" + pair[1] + "config.txt"
    fout = open(specpath, "w")
    for o in outlines:
        fout.write(o + "\n")
    fout.close()

    model = "~/Applications/paml44/dat/lg.dat"
    #model = "/common/paml42/dat/lg.dat"
    
    [startsite,endsite] = get_struct_sites()

    compare_commands.append("python ~/Documents/SourceCode/anccomp/compare_ancs.py --specpath " + specpath + " --modelpath " + model + " --window_sizes 1 --metrics k p hb --runid " + pair[0] + "to" + pair[1] + " --restrict_to_seed True --renumber_sites True --force_bin_width 0.25 --highlight_sites " + startsite.__str__() + "-" + endsite.__str__() + " --pdbtoolsdir /Users/victor/Applications/pdbTools_0.2.1 --pdb_path ../homology_models2/" + pair[1] + ".pdb")
    compare_commands.append("mv /Users/victor/Documents/GENES/ime2/2013-Oct15/pymol_script.k.p /Users/victor/Documents/GENES/ime2/2013-Oct15/" + pair[0] + "to" + pair[1] + ".p")

fout = open("SCRIPTS/compareanc_commands.sh", "w")
for c in compare_commands:
    fout.write(c + "\n")
fout.close()
"""
