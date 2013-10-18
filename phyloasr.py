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
from alrt2alr import *


#
# Create RAxML commands for each alignment and each model.
#
def write_raxml_commands(ap):
    here = os.getcwd()
    commands = []
    for msa in ap.params["msa_algorithms"]:
        phypath = get_phylip_path(msa, ap)
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model) 
            command = ap.params["raxml_exe"]
            command += "  -s " + phypath
            command += " -n " + runid
            command += " -w " + here + "/OUT." + msa
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

def post_raxml(ap):
    #
    # Fetch ML results,
    # Calculate min, max likelihoods for each model, 
    # and find the best-fitting model.
    #
    runid_alpha = {}
    for msa in ap.params["msa_algorithms"]:
        runid_lnl = {} 
        maxl = None
        minl = None
        suml = 0
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model)
            spath = get_statspath(msa, model)
            if False == os.path.exists(spath):
                print "Something is wrong. I can't find the output from RAxML:", spath
                exit()
            fin = open(spath, "r")
            lines = fin.readlines()
            for l in lines:
                if l.startswith("Final GAMMA-based Score of best tree"):
                    tokens = l.split()
                    lnl = float(tokens[6])
                    runid_lnl[runid] = lnl
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
                    runid_alpha[runid] = this_alpha
            fin.close()
    
        for runid in runid_lnl:
            suml += (runid_lnl[runid] - minl)
            
        runid_pp = {}
        fout = open("OUT." + msa + "/raxml.lnl.summary.txt", "w")
        for runid in runid_lnl:
            lnl = runid_lnl[runid]
            pp = (lnl-minl)/suml
            runid_pp[runid] = pp
            special = ""
            if lnl == maxl:
                special = " (ML) "
            line = runid + "\t" + lnl.__str__() + "\t" + "%.4f"%pp + "\t" + special
            fout.write(line + "\n")
        fout.close()


    alrt_commands = []
    for msa in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model)
            mltreepath = get_ml_treepath(msa, runid)
            phylippath = get_phylip_path(msa, ap)
            
            if False == os.path.exists(mltreepath):
                print "Something is wrong. I can't find the ML tree output from RAxML:", mltreepath
                exit()
            if False == os.path.exists(phylippath):
                print "Something is wrong. I can't find the trimmed Phylip alignment", phylippath
                exit()
            
            modelstr = ""
            if runid.__contains__("JTT"):
                modelstr = "JTT"
            elif runid.__contains__("WAG"):
                modelstr = "WAG"
            elif runid.__contains__("LG"):
                modelstr = "LG"
            gammastr = ""
            if runid.__contains__("GAMMA"):
                gammastr = " --nclasses 8 --alpha " + runid_alpha[runid].__str__()
            command = ap.params["phyml_exe"]
            command += " --input " + phylippath
            command += " -u " + mltreepath
            command += " --model " + modelstr
            command += gammastr
            command += " -f m"  # use fixed pi values
            command += " -o n"  # don't optimize anything
            command += " -b -1" # calculate ALRTs
            command += " --run_id " + model + ".alrt"
            command += " --datatype aa"
            command += " > OUT." + msa + "/catch.phyml." + model + ".txt"
            alrt_commands.append( command )
    fout = open("SCRIPTS/alrt_commands.sh", "w")
    for c in alrt_commands:
        fout.write(c + "\n")
    fout.close()
    return "SCRIPTS/alrt_commands.sh"     

def calc_alr(ap):
    #
    # Convert ALRTs to ALRs
    #
    alr_commands = []
    for msa in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model)
            alrt_treepath = get_alrt_treepath(msa, model, ap)
            alr_treepath = get_alr_treepath(msa, model, ap)
            
            if False == os.path.exists(alrt_treepath):
                print "Something is wrong. I can't find the ML tree with aLRT values at", alrt_treepath
                exit()
            
            alrt_to_alr( alrt_treepath, alr_treepath )

            if False == os.path.exists(alr_treepath):
                print "Something is wrong. I can't find converted ML tree with aLR values at", alr_treepath
                exit()

"""

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
"""

def get_asr_commands(ap):
    asr_commands = []
    for msa in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model)
            
            fastapath = get_fasta_path(msa, ap)
            
            fin = open(get_fasta_path(msa, ap), "r")
            lines = fin.readlines()
            fin.close()
            innames = [] # array of taxa in the trimmed alignment:
            for i in range(0, lines.__len__()):
                if lines[i].startswith(">"):
                    l = lines[i].strip()
                    name = re.sub(">", "", l)
                    if name in innames:
                        print ". Your alignment" + fastapath + " has two copies of", name
                        exit()
                    innames.append(name)
        
            fin = open(fastapath, "r")
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
        
            if False == os.path.exists("OUT." + msa + "/asr." + model):
                os.system("mkdir OUT." + msa + "/asr." + model)
            modelstr = "~/Applications/paml44/dat/lg.dat"
            if runid.__contains__("JTT"):
                modelstr = "~/Applications/paml44/dat/jones.dat"
            elif runid.__contains__("WAG"):
                modelstr = "~/Applications/paml44/dat/wag.dat"
            elif runid.__contains__("LG"):
                modelstr = "~/Applications/paml44/dat/lg.dat"
            asrtreepath = get_asr_treepath(msa, runid)
            asr_commands.append(ap.params["lazarus_exe"] + " --alignment " + fastapath + " --tree " + asrtreepath + " --model " + modelstr + " --outputdir OUT." + msa + "/asr." + model + " --branch_lengths fixed --asrv 8 --codeml --gapcorrect True --outgroup " + ap.params["outgroup"] + " --cleanup True")

    fout = open("SCRIPTS/asr_commands.sh", "w")
    for a in asr_commands:
        fout.write(a + "\n")
    fout.close()
    return "SCRIPTS/asr_commands.sh"
    #os.system(MPIRUN + " asr_commands.sh")
    #exit()

"""
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
