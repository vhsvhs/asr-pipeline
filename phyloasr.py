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
    #here = os.path.abspath("./")
    here = os.popen('pwd').read().strip()
    commands = []
    for msa in ap.params["msa_algorithms"]:
        faspath = get_fastapath(msa)
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model) 
            if os.path.exists(here + "/" + msa + "/RAxML_info." + runid): # Remove dirty RAxML data.
                os.system("rm " + here + "/" + msa + "/RAxML_*" + runid)
            command = ap.params["raxml_exe"]
            command += "  -s " + faspath
            command += " -n " + runid
            command += " -w " + here + "/" + msa
            command += " -e 0.001"
            command += " -m " + model
            command += " -p 12345"
            command += " > " + here + "/" + msa + "/catch." + runid + ".txt" 
            commands.append(command)
    p = "SCRIPTS/raxml.commands.sh"
    fout = open(p, "w")
    for c in commands:
        fout.write(c + "\n")
    fout.close()
    return p

def get_mlalpha_pp(ap):
    #
    # Fetch ML resultse from RAxML
    # Calculate min, max likelihoods for each model, 
    # and find the best-fitting model.
    #
    # The ML alpha and the PP for each msa/model
    # are written to ap.params["runid_alpha"]
    # and ap.params["runid_pp"]
    #
    ap.params["runid_alpha"] = {}
    ap.params["runid_pp"] = {}
    for msa in ap.params["msa_algorithms"]:
        runid_lnl = {} 
        maxl = None
        minl = None
        suml = 0
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model)
            lpath = get_raxml_logpath(msa, model)
            if False == os.path.exists(lpath):
                print "Something is wrong. I can't find the log file from RAxML:", lpath
                exit()
            fin = open(lpath, "r")
            lines = fin.readlines()
            fin.close()
            lastline = lines[ lines.__len__()-1 ]
            lnl = float(lastline.split()[1])
            runid_lnl[runid] = lnl
            if maxl == None:
                maxl = lnl
            if minl == None:
                minl = lnl
            if lnl > maxl:
                maxl = lnl
            if lnl < minl:
                minl = lnl

            lpath = get_raxml_infopath(msa, model)
            if False == os.path.exists(lpath):
                print "Something is wrong. I can't find the info file from RAxML:", lpath
                exit()
            fin = open(lpath, "r")
            for l in fin.xreadlines():
                if l.startswith("alpha[0]:"):
                    l = l.strip()
                    tokens = l.split()
                    this_alpha = float(tokens[1])
                    ap.params["runid_alpha"][runid] = this_alpha
            fin.close()
    
        for runid in runid_lnl:
            suml += (runid_lnl[runid] - minl)
            
        fout = open(msa + "/raxml.lnl.summary.txt", "w")
        for runid in runid_lnl:
            lnl = runid_lnl[runid]
            if runid_lnl.keys().__len__() <= 1:
                pp = 1.0
            else:
                pp = (lnl-minl)/suml
            ap.params["runid_pp"][runid] = pp
            alpha = ap.params["runid_alpha"][runid]
            special = ""
            if lnl == maxl:
                special = " (ML) "
            line = runid + "\t" + lnl.__str__() + "\t" + "%.4f"%pp + "\t%.4f"%alpha + special
            fout.write(line + "\n")
        fout.close()
        
        
def calc_alrt(ap):
    alrt_commands = []
    for msa in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model)
            mltreepath = get_raxml_treepath(msa, runid)
            phylippath = get_phylippath(msa)
            
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
                gammastr = " --nclasses 8 --alpha " + ap.params["runid_alpha"][runid].__str__()
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
            command += " > " + msa + "/catch.phyml." + model + ".txt"
            alrt_commands.append( command )
    fout = open("SCRIPTS/alrt_commands.sh", "w")
    for c in alrt_commands:
        fout.write(c + "\n")
    fout.close()
    return "SCRIPTS/alrt_commands.sh"     


def calc_alr(ap):
    """Convert aLRTs to aLRs"""
    alr_commands = []
    for msa in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            alrt_treepath = get_alrt_treepath(msa, model)
            alr_treepath = get_alr_treepath(msa, model)
            
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
        mltreepath = get_raxml_treepath(DIR, runid)
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
            
            fastapath = get_asr_fastapath(msa)
            
            fin = open(fastapath, "r")
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
        
            if False == os.path.exists(msa + "/asr." + model):
                os.system("mkdir " + msa + "/asr." + model)
            modelstr = get_model_path(model, ap)
            asrtreepath = get_raxml_treepath(msa, runid)
            asr_commands.append(ap.params["lazarus_exe"] + " --alignment " + fastapath + " --tree " + asrtreepath + " --model " + modelstr + " --outputdir " + msa + "/asr." + model + " --branch_lengths fixed --asrv 8 --codeml --gapcorrect True --outgroup " + ap.params["outgroup"] + " --cleanup True")

    fout = open("SCRIPTS/asr_commands.sh", "w")
    for a in asr_commands:
        fout.write(a + "\n")
    fout.close()
    return "SCRIPTS/asr_commands.sh"


#
# Get ancestors
#
def get_getanc_commands(ap):
    getanc_commands = []
    for msa in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model)
            here = os.getcwd()
            asrmsa = get_asr_fastapath(msa)
            asrtree = get_raxml_treepath(msa, runid)
            modelstr = ap.params["mmfolder"]
            if runid.__contains__("JTT"):
                modelstr += "/jones.dat"
            elif runid.__contains__("WAG"):
                modelstr += "/wag.dat"
            elif runid.__contains__("LG"):
                modelstr += "/lg.dat"
            for ing in ap.params["ingroup"]:
                getanc_commands.append(ap.params["lazarus_exe"] + " --alignment " + asrmsa + " --tree " + asrtree + " --model " + modelstr + " --outputdir " + here + "/" + msa + "/asr." + model + " --outgroup " + ap.params["outgroup"] + " --ingroup " + ap.params["ingroup"][ing] + " --getanc True")
                getanc_commands.append("mv " + msa + "/asr." + model + "/ancestor-ml.dat " + msa + "/asr." + model + "/" + ing + ".dat")
                getanc_commands.append("mv " + msa + "/asr." + model + "/ancestor.out.txt " + msa + "/asr." + model + "/" + ing + ".txt")
    
    fout = open("SCRIPTS/getanc_commands.txt", "w")
    for a in getanc_commands:
        os.system(a)
        fout.write(a + "\n")
    fout.close()
    return "SCRIPTS/getanc_commands.txt"

#exit()

#ap.params["compareanc"]
def get_compareanc_commands(ap):
    compare_commands = []
    for pair in ap.params["compareanc"]:    
        #msapathlines = "msapaths "
        msanamelines = ""
        comparelines = ""
        weightlines = ""
        for msa in ap.params["msa_algorithms"]:
            for model in ap.params["raxml_models"]:
                runid = get_runid(msa, model)
                pp = ap.params["runid_pp"][runid]
                if pp > 0.0:
                    msapath = msa + "/asr." + model + "/reformatted_alignment.phy"
                    
                    #msapathlines += msapath + " "
                    msanamelines += "msaname " + msapath + " " + runid + "\n"
                    comparelines += "compare " + msa + "/asr." + model + "/" + pair[0] + ".dat " + msa + "/asr." + model + "/" + pair[1] + ".dat " + runid + "\n"
                    weightlines += "msaweight " + runid + " " + pp.__str__() + "\n"
            
        specpath = "compare_ancs." + pair[0] + "-" + pair[1] + ".config.txt"
        fout = open(specpath, "w")
        fout.write("seed " + ap.params["seedtaxa"][ pair[0] ] + "\n")
        #fout.write(msapathlines)
        fout.write(msanamelines)
        fout.write(comparelines)
        fout.write(weightlines)
        fout.close()
    
        modelstr = get_model_path(model, ap)
        
        #
        # At this point,the call to get_bound.. is unnecessary
        # because we already trimmed the alignment to start/stop motif
        # boundaries earlier, in the trim_alignment method
        #
        #[startsite,endsite] = get_boundary_sites(  get_phylipfull_path("msaprobs",ap))
            
        c = ap.params["anccomp"]
        c += " --specpath " + specpath
        c += " --modelpath " + modelstr
        c += " --window_sizes 1"
        c += " --metrics k p hb"
        c += " --runid " + pair[0] + "to" + pair[1]
        c += " --restrict_to_seed True"
        c += " --renumber_sites True"
        #c += " --force_bin_width 0.25"
        #if startsite != None and endsite != None:
        #    c += " --highlight_sites " + startsite.__str__() + "-" + endsite.__str__()
        compare_commands.append(c)
    
    fout = open("SCRIPTS/compareanc_commands.sh", "w")
    for c in compare_commands:
        fout.write(c + "\n")
    fout.close()
    return "SCRIPTS/compareanc_commands.sh"

