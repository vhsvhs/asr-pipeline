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

from log import *
from tools import *
from alrt2alr import *
from dendropy import Tree


#
# Create RAxML commands for each alignment and each model.
#
def write_raxml_commands(ap):
    #here = os.path.abspath("./")
    here = os.popen('pwd').read().strip()
    commands = []
    for msa in ap.params["msa_algorithms"]:
        phypath = get_trimmed_phylippath(msa)
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model) 
            if os.path.exists(here + "/" + msa + "/RAxML_info." + runid): # Remove dirty RAxML data.
                rmcmd = "rm " + here + "/" + msa + "/RAxML_*"
                print rmcmd
                #p = run_subprocess(rmcmd)
                os.system(rmcmd)
            command = ap.params["raxml_exe"]
            command += " -s " + phypath
            command += " -n " + runid
            command += " -w " + here + "/" + msa
            command += " -e 0.001"
            command += " -m " + model
            command += " -p 12345"
            if ap.params["constraint_tree"] != None:
                command += " -g " + ap.params["constraint_tree"]
            command += " > " + here + "/" + msa + "/catch." + runid + ".txt" 
            commands.append(command)
    p = "SCRIPTS/raxml.commands.sh"
    fout = open(p, "w")
    for c in commands:
        fout.write(c + "\n")
    fout.close()
    return p

"""Depricated."""
def make_raxml_quick_command(con, ap, outdir, phylippath, runid):
    """This is intended primarily for the quick RAxML runs called within the ZORRO analysis."""
    here = os.popen('pwd').read().strip()
    if os.path.exists(here + "/" + outdir + "/RAxML_info." + runid): # Remove dirty RAxML data.
        rmcmd = "rm " + here + "/" + outdir + "/RAxML_*"
        print rmcmd
        os.system(rmcmd)
    command = ap.params["raxml_exe"]
    command += " -s " + phylippath
    command += " -n " + runid
    command += " -w " + here + "/" + outdir
    command += " -e 0.1"
    command += " -m PROTGAMMALG -c 1"
    command += " -p 12345"
    command += " -x 12345 -N 100 -f a"
    if ap.params["constraint_tree"] != None:
        command += " -g " + ap.params["constraint_tree"]
    command += " > " + here + "/" + outdir + "/catch." + runid + ".txt" 
    return command

def make_fasttree_command(con, ap, outdir, phylippath):
    newtree = get_fasttree_path(phylippath)
    
    c = ap.params["fasttree_exe"]
    c += " -wag "
    c += " < " + phylippath + " > " + newtree
    return c
    
    
def check_raxml_output(ap):
    here = os.popen('pwd').read().strip()
    commands = []
    for msa in ap.params["msa_algorithms"]:
        phypath = get_trimmed_phylippath(msa)
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model) 
            raxml_treepath = get_raxml_treepath(msa, runid)
            if False == os.path.exists(raxml_treepath):
                print "I can't find the expected result from RAxML at " + here + "/" + msa + "/RAxML_bestTree." + runid
                write_error(ap, "I can't find the expected result from RAxML at " + here + "/" + msa + "/RAxML_bestTree." + runid)
                exit()


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

def read_lnl_summary(ap):
    if "runid_alpha" not in ap.params:
        ap.params["runid_alpha"] = {}
    if "runid_pp" not in ap.params:
        ap.params["runid_pp"] = {}
    
    for msa in ap.params["msa_algorithms"]:
        lnlpath = msa + "/raxml.lnl.summary.txt"
        fin = open(lnlpath, "r")
        for l in fin.xreadlines():
            if l.__len__() > 2:
                tokens = l.split()
                runid = tokens[0]
                lnl = float( tokens[1] )
                pp = float( tokens[2] )
                alpha = float( tokens[3] )
                ap.params["runid_alpha"][runid] = alpha
                ap.params["runid_pp"][runid] = pp 
        fin.close()
        
def calc_alrt(ap):
    alrt_commands = []
    for msa in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model)
            mltreepath = get_raxml_treepath(msa, runid)
            phylippath = get_trimmed_phylippath(msa)
            
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
            elif runid.__contains__("CPREV"):
                modelstr = "CpREV"
            elif runid.__contains__("RTREV"):
                modelstr = "RtREV"
            elif runid.__contains__("DAYHOFF"):
                modelstr = "Dayhoff"
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
                run_subprocess("mkdir " + msa + "/asr." + model)
            modelstr = get_model_path(model, ap)
            asrtreepath = get_raxml_treepath(msa, runid)
            asr_commands.append(ap.params["lazarus_exe"] + " --alignment " + fastapath + " --tree " + asrtreepath + " --model " + modelstr + " --outputdir " + msa + "/asr." + model + " --branch_lengths fixed --asrv 8 --codeml --gapcorrect True --outgroup " + ap.params["outgroup"])

    fout = open("SCRIPTS/asr_commands.sh", "w")
    for a in asr_commands:
        fout.write(a + "\n")
    fout.close()
    return "SCRIPTS/asr_commands.sh"

#
# Get ancestors
#ancseqs
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
        run_subprocess(a)
        fout.write(a + "\n")
    fout.close()
    return "SCRIPTS/getanc_commands.txt"

#exit()

def check_getanc_output(ap):
    # Returns (bool, message)
    #
    # Read the Lazarus job status log file, look for errors.
    #
    for msa in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            logpath = msa + "/asr." + model + "/lazarus_job_status.log"
            if os.path.exists(logpath):
                fin = open(logpath, "r")
                line1 = fin.readline()
                if line1.__contains__("error") or line1.__contains__("Error"):
                    return (False, line1)
    return (True, None)

def setup_pdb_maps(ap):
    """The goal of this method is to build ap.params["map2pdb"][ anc ] from the PHYRE output folder."""
        
    if "phyre_out" not in ap.params:
        return
        #print "\n. I can't find any PHYRE output."
    if False == os.path.exists(os.getcwd() + "/" + ap.params["phyre_out"]):
        print "\n. I can't find your PHYRE output folder at", os.getcwd() + "/" + ap.params["phyre_out"]
        print "\n. Check your configuration file and try again.\n"
        exit()
    if False == os.path.exists(os.getcwd() + "/" + ap.params["phyre_out"] + "/summaryInfo"):
        print "\n. I can't find the file 'summaryInfo' in your PHYRE output folder at", ap.params["phyre_out"]
        print "\n. Did your PHYRE output get corrupted or deleted?n"
        exit()
        
    fin = open( os.getcwd() + "/" + ap.params["phyre_out"] + "/summaryInfo", "r" )
    lines = fin.readlines()[1:]
    for l in lines:
        if l.__len__() > 2:
            tokens = l.split()
            datpath = tokens[0]
            this_msa = datpath.split("/")[0]
            this_model = datpath.split("/")[1].split(".")[1]
            this_anc = datpath.split("/")[2].split(".")[0] 
            
            pdbpath = os.getcwd() + "/" + ap.params["phyre_out"] + "/" + tokens[2] + ".final.pdb"
            
            ap.params["map2pdb"][ this_anc ] = pdbpath
    fin.close()

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

        fout.write("seed " + ap.params["seedtaxa"][ pair[1] ] + "\n")
        if pair[1] in ap.params["map2pdb"]: # there was a definition to map scores on this PDB:
            fout.write("pdb " + ap.params["map2pdb"][pair[1]] + "\n")
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
        c += " --metrics k p Df"
        c += " --runid " + pair[0] + "to" + pair[1]
        c += " --restrict_to_seed True"
        c += " --renumber_sites True"
        if ap.params["do_pdb_analysis"]:
            c += " --pdbtoolsdir " + ap.params["pdbtoolsdir"]
            c += " --pymol_exe " + ap.params["pymol_exe"]
        #c += " --force_bin_width 0.25"
        #if startsite != None and endsite != None:
        #    c += " --highlight_sites " + startsite.__str__() + "-" + endsite.__str__()
        compare_commands.append(c)
        compare_commands.append("source run_rscripts.sh")
    
    fout = open("SCRIPTS/compareanc_commands.sh", "w")
    for c in compare_commands:
        fout.write(c + "\n")
    fout.close()
    return "SCRIPTS/compareanc_commands.sh"

def get_branch_supports( treepath ):
    supports = []
    fin = open(treepath, "r")
    for line in fin.readlines():
        tokens = line.split(":")
        for t in tokens:
            if t.__contains__(")") and False == t.__contains__(";"):
                ts = t.split(")")
                if ts[1] != '':
                    support = float(ts[1])
                    supports.append( support )
    fin.close()
    return supports

def get_sum_of_branches( treepath ):
    fin = open(treepath, "r")
    newick = fin.readline().strip()
    t = Tree()
    t.read_from_string(newick, "newick")
    fin.close()
    return t.length()