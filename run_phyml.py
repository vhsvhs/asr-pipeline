#
# run_phyml.py
#
# Version 1.0
#

import os, re, sys


###################################
gene = "cmgc"
outgroup = "[Homo.sapiens.DYR1A]"
ingroups = {}
ingroups["cdkime2"] = "[Saccharomyces.cerevisiae.CDK1,Saccharomyces.cerevisiae]"
ingroups["mok-ime2"] = "[Saccharomyces.cerevisiae,Homo.sapiens.MOK]"
ingroups["mak-ime2"] = "[Saccharomyces.cerevisiae,Homo.sapiens.MAK]"
ingroups["ime2-cneo"] = "[Saccharomyces.cerevisiae,Cryptococcus.neoformans]"
ingroups["ime2-spom"] = "[Saccharomyces.cerevisiae,Schizosaccharomyces.pombe]"
ingroups["ime2-ncras"] = "[Saccharomyces.cerevisiae,Neurospora.crassa]"
ingroups["ime2-ylip"] = "[Saccharomyces.cerevisiae,Yarrowia.lipolytica]"

comparisons = []
comparisons.append( ["cdkime2","mok-ime2"] )
comparisons.append( ["mok-ime2","mak-ime2"] )
comparisons.append( ["mak-ime2","ime2-cneo"] )
comparisons.append( ["ime2-cneo","ime2-spom"] )
comparisons.append( ["ime2-spom","ime2-ncras"] )
comparisons.append( ["ime2-ncras","ime2-ylip"] )

seed = "Saccharomyces.cerevisiae"

PHYML = "/common/REPOSITORY/phyml-20120412/src/phyml "
CLEAN_TREES = "python ~/Documents/EclipseWork/PhyloWidgets/clean_trees.py "
ALR2 = "python ~/Documents/EclipseWork/PhyloWidgets/alrt2alr.py "
SEP = "."

#DIRS = ["MUSCLE","PRANK", "MSAPROBS"]
DIRS = ["MSAPROBS", "MUSCLE"]
#DIRS = ["MUSCLE"]
DIR_nick = {}
DIR_nick["MUSCLE"] = "muscle"
DIR_nick["PRANK"] = "prank"
DIR_nick["MSAPROBS"] = "msaprobs"

models = ["LG", "WAG", "JTT"]
#models = ["LG"]
ng = [1,8]
nf = ["e", "m"]
nb = []

######################################

commands = []
DIR_runids = {}
DIR_runid_statspath = {}

def get_statspath(DIR, runid):
    keyword = DIR_nick[DIR]
    statspath = DIR + "/" + gene + SEP + keyword + ".phylip_phyml_stats_" + runid + ".txt"
    return statspath

def get_treepath(DIR, runid):
    keyword = DIR_nick[DIR]
    if runid.__contains__("-g"):
        return DIR + "/" + gene + SEP + keyword + ".phylip_phyml_tree_" + runid + ".txt"
    else:
        return DIR + "/" + gene + SEP + keyword + ".phylip_phyml_mixtrees_" + runid + ".txt"

def get_cleantreepath(DIR, runid):
    keyword = DIR_nick[DIR]
    return DIR + "/" + gene + SEP + keyword + ".phylip_phyml_tree_" + runid + ".clean.tre"

def get_alrtreepath(DIR, runid):
    keyword = DIR_nick[DIR]
    return DIR + "/" + gene + SEP + keyword + ".phylip_phyml_tree_" + runid + ".alr.tre"

def get_asrtreepath(DIR, runid):
    keyword = DIR_nick[DIR]
    return DIR + "/" + gene + SEP + keyword + ".phylip_phyml_tree_" + runid + ".asr.tre"

for DIR in DIRS:
    DIR_runid_statspath[DIR] = {}
    DIR_runids[DIR] = []
    keyword = DIR_nick[DIR]
    for model in models:
        for f in nf:
            for g in ng:
                runid = model + "-g" + g.__str__() + "-f" + f
                DIR_runids[DIR].append(runid)
                statspath = get_statspath( DIR, runid )
                DIR_runid_statspath[DIR][runid] = statspath
                command = PHYML + "  --input " + DIR + "/" + gene + SEP + keyword + ".phylip"
                command += " --run_id " + runid
                command += " --datatype aa "
                command += " --nclasses " + g.__str__()
                command += " -f " + f
                command += " -b -1"
                command += " --model " + model
                command += " -o tlr"
                #command += " --search NNI "
                command += " > " + DIR + "/catch." + runid + ".txt" 
                commands.append(command)
            for b in nb:
                runid = model + "-b" + b.__str__() + "-f" + f
                DIR_runids[DIR].append(runid)
                statspath = get_statspath( DIR, runid )
                DIR_runid_statspath[DIR][runid] = statspath
                command = PHYML + "  --input " + DIR + "/" + gene + SEP + keyword + ".phylip"
                command += " --run_id " + runid
                command += " --datatype aa "
                command += " --nclasses 1 "
                command += " --bl_mixtures " + b.__str__()
                command += " -f " + f
                command += " -b -1"
                command += " --model " + model
                command += " -o tlr"
                #command += " --search NNI "
                command += " > " + DIR + "/catch." + runid + ".txt" 
                commands.append(command)
fout = open("phyml_commands.txt", "w")
for c in commands:
    fout.write( c + "\n" )
fout.close() 
#exit()
os.system("rm MSAPROBS/*phyml*")
os.system("rm MUSCLE/*phyml*")
os.system("rm PRANK/*phyml*")
os.system("mpirun -np 29 --machinefile hosts.txt /common/bin/mpi_dispatch phyml_commands.txt")
#exit()

DIR_runid_lnl = {}
DIR_suml = {}
DIR_runid_pp = {}
for DIR in DIR_runids:
    DIR_runid_lnl[DIR] = {}
    DIR_runid_pp[DIR] = {}
    maxl = None
    minl = None
    for runid in DIR_runids[DIR]:
        statspath = DIR_runid_statspath[DIR][runid]
        if False == os.path.exists(statspath):
            print "I can't find:", statspath
            continue
        fin = open(statspath, "r")
        lines = fin.readlines()
        for l in lines:
            if l.startswith(". Log-likelihood:"):
                tokens = l.split()
                lnl = float(tokens[2])
                DIR_runid_lnl[DIR][runid] = lnl
                if maxl == None:
                    maxl = lnl
                if minl == None:
                    minl = lnl
                if lnl > maxl:
                    maxl = lnl
                if lnl < minl:
                    minl = lnl
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

for DIR in DIR_runids:
    for runid in DIR_runid_lnl[DIR]:
        treepath = get_treepath(DIR, runid)
        print treepath
        cleantreepath = get_cleantreepath(DIR, runid)
        asrtreepath = get_asrtreepath(DIR, runid)
        alrtreepath = get_alrtreepath(DIR, runid)
        #os.system(CLEAN_TREES + " " + treepath + " " + cleantreepath)
        #os.system(CLEAN_TREES + " " + treepath + " " + asrtreepath + " --strip_support True")
        #os.system(ALR2 + " " + cleantreepath + " > " + alrtreepath )

#exit()

#
# Run ASR
#
asr_commands = []
for DIR in DIR_runid_lnl:
    fin = open(DIR + "/" + gene + "." + DIR_nick[DIR] + ".fasta", "r")
    lines = fin.readlines()
    fin.close()
    innames = [] # array of taxa in the trimmed alignment:
    for i in range(0, lines.__len__()):
        if lines[i].startswith(">"):
            l = lines[i].strip()
            name = re.sub(">", "", l)
            if name in innames:
                print ". Your alignment has two copies of", name
                exit()
            print "182:", name
            innames.append(name)

    fin = open(DIR + "/" + gene + "." + DIR_nick[DIR] + ".full.fasta", "r")
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

    #os.system("perl ~/Applications/seqConverter.pl -d" + DIR + "/" + gene + "." + DIR_nick[DIR] + ".asr.fasta -ope")

    for runid in DIR_runid_lnl[DIR]:
        os.system("mkdir " + DIR + "/asr." + runid)
        asrv = "1"
        if runid.__contains__("g8"):
            asrv = "8"
        model = "/common/REPOSITORY/paml44/dat/lg.dat"
        if runid.__contains__("JTT"):
            model = "/common/REPOSITORY/paml44/dat/jones.dat"
        if runid.__contains__("WAG"):
            model = "/common/REPOSITORY/paml44/dat/wag.dat"
        asrtreepath = get_asrtreepath(DIR, runid)
        asr_commands.append("python ~/Documents/EclipseWork/Lazarus/lazarus.py --alignment " + DIR + "/" + gene + "." + DIR_nick[DIR] + ".asr.fasta --tree " + asrtreepath + " --model " + model + " --outputdir " + DIR + "/asr." + runid + " --branch_lengths fixed --asrv " + asrv + " --codeml --gapcorrect True --outgroup " + outgroup + " --cleanup True")

fout = open("asr_commands.sh", "w")
for a in asr_commands:
    fout.write(a + "\n")
fout.close()
os.system("mpirun -np 29 --machinefile hosts.txt /common/bin/mpi_dispatch asr_commands.sh")
#for a in asr_commands:
#    os.system(a)
#exit()

#
# Get ancestors
#
#
getanc_commands = []
for DIR in DIR_runid_lnl:
    for runid in DIR_runid_lnl[DIR]:
        here = os.getcwd()
        asrmsa = DIR + "/" + gene + "." + DIR_nick[DIR] + ".asr.fasta"
        asrtree = get_asrtreepath(DIR, runid)
        model = "/common/paml42/dat/lg.dat"
        if runid.__contains__("JTT"):
            model = "/common/paml42/dat/jones.dat"
        if runid.__contains__("WAG"):
            model = "/common/paml42/dat/wag.dat"
        for ing in ingroups:
            ingroup = ingroups[ing]
            getanc_commands.append("python ~/Documents/EclipseWork/Lazarus/lazarus.py --alignment " + asrmsa + " --tree " + asrtree + " --model " + model + " --outputdir " + here + "/" + DIR + "/asr." + runid + " --outgroup " + outgroup + " --ingroup " + ingroup + " --getanc True")
            getanc_commands.append("mv " + DIR + "/asr." + runid + "/ancestor-ml.dat " + DIR + "/asr." + runid + "/anc." + ing + ".dat")
            getanc_commands.append("mv " + DIR + "/asr." + runid + "/ancestor.out.txt " + DIR + "/asr." + runid + "/anc." + ing + ".txt")

fout = open("getanc_commands.txt", "w")
for a in getanc_commands:
    os.system(a)
    fout.write(a + "\n")
fout.close()

#
# Compare ancestors
#
#
compare_commands = []
for pair in comparisons:
    #print pair
    outlines = []
    outlines.append("seed " + seed)

    msapaths = "msapaths "
    for DIR in DIR_runid_lnl:
        for runid in DIR_runid_lnl[DIR]:
            msapaths += DIR + "/asr." + runid + "/reformatted_alignment.phy "
        outlines.append(msapaths)

    for DIR in DIR_runid_lnl:
        for runid in DIR_runid_lnl[DIR]:
            msapath = DIR + "/asr." + runid + "/reformatted_alignment.phy "
            outlines.append("msaname " + msapath + " " + DIR + "-" + runid)

    for DIR in DIR_runid_lnl:
        for runid in DIR_runid_lnl[DIR]:
            msapath = DIR + "/asr." + runid + "/reformatted_alignment.phy "
            outlines.append("compare " + DIR + "/asr." + runid + "/anc." + pair[0]+ ".dat " + DIR + "/asr." + runid + "/anc." + pair[1] + ".dat " + DIR + "-" + runid)

    for DIR in DIR_runid_lnl:
        for runid in DIR_runid_lnl[DIR]:
            msapath = DIR + "/asr." + runid + "/reformatted_alignment.phy "
            outlines.append("msaweight " + DIR + "-" + runid + " " + DIR_runid_pp[DIR][runid].__str__())

    specpath = "compare_ancs." + pair[0] + "-" + pair[1] + "config.txt"
    fout = open(specpath, "w")
    for o in outlines:
        fout.write(o + "\n")
    fout.close()

    #model = "/common/REPOSITORY/paml44/dat/lg.dat"
    model = "/common/paml42/dat/lg.dat"

    compare_commands.append("python ~/Documents/EclipseWork/anccomp/compare_ancs.py --specpath " + specpath + " --modelpath " + model + " --window_sizes 1 5 --metrics h hb --runid " + pair[1] + " --renumber_sites True" )

fout = open("compareanc_commands.sh", "w")
for c in compare_commands:
    fout.write(c + "\n")
fout.close()
