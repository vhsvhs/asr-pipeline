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
import dendropy
from dendropy import Tree
from asrpipelinedb_api import *


#
# Create RAxML commands for each alignment and each model.
#
def write_raxml_commands(con):
    #here = os.path.abspath("./")
    here = os.popen('pwd').read().strip()
    commands = []
    for msa in get_alignment_method_names(con):
        phypath = get_raxml_phylippath(msa)
        for model in get_phylo_modelnames(con):
            runid = get_runid(msa, model) 
            """Remove residual files from previous RAxML job with the same runID"""
            if os.path.exists(here + "/" + msa + "/RAxML_info." + runid): # Remove dirty RAxML data.
                rmcmd = "rm " + here + "/" + msa + "/RAxML_*"
                os.system(rmcmd)
            """Verify the raxml executable"""
            cv = get_setting_values(con, "raxml_exe")
            if cv == None:
                write_error(con, "I cannot find the executable path for raxml.")
                exit()
            command = cv[0] # i.e., raxml_exe
            command += " -s " + phypath
            command += " -n " + runid
            command += " -w " + here + "/" + msa
            command += " -e 0.001"
            command += " -m " + model
            command += " -p 12345"
            constraint_tree = get_setting_values(con, "constraint_tree")
            if constraint_tree != None:
                command += " -g " + constraint_tree[0]
            command += " > " + here + "/" + msa + "/catch." + runid + ".txt" 
            commands.append(command)
    p = "SCRIPTS/raxml.commands.sh"
    fout = open(p, "w")
    for c in commands:
        fout.write(c + "\n")
    fout.close()
    return p


def make_fasttree_command(con, ap, outdir, phylippath):
    newtree = get_fasttree_path(phylippath)
    
    c = get_setting_values(con, "fasttree_exe")[0]
    c += " -wag "
    c += " < " + phylippath + " > " + newtree
    return c
    
    
def check_raxml_output(con):   
    """Checks RAXmL's output, and imports trees into UnsupportedMlPhylogenies"""
    cur = con.cursor()
    
    sql = "delete from UnsupportedMlPhylogenies"
    cur.execute(sql)
    con.commit()
    
    here = os.popen('pwd').read().strip()
    commands = []
    for msa in get_alignment_method_names(con):
        sql = "select id from AlignmentMethods where name='" + msa + "'"
        cur.execute(sql)
        msaid = cur.fetchone()[0]
        
        phypath = get_raxml_phylippath(msa)
        for model in get_phylo_modelnames(con):
            runid = get_runid(msa, model) 
            raxml_treepath = get_raxml_treepath(msa, runid)
            if False == os.path.exists(raxml_treepath):
                print "I can't find the expected result from RAxML at " + here + "/" + msa + "/RAxML_bestTree." + runid
                write_error(con, "I can't find the expected result from RAxML at " + here + "/" + msa + "/RAxML_bestTree." + runid)
                
                sql = "delete from PhyloModels where name='" + model + "'"
                cur.execute(sql)
                con.commit()
                exit()
            
            # get the modelid
            sql = "select modelid from PhyloModels where name='" + model + "'"
            cur.execute(sql)
            modelid = cur.fetchone()[0]
                        
            # get the Newick string from the tree
            fin = open(raxml_treepath, "r")
            treestring = fin.readline()
            fin.close()
            
            """Re-root the tree based on the user-supplied outgroup"""
            rooted_treestring = reroot_newick(con, treestring)
            
            """Continue here -- check this -- new for December 2014"""
            rooted_treestring = re.sub("'", "", rooted_treestring)
                    
            sql = "insert into UnsupportedMlPhylogenies (almethod,phylomodelid,newick) VALUES("
            sql += msaid.__str__() + "," + modelid.__str__() + ",'" + rooted_treestring + "')"
            cur.execute(sql)
            con.commit()


def get_mlalpha_pp(con):
    """Fetch ML resultse from RAxML
    Calculate min, max likelihoods for each model, 
    and find the best-fitting model.
    
    The ML alpha and the PP for each msa/model
    are written to runid_alpha
    and runid_pp"""
    
    cur = con.cursor()
    sql = "delete from TreeMl"
    cur.execute(sql)
    con.commit()
    sql = "delete from TreeAlpha"
    cur.execute(sql)
    con.commit()
    sql = "delete from TreePP"
    cur.execute(sql)
    con.commit()
    

    for msa in get_alignment_method_names(con):
        sql = "select id from AlignmentMethods where name='" + msa + "'"
        cur.execute(sql)
        msaid = cur.fetchone()[0]

        runid_alpha = {}        
        runid_lnl = {}
        treeid_lnl = {}
        treeid_pp = {} 
        maxl = None # maximum log-L
        minl = None # minimum log-L
        suml = 0    # sum of log-L for all ML trees from this alignment
        for model in get_phylo_modelnames(con):
            sql = "select modelid from PhyloModels where name='" + model + "'"
            cur.execute(sql)
            modelid = cur.fetchone()[0]

            sql = "select id from UnsupportedMlPhylogenies where almethod=" + msaid.__str__() + " and phylomodelid=" + modelid.__str__()
            cur.execute(sql)
            treeid = cur.fetchone()[0]
            
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
            
            sql = "insert into TreeMl(mltreeid, likelihood) values("
            sql += treeid.__str__() + "," + lnl.__str__() + ")"
            cur.execute(sql)
            con.commit()
                        
            runid_lnl[runid] = lnl
            treeid_lnl[treeid] = lnl
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
                    runid_alpha[runid] = this_alpha
            fin.close()
            
            sql = "insert into TreeAlpha(mltreeid, alpha) values("
            sql += treeid.__str__() + "," + this_alpha.__str__() + ")"
            cur.execute(sql)
            con.commit()
    
        for runid in runid_lnl:
            suml += (runid_lnl[runid] - minl)
        
        """Compute Posterior Probability of each Tree"""
        for treeid in treeid_lnl:
            lnl = treeid_lnl[treeid]
            if treeid_lnl.keys().__len__() <= 1:
                treeid_pp[treeid] = 1.0
            else:
                treeid_pp[treeid] = (lnl-minl)/suml
                
            sql = "insert into TreePP (mltreeid, pp) values(" + treeid.__str__() + "," + treeid_pp[treeid].__str__() + ")"
            cur.execute(sql)
            con.commit() 
        
        """Write the results to raxml.lnl.summary.txt"""
        fout = open(msa + "/raxml.lnl.summary.txt", "w")
        for runid in runid_lnl:
            lnl = runid_lnl[runid]
            if runid_lnl.keys().__len__() <= 1:
                pp = 1.0
            else:
                pp = (lnl-minl)/suml

            alpha = runid_alpha[runid]
            special = ""
            if lnl == maxl:
                special = " (ML) "
            line = runid + "\t" + lnl.__str__() + "\t" + "%.4f"%pp + "\t%.4f"%alpha + special
            fout.write(line + "\n")
        fout.close()
                
  
def calc_alrt(con):
    cur = con.cursor()
    alrt_commands = []
    for msa in get_alignment_method_names(con):
    
        sql = "select id from AlignmentMethods where name='" + msa + "'"
        cur.execute(sql)
        msaid = cur.fetchone()[0]
        
        for model in get_phylo_modelnames(con):
            runid = get_runid(msa, model)
            mltreepath = get_raxml_treepath(msa, runid)
            phylippath = get_raxml_phylippath(msa)
                        
            sql = "select modelid from PhyloModels where name='" + model + "'"
            cur.execute(sql)
            modelid = cur.fetchone()[0]
            
            sql = "select id from UnsupportedMlPhylogenies where almethod=" + msaid.__str__() + " and phylomodelid=" + modelid.__str__()
            cur.execute(sql)
            treeid = cur.fetchone()[0]
            
            if False == os.path.exists(mltreepath):
                print "Something is wrong. I can't find the ML tree output from RAxML:", mltreepath
                exit()
            if False == os.path.exists(phylippath):
                print "Something is wrong. I can't find the trimmed Phylip alignment", phylippath
                exit()
            
            modelstr = ""
            if model.__contains__("JTT"):
                modelstr = "JTT"
            elif model.__contains__("WAG"):
                modelstr = "WAG"
            elif model.__contains__("LG"):
                modelstr = "LG"
            elif model.__contains__("CPREV"):
                modelstr = "CpREV"
            elif model.__contains__("RTREV"):
                modelstr = "RtREV"
            elif model.__contains__("DAYHOFF"):
                modelstr = "Dayhoff"
            gammastr = ""
            if model.__contains__("GAMMA"):
                """Get the alpha value for the model."""
                sql = "select alpha from TreeAlpha where mltreeid=" + treeid.__str__()
                cur.execute(sql)
                alpha = cur.fetchone()[0]
                gammastr = " --nclasses 8 --alpha " + alpha.__str__()
            
            x = get_setting_values(con, "phyml_exe")
            if x == None:
                write_error(con, "I cannot find the path to the PhyML in your settings. Error 298")
                exit()
            
            command = x[0]
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


def calc_alr(con):
    """Convert aLRTs to aLRs"""
    alr_commands = []
    for msa in get_alignment_method_names(con):
        for model in get_phylo_modelnames(con):
            alrt_treepath = get_alrt_treepath(msa, model)
            alr_treepath = get_alr_treepath(msa, model)
            
            if False == os.path.exists(alrt_treepath):
                print "Something is wrong. I can't find the ML tree with aLRT values at", alrt_treepath
                exit()
            
            alrt_to_alr( alrt_treepath, alr_treepath )

            if False == os.path.exists(alr_treepath):
                print "Something is wrong. I can't find converted ML tree with aLR values at", alr_treepath
                exit()


def import_supported_trees(con):
    cur = con.cursor()
    sql = "delete from BranchSupportMethods"
    cur.execute(sql)
    sql = "delete from SupportedMlPhylogenies"
    cur.execute(sql)
    con.commit()
    
    sql = "insert into BranchSupportMethods (name) values('aLRT')"
    cur.execute(sql)
    con.commit()
    sql = "select id from BranchSupportMethods where name='aLRT'"
    cur.execute(sql)
    alrt_id = cur.fetchone()[0]

    sql = "insert into BranchSupportMethods (name) values('aLR')"
    cur.execute(sql)
    con.commit()
    sql = "select id from BranchSupportMethods where name='aLR'"
    cur.execute(sql)
    alr_id = cur.fetchone()[0]
    
    sql = "select id, almethod, phylomodelid from UnsupportedMlPhylogenies"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        treeid = ii[0]
        msaid = ii[1]
        modelid = ii[2]
    
        sql = "select name from AlignmentMethods where id=" + msaid.__str__()
        cur.execute(sql)
        msaname = cur.fetchone()[0]
        
        sql = "select name from PhyloModels where modelid=" + modelid.__str__()
        cur.execute(sql)
        modelname = cur.fetchone()[0]
        
        alrt_treepath = get_alrt_treepath(msaname, modelname)
        alr_treepath = get_alr_treepath(msaname, modelname)
        if False == os.path.exists(alrt_treepath):
            write_error(con, "Oops, I cannot find the aLRT-supported tree at " + alrt_treepath)
            exit()
        if False == os.path.exists(alr_treepath):
            write_error(con, "Oops, I cannot find the aLR-supported tree at " + alr_treepath)
            exit()       
        alrt_string = open(alrt_treepath,"r").readline()
        alr_string = open(alr_treepath,"r").readline()
                
        """Ensure the trees are rooted according go the outgroup."""
        alrt_string = reroot_newick(con, alrt_string)
        alrt_string = reroot_newick(con, alr_string)
        
        """Save the aLRT tree"""
        sql = "insert into SupportedMlPhylogenies (unsupportedmltreeid,newick,supportmethodid)"
        sql += " values(" + treeid.__str__() + ",\"" + alrt_string + "\"," + alrt_id.__str__() + ")"
        cur.execute(sql)
        con.commit()
        
        """Save the aLR tree"""
        sql = "insert into SupportedMlPhylogenies (unsupportedmltreeid,newick,supportmethodid)"
        sql += " values(" + treeid.__str__() + ",\"" + alr_string + "\"," + alr_id.__str__() + ")" 
        cur.execute(sql)
        con.commit()
            

def get_asr_commands(con):
    cur = con.cursor()
    
    asr_commands = []
    for msa in get_alignment_method_names(con):
        for model in get_phylo_modelnames(con):
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
            modelstr = get_model_path(model, con)
            asrtreepath = get_raxml_treepath(msa, runid)
            lazaarus_exe = get_setting_values(con, "lazarus_exe")[0]
            
            sql = "select id from TaxaGroups where name='outgroup'"
            cur.execute(sql)
            outgroup_id = cur.fetchone()[0]
            #print "456:", outgroup_id
            outgroup_list = get_taxaid_in_group(con, outgroup_id)
            #print "457:", outgroup_list
            outgroup_list = [get_taxon_name(con, i) for i in outgroup_list]
            #print "458:", outgroup_list
            outgroup_string = "[" + ",".join( outgroup_list ) + "]"
            #print "460:", outgroup_string
            asr_commands.append(lazaarus_exe + " --alignment " + fastapath + " --tree " + asrtreepath + " --model " + modelstr + " --outputdir " + msa + "/asr." + model + " --branch_lengths fixed --asrv 8 --codeml --gapcorrect True --outgroup " + outgroup_string)

    fout = open("SCRIPTS/asr_commands.sh", "w")
    for a in asr_commands:
        fout.write(a + "\n")
    fout.close()
    return "SCRIPTS/asr_commands.sh"

def check_asr_output(con):
    """Verifies that ASR occurred correctly, and if it did, imports the results into the database.
        This method fills the tables Ancestors, AncestralStates, and AncestralCladogram"""
    cur = con.cursor()
    
    sql = "delete from Ancestors"
    cur.execute(sql)
    sql = "delete from AncestralStates"
    cur.execute(sql)
    sql = "delete from AncestralCladogram"
    cur.execute(sql)
    con.commit()

    
    for msa in get_alignment_method_names(con):
        sql = "select id from AlignmentMethods where name='" + msa + "'"
        cur.execute(sql)
        msaid = cur.fetchone()[0]
        
        for model in get_phylo_modelnames(con):
            sql = "select modelid from PhyloModels where name='" + model + "'"
            cur.execute(sql)
            modelid = cur.fetchone()[0]
            
            """What is the reference tree ID for this msa-model pair?"""
            sql = "select id from UnsupportedMlPhylogenies where almethod=" + msaid.__str__() + " and phylomodelid=" + modelid.__str__()
            cur.execute(sql)
            treeid = cur.fetchone()[0]
            
            """Where does the ASR output live?"""
            outputdir = msa + "/asr." + model
            
            """Get the post-ASR cladogram with ancestral node numbers."""
            cladopath = outputdir + "/tree1/tree1.txt" 
            if False == os.path.exists( cladopath ):
                write_error(con, "I cannot find the expected cladogram at " + cladopath)
                exit()
            fin = open(cladopath, "r")
            cladostring = fin.readlines()[3]
            fin.close()
            
            """Ensure the cladogram is rooted according to the outgroup."""
            cladostring = reroot_newick(con, cladostring)
            cladostring = re.sub("\n", "", cladostring)
            
            """Save the cladogram to the database."""
            sql = "insert into AncestralCladogram (unsupportedmltreeid,newick) values(" + treeid.__str__() + ",\"" + cladostring + "\")"
            #print "528", sql
            cur.execute(sql)
            con.commit()
            
            """Where do the ancestral DAT files live? (i.e. the state posterior probability distributions for each ancestor)"""
            datpath = outputdir + "/tree1"
             
            if False == os.path.exists(datpath):
                write_error(con, "I cannot find reconstructed ancestors in the folder " + datpath)
                exit()
            
            """Parse each ancestral DAT file."""
            for d in os.listdir(datpath):
                if False == d.endswith(".dat"):
                    continue
                nodenum = int( re.sub("node", "", d.split(".")[0]) )
                ssp = get_site_state_pp( datpath + "/" + d)
                
                sql = "insert into Ancestors (almethod, phylomodel, name) values("
                sql += msaid.__str__() + "," + modelid.__str__() + ",'Node" + nodenum.__str__() + "')"
                cur.execute(sql) 
                con.commit()
                
                sql = "select id from Ancestors where almethod=" + msaid.__str__() + " and phylomodel=" + modelid.__str__() + " and name='Node" + nodenum.__str__() + "'"
                cur.execute(sql)
                ancid = cur.fetchone()[0]
                
                for site in ssp:
                    for state in ssp[site]:
                        pp = ssp[site][state]
                        sql = "insert into AncestralStates (ancid, site, state, pp) values("
                        sql += ancid.__str__() + "," + site.__str__() + ",'" + state + "',"
                        sql += pp.__str__() + ")"
                        cur.execute(sql)
                con.commit()
                
                
def get_getanc_commands(con):
    """Returns the path to a script containing the commands for comparing ancestors."""
    cur = con.cursor()
    getanc_commands = []
    for msa in get_alignment_method_names(con):
        for model in get_phylo_modelnames(con):
            runid = get_runid(msa, model)
            here = os.getcwd()
            asrmsa = get_asr_fastapath(msa)
            asrtree = get_raxml_treepath(msa, runid)
            modelstr = get_setting_values(con, "mmfolder")[0]
            
            if runid.__contains__("JTT"):
                modelstr += "/jones.dat"
            elif runid.__contains__("WAG"):
                modelstr += "/wag.dat"
            elif runid.__contains__("LG"):
                modelstr += "/lg.dat"

            sql = "select id from TaxaGroups where name='outgroup'"
            cur.execute(sql)
            outgroupid = cur.fetchone()[0]
            taxa = get_taxaid_in_group(con, outgroupid)
            taxa = [get_taxon_name(con,i) for i in taxa]
            outgroup_string = "[" + ",".join(taxa) + "]"
                
            ingroup_ids = get_ingroup_ids(con)  
            for ing in ingroup_ids:
                sql = "select name from TaxaGroups where id=" + ing.__str__()
                cur.execute(sql)
                ingroup_name = cur.fetchone()[0]
                taxa = get_taxaid_in_group(con, ing)
                taxa = [get_taxon_name(con,i) for i in taxa]
                ingroup_string = "[" + ",".join(taxa) + "]"
                
                lazarus_exe = get_setting_values(con, "lazarus_exe")[0]
                getanc_commands.append(lazarus_exe + " --alignment " + asrmsa + " --tree " + asrtree + " --model " + modelstr + " --outputdir " + here + "/" + msa + "/asr." + model + " --outgroup " + outgroup_string + " --ingroup " + ingroup_string + " --getanc True")
                getanc_commands.append("mv " + msa + "/asr." + model + "/ancestor-ml.dat " + msa + "/asr." + model + "/" + ingroup_name + ".dat")
                getanc_commands.append("mv " + msa + "/asr." + model + "/ancestor.out.txt " + msa + "/asr." + model + "/" + ingroup_name + ".txt")
    
    fout = open("SCRIPTS/getanc_commands.txt", "w")
    for a in getanc_commands:
        run_subprocess(a)
        fout.write(a + "\n")
    fout.close()
    return "SCRIPTS/getanc_commands.txt"


def check_getanc_output(con):
    """This method fills the tables AncestorsAlias, AncestorsGroups"""
    
    cur = con.cursor()
    
    sql = "delete from AncestorsAlias"
    cur.execute(sql)
    sql = "delete from AncestorsGroups"
    cur.execute(sql)
    con.commit()
    
    sql = "select id, name from TaxaGroups"
    cur.execute(sql)
    x = cur.fetchall()
    groupid_name = {}
    for ii in x:
        if ii[1] != "outgroup":
            groupid_name[ ii[0] ] = ii[1]
    
    sql = "select id from TaxaGroups where name='outgroup'"
    cur.execute(sql)
    outgroupid = cur.fetchone()[0]
    
    for msa in get_alignment_method_names(con):
        sql = "select id from AlignmentMethods where name='" + msa + "'"
        cur.execute(sql)
        msaid = cur.fetchone()[0]
        
        for model in get_phylo_modelnames(con):
            sql = "select modelid from PhyloModels where name='" + model + "'"
            cur.execute(sql)
            modelid = cur.fetchone()[0]
            
            """Check the Lazarus job status log file for errors."""
            logpath = msa + "/asr." + model + "/lazarus_job_status.log"
            if os.path.exists(logpath):
                fin = open(logpath, "r")
                line1 = fin.readline()
                if line1.__contains__("error") or line1.__contains__("Error"):
                    write_error(con, "There was an error extracting ancestors:")
                    write_error(con, line1 )
                    exit()
            
            """For each named ancestor, read its special text file
                and determine its node number."""
            for ingroupid in groupid_name:
                ancpath = msa + "/asr." + model + "/" + groupid_name[ ingroupid ] + ".txt"
                if False == os.path.exists( ancpath ):
                    write_error(con, "I cannot find the expected ancestor " + ancpath)
                    exit()
                
                """this_node will be ID of the node corresponding to the ancestor."""
                this_node = None
                fin = open(ancpath, "r")
                for l in fin.readlines():
                    if l.__contains__("On the ML tree"):
                        this_node = int(  re.sub("#", "", l.split()[8])  )
                fin.close()
                                
                sql = "select id from Ancestors where name='Node" + this_node.__str__() + "' and almethod=" + msaid.__str__() + " and phylomodel=" + modelid.__str__()
                cur.execute(sql)
                x = cur.fetchall()
                if x.__len__() == 0:
                    write_error(con, "Error 596 - I cannot find the entry in Ancestors for the ingroup ID " + ingroupid.__str__())
                    exit()
                ancid = x[0][0]
                
                sql = "insert into AncestorsAlias (ancid, alias) values(" + ancid.__str__() + ",'" + groupid_name[ ingroupid ] + "')"
                sql
                cur.execute(sql)
                con.commit()

                sql = "insert into AncestorsGroups (ancid, ingroupid, outgroupid) values(" + ancid.__str__() + "," + ingroupid.__str__() + "," + outgroupid.__str__() + ")"
                cur.execute(sql)
                con.commit()


def setup_pdb_maps(ap):
    """This method is experimental"""
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


def get_compareanc_commands(con):
    cur = con.cursor()
    
    sql = "delete from FScore_Tests"
    cur.execute(sql)
    con.commit()
    
    compare_commands = []
    
    """pairs is a list of tuples, each containing alias names for ancestors."""    
    pairs = get_ancestral_comparison_pairs(con)
        
    for pair in pairs:   
        #msapathlines = "msapaths "
        msanamelines = ""
        comparelines = ""
        weightlines = ""
        for msa in get_alignment_method_names(con):
            
            sql = "select id from AlignmentMethods where name='" + msa + "'"
            cur.execute(sql)
            msaid = cur.fetchone()[0]
            
            for model in get_phylo_modelnames(con):
                runid = get_runid(msa, model)
                
                sql = "select modelid from PhyloModels where name='" + model + "'"
                cur.execute(sql)
                modelid = cur.fetchone()[0]

                ancid1 = get_ancestorid(con, pair[0], msaid, modelid)
                ancid2 = get_ancestorid(con, pair[1], msaid, modelid)
                
                if ancid1 == ancid2:
                    """Potentially skip this analysis if anc1 == anc2"""
                    write_log(con, "I'm skipping the Fscore analysis for " + msa + "." + model + " " + pair.__str__() + " because the ancestors are the same node.")
                    continue

                sql = "select id from UnsupportedMlPhylogenies where almethod=" + msaid.__str__() + " and phylomodelid=" + modelid.__str__()
                cur.execute(sql)
                treeid = cur.fetchone()[0]

                sql = "select pp from TreePP where mltreeid in (select id from UnsupportedMlPhylogenies where almethod=" + msaid.__str__() + " and phylomodelid=" + modelid.__str__() + ")"
                cur.execute(sql)
                pp = cur.fetchone()[0]
            
                msapath = msa + "/asr." + model + "/reformatted_alignment.phy"

                msanamelines += "msaname " + msapath + " " + runid + "\n"
                comparelines += "compare " + msa + "/asr." + model + "/" + pair[0] + ".dat " + msa + "/asr." + model + "/" + pair[1] + ".dat " + runid + "\n"
                weightlines += "msaweight " + runid + " " + pp.__str__() + "\n"
                                
                sql = "insert into FScore_Tests (almethod, phylomodel, ancid1, ancid2)"
                sql += " values(" + msaid.__str__() + "," + modelid.__str__() + "," + ancid1.__str__() + "," + ancid2.__str__() + ")"
                cur.execute(sql)
                con.commit()
                
        specpath = "compare_ancs." + pair[0] + "-" + pair[1] + ".config.txt"
        if msanamelines.__len__() > 2 and comparelines.__len__() > 2 and weightlines.__len__() > 2:      
            fout = open(specpath, "w")
            
            seed_taxonname = get_setting_values(con, "seedtaxa")[0]
            fout.write("seed " + seed_taxonname + "\n")
            #if pair[1] in ap.params["map2pdb"]: # there was a definition to map scores on this PDB:
            #    fout.write("pdb " + ap.params["map2pdb"][pair[1]] + "\n")
            fout.write(msanamelines)
            fout.write(comparelines)
            fout.write(weightlines)
            fout.close()
    
            #
            # to-do continue here -- fix this so we're not using just one model.
            #
            modelstr = get_model_path(   get_phylo_modelnames(con)[0] , con)
                        
            c = get_setting_values(con, "anccomp_exe")[0]
            c += " --specpath " + specpath
            c += " --modelpath " + modelstr
            c += " --window_sizes 1"
            c += " --metrics k p Df"
            c += " --runid " + pair[0] + "to" + pair[1]
            c += " --restrict_to_seed True"
            c += " --renumber_sites True"
            #if ap.params["do_pdb_analysis"]:
            #    c += " --pdbtoolsdir " + ap.params["pdbtoolsdir"]
            #    c += " --pymol_exe " + get_setting_values(con, "pymol_exe")
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

def parse_compareanc_results(con):
    cur = con.cursor()
    
    sql = "delete from FScore_Sites"
    cur.execute(sql)
    con.commit()
        
    pairs = get_ancestral_comparison_pairs(con)    
    for pair in pairs:
        outdir = pair[0] + "to" + pair[1]
    
        for msa in get_alignment_method_names(con):        
            sql = "select id from AlignmentMethods where name='" + msa + "'"
            cur.execute(sql)
            msaid = cur.fetchone()[0]
                
            for modelname in get_phylo_modelnames(con):                
                sql = "select modelid from PhyloModels where name='" + modelname + "'"
                cur.execute(sql)
                modelid = cur.fetchone()[0]
                
                ancid1 = get_ancestorid(con, pair[0], msaid, modelid)
                ancid2 = get_ancestorid(con, pair[1], msaid, modelid)
                
                if ancid1 == ancid2:
                    """Skip this analysis, because it probably wasn't computed."""
                    continue
                
                sql = "select id from FScore_Tests where almethod=" + msaid.__str__() + " and phylomodel=" + modelid.__str__()
                sql += " and ancid1=" + ancid1.__str__() + " and ancid2=" + ancid2.__str__()
                cur.execute(sql)
                x = cur.fetchall()
                if x.__len__() == 0:
                    continue
                testid = x[0][0]
                
                site_df = {}
                summary_path = outdir + "/Df." + msa + "." + modelname + ".summary.txt"
                if False == os.path.exists(summary_path):
                    write_error(con, "I cannot find your Df output file at " + summary_path)
                    continue
                fin = open(summary_path, "r")
                lines = fin.readlines()
                for l in lines:
                    if l.__len__() > 2:
                        tokens = l.split()
                        site = int(tokens[0])
                        site_df[site] = float(tokens[3])
                fin.close()

                site_p = {}
                summary_path = outdir + "/p." + msa + "." + modelname + ".summary.txt"
                if False == os.path.exists(summary_path):
                    write_error(con, "I cannot find your p output file at " + summary_path)
                    continue
                fin = open(summary_path, "r")
                lines = fin.readlines()
                for l in lines:
                    if l.__len__() > 2:
                        tokens = l.split()
                        site = int(tokens[0])
                        site_p[site] = float(tokens[3])
                fin.close()

                site_k = {}
                summary_path = outdir + "/k." + msa + "." + modelname + ".summary.txt"
                if False == os.path.exists(summary_path):
                    write_error(con, "I cannot find your k output file at " + summary_path)
                    continue
                fin = open(summary_path, "r")
                lines = fin.readlines()
                for l in lines:
                    if l.__len__() > 2:
                        tokens = l.split()
                        site = int(tokens[0])
                        site_k[site] = float(tokens[3])                
                fin.close()
                
                for site in site_df:
                    df = site_df[site]
                    p = site_p[site]
                    k = site_k[site]
                    sql = "insert into FScore_Sites (testid,site,df,k,p)"
                    sql += " values(" + testid.__str__() + "," + site.__str__() + "," + df.__str__() + "," + k.__str__() + "," + p.__str__() + ")"
                    cur.execute(sql)

                con.commit()
                
                sql = "select count(*) from FScore_Sites where testid=" + testid.__str__()
                cur.execute(sql)
                count = cur.fetchone()[0]
                
                write_log(con, "I found " + count.__str__() + " F-scored sites for test ID " + testid.__str__() )                
        
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
    t.read_from_string(newick.__str__(), "newick")
    fin.close()
    return t.length()

def match_ancestors_across_models(con):
    """This method fills data in the table AncestorsAcrossModels"""
    cur = con.cursor()
    
    modelids = get_phylo_modelids(con)
    msaids = get_alignment_method_ids(con)
    
    ancid_childrenids = {} # key = Ancestor ID, value = list of Taxa IDs
    
    """Pull the map of taxon names to IDs from the database.
        We'll access this information a lot, so let's save it in a separate hashtable
        rather than repeatedly querying the databse."""
    taxonname_id = {}
    sql = "select id, shortname from Taxa"
    cur.execute(sql)
    for ii in cur.fetchall():
        id = ii[0]
        name = ii[1]
        taxonname_id[ name ] = id
    
    for modelid in modelids:
        for msaid in msaids:
            sql = "select newick from AncestralCladogram where unsupportedmltreeid in (select id from UnsupportedMlPhylogenies where almethod=" + msaid.__str__() + " and phylomodelid=" + modelid.__str__() + ")"
            cur.execute(sql)
            xx = cur.fetchone()
            if xx == None:
                write_error(con, "I cannot find the ancestral Newick cladogram for almethod=" + msaid.__str__() + " and phylomodelid=" + modelid.__str__())
            cladonewick = xx[0].__str__()
            
            t = Tree()
            t.read_from_string(cladonewick, "newick")
            
            for node in t.nodes():
                if node.is_leaf() == False and node.level() > 0:                    
                    sql = "select id from Ancestors where name='Node" + node.label + "' and almethod=" + msaid.__str__() + " and phylomodel=" + modelid.__str__()
                    print "976:", modelid, msaid, node.label
                    cur.execute(sql)
                    ancid = cur.fetchone()[0]
                    ancid_childrenids[ancid] = []
                    
                    for l in node.leaf_iter():
                        #print "978:", l
                        taxonname = l.as_newick_string()
                        #print "980:", taxonname
                        taxonname = re.sub("'", "", taxonname)
                        ancid_childrenids[ancid].append( taxonname_id[taxonname] )
    
    ancid_matches = {} # key = Ancestor ID, value = list of other ancestor IDs with the same children.
    for anc1 in ancid_childrenids:
        ancid_matches[anc1] = []
        mychildren = ancid_childrenids[anc1]
        mychildren.sort()
        for anc2 in ancid_childrenids:
            if anc1 == anc2:
                """Skip the self comparison."""
                continue
            theirchildren = ancid_childrenids[anc2]
            theirchildren.sort()
            if mychildren == theirchildren:
                ancid_matches[anc1].append(anc2)
    
    sql = "delete from AncestorsAcrossModels"
    cur.execute(sql)
    con.commit()
    
    for anc1 in ancid_matches:
        for anc2 in ancid_matches[anc1]:
            sql = "insert into AncestorsAcrossModels (ancid, same_ancid) values(" + anc1.__str__() + "," + anc2.__str__() + ")"
            cur.execute(sql)
    con.commit()
                        
                        
            
    