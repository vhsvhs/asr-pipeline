import re,sys,os
from tools import *
from log import *
from asrpipelinedb import *
from asrpipelinedb_api import *

def read_config_file(con, ap):
    """This method reads the user-written configuration file.
    The format of the file is generally KEYWORD = VALUE,
    where KEYWORD can be a limited set of parameters, such as GENEID or RAXML,
    and VALUE is the corresponding parameter value.
    
    For more information about valid KEYWORDs and VALUEs, see the user documentation,
    or simply read this method and look for lines containing if tokens[0].startswith("...")
    
    """
    cpath = ap.getArg("--configpath")
    if False == os.path.exists("./" + cpath):
        print "ERROR: I can't find your configfile at", cpath
        write_error(con, "I can't find your configfile at " + cpath)
        exit()

    cur = con.cursor()
    sql = "delete from Settings"
    cur.execute(sql)
    sql = "delete from PhyloSoftwares"
    cur.execute(sql)
    con.commit()
    
    fin = open(cpath, "r")

    ap.params["end_motif"] = None
    ap.params["start_motif"] = None

    for l in fin.xreadlines():
        l = l.strip()
        if l.startswith("#"):
            continue
        if l.__len__() < 2:
            continue
        tokens = l.split("=")
        if tokens.__len__() < 1: 
            continue
        
        if tokens[0].startswith("GENE_ID"):
            ap.params["geneid"] = re.sub(" ", "", tokens[1])
        
        if tokens[0].startswith("PROJECT_TITLE"):
            ap.params["project_title"] = tokens[1]

        if tokens[0].startswith("SEQUENCES"):
            ergseqpath = re.sub(" ", "", tokens[1])
            ap.params["ergseqpath"] = ergseqpath
            add_setting_value(con, "ergseqpath", ergseqpath)
        
        elif tokens[0].startswith("RAXML"):
            exe = tokens[1].strip() 
            ap.params["raxml_exe"] = exe
            add_setting_value(con, "raxml_exe", exe)
            
            sql = "select count(*) from PhyloSoftwares where name='raxml'"
            cur.execute(sql)
            count = cur.fetchone()[0]
            if count == 0:
                sql = "insert into PhyloSoftwares (name, exe_path) VALUES('raxml','" + exe + "')"
                cur.execute(sql)
                con.commit()
        
        elif tokens[0].startswith("FASTTREE"):
            exe = tokens[1].strip() 
            ap.params["fasttree_exe"] = exe
            add_setting_value(con, "fasttree_exe", exe)
        
        elif tokens[0].startswith("PHYML"):
            exe = tokens[1].strip() 
            ap.params["phyml_exe"] = exe
            add_setting_value(con, "phyml_exe", exe)
            
            sql = "select count(*) from PhyloSoftwares where name='phyml_exe'"
            cur.execute(sql)
            count = cur.fetchone()[0]
            if count == 0:
                sql = "insert into PhyloSoftwares (name, exe_path) VALUES('phyml','" + exe + "')"
                cur.execute(sql)
                con.commit()

        elif tokens[0].startswith("LAZARUS"):
            exe = tokens[1].strip() 
            ap.params["lazarus_exe"] = exe
            add_setting_value(con, "lazarus_exe", exe)

        elif tokens[0].startswith("MARKOV_MODEL_FOLDER"):
            ap.params["mmfolder"] = tokens[1].strip()
        
        elif tokens[0].startswith("MPIRUN"):
            exe = tokens[1].strip() 
            ap.params["mpirun_exe"] = exe
            add_setting_value(con, "mpirun_exe", exe)
        
        elif tokens[0].startswith("RUN"):
            ap.params["run_exe"] = tokens[1]
        
        elif tokens[0].startswith("MSAPROBS"):
            exe = tokens[1].strip()
            ap.params["msaprobs_exe"] = exe
            add_setting_value(con, "msaprobs_exe", exe)

        elif tokens[0].startswith("MUSCLE"):
            exe = tokens[1].strip()
            ap.params["muscle_exe"] = exe
            add_setting_value(con, "muscle_exe", exe)

        elif tokens[0].startswith("PRANK"):
            exe = tokens[1].strip()
            ap.params["prank_exe"] = exe
            add_setting_value(con, "prank_exe", exe)

        elif tokens[0].startswith("MAFFT"):
            exe = tokens[1].strip()
            ap.params["mafft_exe"] = exe
            add_setting_value(con, "mafft_exe", exe)
        
        elif tokens[0].startswith("ANCCOMP"):
            ap.params["anccomp"] = tokens[1].strip()
        
        elif tokens[0].startswith("PYMOL"):
            exe = tokens[1].strip()
            ap.params["pymol_exe"] = exe
            add_setting_value(con, "pymol_exe", exe)

        elif tokens[0].startswith("ZORRO"):
            exe = tokens[1].strip()
            ap.params["zorro_exe"] = exe
            add_setting_value(con, "zorro_exe", exe)

        elif tokens[0].startswith("ALIGNMENT_ALGORITHMS"):
            x = tokens[1].split()
            ap.params["msa_algorithms"] = []
            for i in x:
                ap.params["msa_algorithms"].append( i )

        elif tokens[0].startswith("THRESHOLDS_ZORRO"):
            x = tokens[1].split()
            ap.params["zorro_thresholds"] = []
            for i in x:
                ap.params["zorro_thresholds"].append( float(i) )
            if 1.0 not in ap.params["zorro_thresholds"]:
                ap.params["zorro_thresholds"].append( 1.0 )
        
        elif tokens[0].startswith("MODELS_RAXML"):
            x = tokens[1].split()
            for i in x:
                sql = "select count(*) from PhyloModels where name='" + i + "'"
                cur.execute(sql)
                count = cur.fetchone()[0]
                if count == 0:               
                    sql = "insert into PhyloModels(name) VALUES('" + i + "')"
                    cur.execute(sql)
                    con.commit()
        
        elif tokens[0].startswith("START_MOTIF"):
            ap.params["start_motif"] = re.sub(" ", "", tokens[1])
            if ap.params["start_motif"].__contains__("None"):
                ap.params["start_motif"] = None 
        
        elif tokens[0].startswith("END_MOTIF"):
            ap.params["end_motif"] = re.sub(" ", "", tokens[1])
            if ap.params["end_motif"].__contains__("None"):
                ap.params["end_motif"] = None 
 
        elif tokens[0].startswith("CONSTRAINT_TREE"):
            constraint_tree = re.sub(" ", "", tokens[1])
            if False == os.path.exists(constraint_tree):
                print "\n. I can't find your constraint tree at", constraint_tree
                exit()
            add_setting_value(con, "constraint_tree", constraint_tree, unique=True)

        elif tokens[0].startswith("FAMILY_DESCRIPTION"):
            ap.params["family_description"] = tokens[1] 
 
# depricated: 
# continue here... and look in the msa get_bondary_sites function      
        elif tokens[0].startswith("SEED"):
            cleaned_taxa_name = re.sub(" ", "", tokens[1])
            ap.params["seed_motif_seq"] = cleaned_taxa_name
            add_setting_value(con, "seedtaxa", cleaned_taxa_name)
        
        elif tokens[0].startswith("N_BAYES_SAMPLES"):
            ap.params["n_bayes_samples"] = int(tokens[1])
        
        elif tokens[0].startswith("OUTGROUP"):
            """Note: If you change code here, also check the method 'verify_erg_seqs'
            in which this data gets written into the SQL DB."""
            ap.params["outgroup"] = re.sub(" ", "", tokens[1])

        
        elif tokens[0].startswith("INGROUP"):
            """Note: If you change code here, also check the method 'verify_erg_seqs'
            in which this data gets written into the SQL DB."""
            if "ingroup" not in ap.params:
                ap.params["ingroup"] = {}
            anc = tokens[0].split()[1]
            ingroup = tokens[0].split()[2:]
            ingroup = "".join(ingroup)
            ingroup = re.sub(" ", "", ingroup) # tolerate and remove whitespace.
            ap.params["ingroup"][ anc ] = ingroup
            
            
            
        elif tokens[0].startswith("ASRSEED"):
            if "seedtaxa" not in ap.params:
                ap.params["seedtaxa"] = {}
            ts = tokens[0].split()
            anc = ts[1]
            seed = ts[2]
            ap.params["seedtaxa"][ anc ] = re.sub(" ", "", seed)

        elif tokens[0].startswith("COMPARE"):
            anc1 = tokens[0].split()[1]
            anc2 = tokens[0].split()[2]
            if "compareanc" not in ap.params:
                ap.params["compareanc"] = []
            ap.params["compareanc"].append( (anc1,anc2) )
            
        elif tokens[0].startswith("USE_MPI"):
            print tokens[1]
            answer = re.sub(" ", "", tokens[1])
            if answer == "on" or answer == "True":
                ap.params["usempi"] = True
            else:
                ap.params["usempi"] = False
        
        elif tokens[0].startswith("PHYRE_OUTPUT"):
            ap.params["phyre_out"] = re.sub("\ ", "", tokens[1])
            ap.params["phyre_out"].strip()

        elif tokens[0].startswith("PDBTOOLSDIR"):
            ap.params["pdbtoolsdir"] = re.sub("\ ", "", tokens[1])
            ap.params["pdbtoolsdir"].strip()

        elif tokens[0].startswith("MAP2PDB"): # map Df scores for this ancestor onto the PDB path
            anc = tokens[0].split()[1]
            pdbpath = tokens[0].split()[2]
            if "map2pdb" not in ap.params:
                ap.params["map2pdb"] = {}
            if anc not in ap.params["map2pdb"]:
                ap.params["map2pdb"][anc] = []
            if pdbpath not in ap.params["map2pdb"][anc]:
                ap.params["map2pdb"][ anc ].append( pdbpath )
        
        elif tokens[0].startswith("HTML_SPECIAL1"):
            thing = ""
            for i in range(1, tokens.__len__()):
                t = tokens[i]
                if i < tokens.__len__()-1:
                    thing += t + "="
                else:
                    thing +=  t
            ap.params["HTML_SPECIAL1"] = thing
    
    fin.close()

def verify_config(con, ap):
    """Will return nothing if the configuration is ok.  Will error and quit if the
    configuration is flawed."""
    
    cur = con.cursor()
    
    x = get_setting_values(con, "ergseqpath")
    if x == None:
        write_error(con, "You did not specify SEQUENCES, a path to the original sequences, in your configuration file.")
        exit()
    
    if "run_exe" not in ap.params:
        ap.params["run_exe"] = "source"
    
    if "ancestors" in ap.params:
        for a in ap.params["ancestors"]:
            print a, ap.params["seedtaxa"]
            if a not in ap.params["seedtaxa"]:
                print "\n. ERROR: You did not specify a SEED for the ancestor", a
                write_error(con, "You did not specify a SEED for the ancestor " + a)
                exit()
        for c in ap.params["compareanc"]:
            a1 = c[0]
            a2 = c[1]
            if a1 not in ap.params["ancestors"]:
                print "\n. ERROR: you specified a comparison between ancestors", a1, "and", a2, "but", a1,"was not defined in the ANCESTORS line."
                write_error(con, "you specified a comparison between ancestors " + a1, " and " + a2 + " but " + a1 + " was not defined in the ANCESTORS line.") 
                exit() 
            if a1 not in ap.params["ancestors"]:
                print "\n. ERROR: you specified a comparison between ancestors", a1, "and", a2, "but", a2,"was not defined in the ANCESTORS line."
                write_error(con, "you specified a comparison between ancestors " + a1, " and " + a2 + " but " + a2 + " was not defined in the ANCESTORS line.")
                exit()
    
    ergseqpath = get_setting_values(con, "ergseqpath")[0]
    if False == os.path.exists(  ergseqpath  ):
        print "\n. I could not find your sequences at", ergseqpath
        write_error(con, "I could not find your sequences at " + ergseqpath  )
        exit()

    if ap.params["start_motif"] == None:
        ap.params["start_motif"] = ""
    if ap.params["end_motif"] == None:
        ap.params["end_motif"] = ""
        
    for msa in ap.params["msa_algorithms"]:
        if msa == "muscle" and "muscle_exe" not in ap.params:
            print "\n. Something is wrong. Your config file doesn't have an executable path for MUSCLE."
            write_error(con, "Something is wrong. Your config file doesn't have an executable path for MUSCLE.")
            exit()
        elif msa == "muscle" and "muscle_exe" in ap.params:
            import_alignment_method(con, msa, ap.params["muscle_exe"])
        
        if msa == "msaprobs" and "msaprobs_exe" not in ap.params:
            print "\n. Something is wrong. Your config file doesn't have an executable path for MSAPROBS."
            write_error(con, "Something is wrong. Your config file doesn't have an executable path for MSAPROBS.")
            exit()
        elif msa == "msaprobs" and "msaprobs_exe" in ap.params:
            import_alignment_method(con, msa, ap.params["msaprobs_exe"])
        
        if msa == "prank" and "prank_exe" not in ap.params:
            print "\n. Something is wrong. Your config file doesn't have an executable path for PRANK."
            write_error(con, "Something is wrong. Your config file doesn't have an executable path for PRANK.")
            exit()
        elif msa == "prank" and "prank_exe" in ap.params:
            import_alignment_method(con, msa, ap.params["prank_exe"])
        
        if msa == "mafft" and "mafft_exe" not in ap.params:
            print "\n. Something is wrong. Your config file doesn't have an executable path for PRANK."
            write_error(con, "Something is wrong. Your config file doesn't have an executable path for PRANK.")
            exit()
        elif msa == "mafft" and "mafft_exe" in ap.params:
            import_alignment_method(con, msa, ap.params["mafft_exe"])
        
    
    if "map2pdb" not in ap.params:
        ap.params["map2pdb"] = {}
        
    if "pdbtoolsdir" in ap.params and "pymol_exe" in ap.params:
        ap.params["do_pdb_analysis"] = True
    else:
        ap.params["do_pdb_analysis"] = False
        
    if "zorro_thresholds" not in ap.params:
        ap.params["zorro_thresholds"] = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.5, 0.6, 0.75, 1.0]
    
    for t in ap.params["zorro_thresholds"]:
        sql = "insert into Settings (keyword, value) VALUES('zorro_threshold'," + t.__str__() + ")"
        cur.execute(sql)
    con.commit()
        
    return ap

def print_config(ap):
    for p in ap.params:
        print p, ":", ap.params[p]

def setup_workspace(con):
    for msa in get_alignment_method_names(con):
        if False == os.path.exists(msa):
            os.system("mkdir " + msa)
        if False == os.path.exists(msa):
            print "I can't make the directory for output of " + msa
            write_error(con, "I can't make the directory for output of " + msa)
            exit()
    if False == os.path.exists("SCRIPTS"):
        os.system("mkdir SCRIPTS")
            