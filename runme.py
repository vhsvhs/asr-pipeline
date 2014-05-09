import random,sys,os
from tools import *
from phyloasr import *

from splash import *
from read_config import *
from run_msas import *
from asr_bayes import *
from html_helper import *
from struct_analysis import *
from log import *

print_splash()

jump = ap.getOptionalArg("--jump") # Will start the script at this point.
if jump == False:
    jump = 0
else:
    jump = float( jump )

stop = ap.getOptionalArg("--stop") # Will stop the script upon reaching, but not completing, this step.
if stop == False:
    stop = 100
else:
    stop = float( stop )

""" Setup """
read_config_file( ap )
print_config(ap)
verify_config(ap)
setup_workspace(ap)

if jump <= 0 and stop > 0:
    init_log(ap)
    print "\n. Reading your FASTA sequences..."
    write_log(ap, 0, "Reading sequences")
    clean_erg_seqs(ap)
else:
    init_log(ap, overwrite=False)

    

""" MSAs """
if jump <= 1 and stop > 1:
    print "\n. Aligning sequences..."
    write_log(ap, 1, "Aligning sequences")
    p = write_msa_commands(ap)
    run_script(p)

if jump <= 1.1 and stop > 1.1:
    print "\n. Converting the alignments to PHYLIP..."
    write_log(ap, 1.1, "Creating Phylip-formatted versions of the alignments.")
    convert_all_fasta_to_phylip(ap)

#if jump <= 2:
#    trim_alignments(ap)

""" ML Trees """
if jump <= 3 and stop > 3:
    print "\n. Inferring ML phylogenies with RAxML..."
    write_log(ap, 3, "Inferring ML phylogenies with RAxML.")
    p = write_raxml_commands(ap)
    run_script(p)

""" Branch Support """
if jump <= 4 and stop > 4:
    print "\n. Calculating aLRT branch support with PhyML..."
    write_log(ap, 4, "Calculating aLRT branch support values with PhyML.")
    get_mlalpha_pp(ap)
    x = calc_alrt(ap)
    run_script(x)
    calc_alr(ap)

""" A.S.R. """
if jump <= 5 and stop > 5:
    print "\n. Reconstructing ancestral sequences..."
    write_log(ap, 5, "Reconstructing ancestral sequences, using Lazarus.")
    x = get_asr_commands(ap)
    run_script(x)

if jump <= 5.1 and stop > 5.1:
    write_log(ap, 5.1, "Extracting relevant ancestors")
    x = get_getanc_commands(ap)
    run_script(x)

if jump <= 5.2 and stop > 5.2:
    run_asr_bayes(ap)

""" Predict sites of functional evolution """
if jump <= 6 and stop > 6:
    if "compareanc" in ap.params:
        if (jump > 4):
            get_mlalpha_pp(ap)
        write_log(ap, 6, "Setting up PDB maps")
        setup_pdb_maps(ap)
        write_log(ap, 6.1, "Screening for functional loci.")
        x = get_compareanc_commands(ap)
        os.system( ap.params["run_exe"] + " " + x)
        #run_script(x)

""" Build an HTML Report """
if jump <= 7 and stop > 7:
    write_log(ap, 7, "Writing an HTML report.")
    write_css()
    write_index()
    write_alignments()
    write_treesancs()
    write_ancestors_indi() # write individual ancestor pages

if jump <= 7.1 and stop > 7.1:
    for pair in ap.params["compareanc"]:
        write_anccomp_indi(pair, ap)
        write_mutations_indi(pair, ap)

if jump <= 7.2 and stop > 7.3:
    write_ancseq_fasta(ap)

write_log(ap, 8, "Complete!")