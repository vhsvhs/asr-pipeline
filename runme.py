import random,sys,os
from tools import *
from phyloasr import *

from splash import *
from read_config import *
from run_msas import *

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
    print "\n. Reading your FASTA sequences..."
    clean_erg_seqs(ap)
    

""" MSAs """
if jump <= 1 and stop > 1:
    print "\n. Aligning sequences..."
    p = write_msa_commands(ap)
    run_script(p)

if jump <= 1.1 and stop > 1.1:
    print "\n. Converting the alignments to PHYLIP..."
    convert_all_fasta_to_phylip(ap)

#if jump <= 2:
#    trim_alignments(ap)

""" ML Trees """
if jump <= 3 and stop > 3:
    print "\n. Inferring ML phylogenies with RAxML..."
    p = write_raxml_commands(ap)
    run_script(p)

""" Branch Support """
if jump <= 4 and stop > 4:
    print "\n. Calculating aLRT branch support with PhyML..."
    get_mlalpha_pp(ap)
    x = calc_alrt(ap)
    run_script(x)
    calc_alr(ap)

""" A.S.R. """
if jump <= 5 and stop > 5:
    if "ancestors" in ap.params:
        print "\n. Reconstructing ancestral sequences..."
        x = get_asr_commands(ap)
        run_script(x)
        
        x = get_getanc_commands(ap)
        run_script(x)

""" Predict sites of functional evolution """
if jump <= 6 and stop > 6:
    if "compareanc" in ap.params:
        if (jump > 4):
            get_mlalpha_pp(ap)
        x = get_compareanc_commands(ap)
        run_script(x)

""" Build an HTML Report """
if jump <= 7 and stop > 7:
    from html_helper import *
    write_css()
    write_index()
    write_alignments()
    write_treesancs()
    write_ancestors_indi() # write individual ancestor pages


"""The maximum jump step index is 100.  See the command-line argument --stop"""

"""
python run_msas.py ime2
source run_msa_cleanup.sh cmgc.erg.raw.fasta cmgc
python run_phyloasr.py
source compareanc_commands.sh
python write_ancestral_summary.py
python plot_pp_distro.2.py > ancestral_summary.txt
"""