import random,sys,os
from argParser import *
from tools import *
from phyloasr import *

from splash import *
from read_config import *
from run_msas import *

ap = ArgParser(sys.argv)
print_splash()

jump = ap.getOptionalArg("--jump")
if jump == False:
    jump = 0
else:
    jump = float( jump )

""" Setup """
read_config_file( ap )
print_config(ap)
verify_config(ap)
setup_workspace(ap)

""" MSAs """
if jump <= 1:
    print "\n. Aligning sequences..."
    p = write_msa_commands(ap)
    run_script(p, ap)

if jump <= 2:
    trim_alignments(ap)

""" ML Trees """
if jump <= 3:
    print "\n. Inferring ML phylogenies with RAxML..."
    p = write_raxml_commands(ap)
    run_script(p, ap)

if jump <= 4:
    print "\n. Calculating aLRT branch support with PhyML..."
    x = post_raxml(ap)
    run_script(x, ap)
    calc_alr(ap)

if jump <= 5:
    if "ancestors" in ap.params:
        print "\n. Reconstructing ancestral sequences..."
        x = get_asr_commands(ap)
        run_script(x, ap)
        
        x = get_getanc_commands(ap)
        run_script(x, ap)


"""
python run_msas.py ime2
source run_msa_cleanup.sh cmgc.erg.raw.fasta cmgc
python run_phyloasr.py
source compareanc_commands.sh
python write_ancestral_summary.py
python plot_pp_distro.2.py > ancestral_summary.txt
"""