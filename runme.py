import random,sys,os
from argParser import *
from tools import *
from phyloasr import *

from read_config import *
from run_msas import *

ap = ArgParser(sys.argv)

""" Setup """
read_config_file( ap )
print_config(ap)
verify_config(ap)
setup_workspace(ap)

""" MSAs """
p = write_msa_commands(ap)
#run_script(p, ap)
trim_alignments(ap)

""" ML Trees """
p = write_raxml_commands(ap)
run_script(p, ap)

#
# continue here - move the raxml output into OUT.* directories.
#


"""
python run_msas.py ime2
source run_msa_cleanup.sh cmgc.erg.raw.fasta cmgc
python run_phyloasr.py
source compareanc_commands.sh
python write_ancestral_summary.py
python plot_pp_distro.2.py > ancestral_summary.txt
"""