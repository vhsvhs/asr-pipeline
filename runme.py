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
from plots import *
print_splash()

jump = ap.getOptionalArg("--jump") # Will start the script at this point.
if jump == None:
    jump = 0
else:
    jump = float( jump )

stop = ap.getOptionalArg("--stop") # Will stop the script upon reaching, but not completing, this step.
if stop == None:
    stop = 100
else:
    stop = float( stop )

"""Restore/build the database."""
con = build_db( dbpath = ap.getOptionalArg("--dbpath") )

""" Setup """
read_config_file( con, ap )
print_config(ap)
verify_config(con, ap)
setup_workspace(ap)

if jump <= 0 and stop > 0:
    ap.params["checkpoint"] = -1
    ap.params["pending_checkpoint"] = 0
    print "\n. Reading your FASTA sequences..."
    write_log(con, "Reading sequences")
    import_and_clean_erg_seqs(con, ap)

verify_erg_seqs(con)

""" MSAs """
if jump <= 1 and stop > 1:
    ap.params["checkpoint"] = 0
    ap.params["pending_checkpoint"] = 1
    print "\n. Aligning sequences..."
    write_log(con, "Aligning sequences")
    p = write_msa_commands(con, ap)
    run_script(p)

if jump <= 1.1 and stop > 1.1:
    ap.params["checkpoint"] = 1
    ap.params["pending_checkpoint"] = 1.1
    check_aligned_sequences(con, ap)

if jump <= 1.9 and stop > 1.91:
    clear_sitesets(con, ap)

if jump <= 2 and stop > 2.1:
    """Trim the alignments to match the seed sequence(s)."""
    trim_alignments(con, ap)

"""Use ZORRO to the find the phylogenetically informative sites."""
if jump <= 2.1 and stop > 2.2:
    p = build_zorro_commands(con, ap)
    run_script(p)
if jump <= 2.2 and stop > 2.3:
    import_zorro_scores(con)
if jump <= 2.3 and stop > 2.4:
    plot_zorro_stats(con)
if jump <= 2.4 and stop > 2.5:
    p = build_raxml4zorro_commands(con)
    run_script(p)
if jump <= 2.5 and stop > 2.6:
    analyze_zorro_raxml(con)
if jump <= 2.6 and stop > 2.7:
    cleanup_zorro_analysis(con)

exit()

if jump <= 2.5 and stop > 2.6:
    convert_all_fasta_to_phylip(ap)

""" ML Trees """
if jump <= 3 and stop > 3:
    ap.params["checkpoint"] = 2
    ap.params["pending_checkpoint"] = 3
    print "\n. Inferring ML phylogenies with RAxML..."
    write_log(con, "Inferring ML phylogenies with RAxML.")
    p = write_raxml_commands(ap)
    run_script(p)
    check_raxml_output(ap)

""" Branch Support """
if jump <= 4 and stop > 4:
    ap.params["checkpoint"] = 3
    ap.params["pending_checkpoint"] = 4
    print "\n. Calculating aLRT branch support with PhyML..."
    write_log(con, "Calculating aLRT branch support values with PhyML.")
    get_mlalpha_pp(ap)
    x = calc_alrt(ap)
    run_script(x)
    calc_alr(ap)

""" A.S.R. """
if jump <= 5 and stop > 5:
    ap.params["checkpoint"] = 4
    ap.params["pending_checkpoint"] = 5
    print "\n. Reconstructing ancestral sequences..."
    write_log(con, "Reconstructing ancestral sequences, using Lazarus.")
    x = get_asr_commands(ap)
    run_script(x)
    #
    # to-do: insert check that ancestors were created
    #

if jump <= 5.1 and stop > 5.1:
    ap.params["checkpoint"] = 5
    ap.params["pending_checkpoint"] = 5.1
    write_log(con, "Extracting relevant ancestors")
    x = get_getanc_commands(ap)
    run_script(x)
    (flag, msg) = check_getanc_output(ap)
    if not flag:
        write_error(ap, msg)
        exit()
    
    #
    # to-do: insert check that ancestors were gotten.
    #

if jump <= 5.2 and stop > 5.2:
    ap.params["checkpoint"] = 5.1
    ap.params["pending_checkpoint"] = 5.2
    run_asr_bayes(ap)

""" Predict sites of functional evolution """
if jump <= 6 and stop > 6:
    if "compareanc" in ap.params:
        if "runid_alpha" not in ap.params or "runid_pp" not in ap.params:
            read_lnl_summary(ap)
        ap.params["checkpoint"] = 5.2
        ap.params["pending_checkpoint"] = 6
        write_log(con, "Setting up PDB maps")
        setup_pdb_maps(ap)
        write_log(con, "Screening for functional loci.")
        x = get_compareanc_commands(ap)
        args = x.split()
        ap.params["run_exe"]
        #proc = subprocess.Popen( args, preexec_fn=os.setsid ) # see http://pymotw.com/2/subprocess/
        #os.system( ap.params["run_exe"] + " " + x)
        run_script(x)

""" Build an HTML Report """
if jump <= 7 and stop > 7:
    ap.params["checkpoint"] = 6
    ap.params["pending_checkpoint"] = 7
    write_log(con, "Writing an HTML report.")
    write_css()
    write_index()
    write_alignments()
    write_treesancs()
    write_ancestors_indi() # write individual ancestor pages

if jump <= 7.1 and stop > 7.1:
    ap.params["checkpoint"] = 7
    ap.params["pending_checkpoint"] = 7.1
    if "compareanc" in ap.params:
        for pair in ap.params["compareanc"]:
            write_anccomp_indi(pair, ap)
            write_mutations_indi(pair, ap)

if jump <= 7.2 and stop > 7.3:
    ap.params["checkpoint"] = 7.1
    ap.params["pending_checkpoint"] = 7.2
    write_ancseq_fasta(ap)

if stop >= 8:
    ap.params["checkpoint"] = 100
    write_log(con, "Done")