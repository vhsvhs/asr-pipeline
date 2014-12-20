import random,sys,os
from tools import *
from phyloasr import *

from splash import *
from read_config import *
from run_msas import *
from fscores import *
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
setup_workspace(con)

if jump <= 0 and stop > 0:
    write_log(con, "Checkpoint: reading sequences.")
    ap.params["checkpoint"] = -1
    ap.params["pending_checkpoint"] = 0
    print "\n. Reading your FASTA sequences..."
    write_log(con, "Reading sequences")
    import_and_clean_erg_seqs(con, ap)


verify_erg_seqs(con, ap)

"""Note: An imported SQL DB could be retrieved right here, avoiding all the prior steps."""

""" MSAs """
if jump <= 1 and stop > 1:
    write_log(con, "Checkpoint: aligning sequences.")
    ap.params["checkpoint"] = 0
    ap.params["pending_checkpoint"] = 1
    print "\n. Aligning sequences..."
    write_log(con, "Aligning sequences")
    p = write_msa_commands(con)
    run_script(p)

if jump <= 1.1 and stop > 1.1:
    write_log(con, "Checkpoint: checking aligned sequences.")
    ap.params["checkpoint"] = 1
    ap.params["pending_checkpoint"] = 1.1
    check_aligned_sequences(con)

if jump <= 1.9 and stop > 1.91:
    clear_sitesets(con)

if jump <= 2 and stop > 2.1:
    write_log(con, "Checkpoint: trimming alignments.")
    """Trim the alignments to match the seed sequence(s)."""
    trim_alignments(con)

if ap.getOptionalToggle("--skip_zorro"):
    write_log(con, "Checkpoint: skipping ZORRO analysis.")
    bypass_zorro(con)

if False == ap.getOptionalToggle("--skip_zorro"):
    write_log(con, "Checkpoint: starting ZORRO analysis")
    """Use ZORRO to the find the phylogenetically informative sites."""
    if jump <= 2.1 and stop > 2.2:
        p = build_zorro_commands(con)
        run_script(p)
    if jump <= 2.2 and stop > 2.3:
        import_zorro_scores(con)
    if jump <= 2.3 and stop > 2.4:
        plot_zorro_stats(con)
    if jump <= 2.4 and stop > 2.5:
        p = build_fasttree4zorro_commands(con)
        run_script(p)
    if jump <= 2.5 and stop > 2.6:
        analyze_zorro_fasttrees(con)
    if jump <= 2.6 and stop > 2.7:
        measure_fasttree_distances(con)
    if jump <= 2.7 and stop > 2.8:
        compare_fasttrees(con)
    if jump <= 2.9 and stop > 2.91:
        cleanup_zorro_analysis(con)
        
    if jump <= 2.91 and stop > 2.92:
        write_alignment_for_raxml(con)

if jump <= 2.99 and stop > 2.99:
    write_log(con, "Checkpoint: converting all FASTA to PHYLIP")
    convert_all_fasta_to_phylip(con)

""" ML Trees """
if jump <= 3 and stop > 3:
    write_log(con, "Checkpoint: finding ML trees with RAxML")
    ap.params["checkpoint"] = 2
    ap.params["pending_checkpoint"] = 3
    p = write_raxml_commands(con)
    run_script(p)

if jump <= 3.1 and stop > 3.1:
    write_log(con, "Checkpoint: checking RAxML output")
    """ML trees, part 2"""
    check_raxml_output(con)
    get_mlalpha_pp(con)

""" Branch Support """
if jump <= 4 and stop > 4:
    write_log(con, "Checkpoint: calculating branch support")
    ap.params["checkpoint"] = 3
    ap.params["pending_checkpoint"] = 4

    x = calc_alrt(con)
    run_script(x)
    calc_alr(con)
    import_supported_trees(con)

""" A.S.R. """
if jump <= 5 and stop > 5:
    write_log(con, "Checkpoint: reconstructing ancestral sequences")
    ap.params["checkpoint"] = 4
    ap.params["pending_checkpoint"] = 5
    x = get_asr_commands(con)
    run_script(x)
if jump <= 5.1 and stop > 5.1:
    check_asr_output(con)
    
if jump <= 5.2 and stop > 5.3:
    write_log(con, "Checkpoint: extracting relevant ancestors")
    ap.params["checkpoint"] = 5
    ap.params["pending_checkpoint"] = 5.1
    x = get_getanc_commands(con)
    run_script(x)
    check_getanc_output(con)

if jump <= 5.9 and stop > 5.9:
    write_log(con, "Checkpoint: calculating Bayesian-sampled alternate ancestors")
    ap.params["checkpoint"] = 5.1
    ap.params["pending_checkpoint"] = 5.2
    run_asr_bayes(con, ap)

""" Predict sites of functional evolution """
if jump <= 6 and stop > 6:
    if "compareanc" in ap.params:
        write_log(con, "Checkpoint: performing Df ranking for functional loci")
        ap.params["checkpoint"] = 5.2
        ap.params["pending_checkpoint"] = 6
        
        #write_log(con, "Setting up PDB maps")
        # continue here -- retool this method to use SQL:
        #setup_pdb_maps(ap)
        
        write_log(con, "Screening for functional loci.")
        x = get_compareanc_commands(con)
        args = x.split()
        run_script(x)
if jump <= 6.1 and stop > 6.1:
    write_log(con, "Checkpoint: checking Df output")
    parse_compareanc_results(con)


sql = "select count(*) from Compare_DNDS_Fscores"
cur = con.cursor()
cur.execute(sql)
if cur.fetchone()[0] > 0:
    if jump <= 6.5 and stop > 6.5:
        """Do dn/ds test."""
        write_log(con, "Checkpoint: performing dN/dS tests")
        x = get_dnds_commands(con)
        run_script(x)
    if jump <= 6.6 and stop > 6.6:
        write_log(con, "Checkpoint: checking dN/dS output")
        parse_dnds_results(con)
        
    if jump <= 6.8 and stop > 6.8:
        write_log(con, "Checkpoint: comparing Df to dN/dS")
        setup_compare_functional_loci(con)
        compare_functional_loci(con)
else:
    write_log(con, "Codon sequences were not given by the user; I will skip dN/dS analysis.")

"""December 2014: The new Django-version of this code should stop here.
    Rather than building static HTML pages (code below), we'll use Django
    to produce dynamic HTML content on the fly."""
    
"""Continue here:

    write a method to delete all the residual files from the above analysis.
"""

if jump <= 7 and stop > 7:
    write_log(con, "Checkpoint: cleaning-up residual files")
    cleanup(con) 


exit()

""" Build an HTML Report """
if jump <= 7 and stop > 7:
    write_log(con, "Checkpoint - writing an HTML report.")
    ap.params["checkpoint"] = 6
    ap.params["pending_checkpoint"] = 7
    write_log(con, "Writing an HTML report.")
    write_css()
    write_index(con)
    write_alignments(con)
    write_treesancs(con)
    write_ancestors_indi(con) # write individual ancestor pages

if jump <= 7.1 and stop > 7.1:
    ap.params["checkpoint"] = 7
    ap.params["pending_checkpoint"] = 7.1
    if "compareanc" in ap.params:
        for pair in ap.params["compareanc"]:
            write_anccomp_indi(pair, con, ap)
            write_mutations_indi(pair, con, ap)

# if jump <= 7.2 and stop > 7.3:
#     ap.params["checkpoint"] = 7.1
#     ap.params["pending_checkpoint"] = 7.2
#     write_ancseq_fasta(con, ap)

if stop >= 8:
    ap.params["checkpoint"] = 100
    write_log(con, "Done")