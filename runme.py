import random,sys,os
from tools import *
import time
from phyloasr import *

from splash import *
from read_config import *
from run_msas import *
from fscores import *
from asr_bayes import *
from aws_tools import *
from struct_analysis import *
from index_mutations import *
#from pythoscape_tools import *
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

enable_aws = ap.getOptionalArg("--enable_aws")
if enable_aws == None:
    enable_aws = False
else:
    enable_aws = True
    
S3_BUCKET = None # the name of a bucket in S3
S3_KEYBASE = None # the base of a key string for this job, presumably this string will be the job ID.
if enable_aws:
    S3_BUCKET = ap.getOptionalArg("--s3_bucket")
    if S3_BUCKET == None:
        print "\n.You enabled AWS with --enable_aws, but you did not specify an S3 bucket with --s3_bucket"
        exit()
    S3_KEYBASE = ap.getOptionalArg("--s3_keybase")
    if S3_KEYBASE == None:
        print "\n.You enabled AWS with --enable_aws, but you did not specify an S3 keybase with --s3_keybase"
        exit()

if enable_aws:
    aws_update_status("Launching the PhyloBot Pipeline", S3_BUCKET, S3_KEYBASE)

"""Restore/build the database."""
dbpath = ap.getOptionalArg("--dbpath") 
if dbpath == None or dbpath == False:
    dbpath = "asr.db"
    print "\n. Creating a new database at", dbpath
else:
    print "\n. Restoring the existing database at", dbpath
con = build_db( dbpath )

""" Setup """
read_config_file( con, ap )
print_config(ap)
verify_config(con, ap)
verify_all_exe(con)
setup_workspace(con)

if enable_aws:
    aws_update_status("Importing Sequences", S3_BUCKET, S3_KEYBASE)
    aws_checkpoint(0, S3_BUCKET, S3_KEYBASE)

if jump <= 0 and stop > 0:
    write_log(con, "Checkpoint: reading sequences.")
    ap.params["checkpoint"] = -1
    ap.params["pending_checkpoint"] = 0
    write_log(con, "Reading input sequences")
    import_and_clean_erg_seqs(con, ap)

verify_erg_seqs(con, ap)
write_log(con, "Checkpoint: configuration is OK.")
if enable_aws:
    aws_update_status("Sequences and Configuration are OK", S3_BUCKET, S3_KEYBASE)
    push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
    aws_checkpoint(1, S3_BUCKET, S3_KEYBASE)

"""Note: Continue Here: An imported SQL DB could be retrieved right here, avoiding all the prior steps,
    and avoiding the need to keep the original FASTA and the configuration file around.
    This would need to be added as a command-line option."""

if jump != None:
    print "\n. Jumping to Step " + jump.__str__()

# if jump <= 0.5 and stop > 0.5:
#     write_log(con, "Checkpoint: building sequence similarity network")
#     run_similarity_network_analysis(con)
# if jump <= 0.6 and stop > 0.6:
#     write_log(con, "Checkpoint: checking the output from the similarity network analysis.")
#     check_similarity_network_analysis(con)
    
""" MSAs """
if jump <= 1 and stop > 1:
    write_log(con, "Checkpoint: aligning sequences.")
    if enable_aws:
        aws_update_status("Aligning Sequences", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
    ap.params["checkpoint"] = 0
    ap.params["pending_checkpoint"] = 1
    p = write_msa_commands(con)
    run_script(p)

if jump <= 1.1 and stop > 1.1:
    write_log(con, "Checkpoint: checking aligned sequences.")
    if enable_aws:
        aws_update_status("Checking the Aligned Sequences", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(1.1, S3_BUCKET, S3_KEYBASE)
    ap.params["checkpoint"] = 1
    ap.params["pending_checkpoint"] = 1.1
    check_aligned_sequences(con)
    
if jump <= 1.11 and stop > 1.11:
    write_log(con, "Building a site map between different alignments.")
    if enable_aws:
        aws_update_status("Building a Site Map Between Different Alignments", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(1.11, S3_BUCKET, S3_KEYBASE)
    build_site_map(con)

if jump <= 1.9 and stop > 1.91:
    clear_sitesets(con)

if jump <= 2 and stop > 2.1:
    write_log(con, "Checkpoint: trimming alignments.")
    """Trim the alignments to match the seed sequence(s)."""
    trim_alignments(con)

if ap.getOptionalToggle("--skip_zorro"):
    if jump <= 2.2:
        write_log(con, "Checkpoint: skipping ZORRO analysis.")
        bypass_zorro(con)

if False == ap.getOptionalToggle("--skip_zorro"):
    """Use ZORRO to the find the phylogenetically informative sites."""
    if jump <= 2.41 and stop > 2.41:
        write_log(con, "Checkpoint: starting ZORRO analysis")
        if enable_aws:
            aws_update_status("Running the ZORRO Analysis", S3_BUCKET, S3_KEYBASE)
            push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
            aws_checkpoint(2.41, S3_BUCKET, S3_KEYBASE)
        p = build_zorro_commands(con)
        run_script(p)
    if jump <= 2.42 and stop > 2.42:
        import_zorro_scores(con)
    if jump <= 2.43 and stop > 2.43:
        plot_zorro_stats(con)
    if jump <= 2.44 and stop > 2.44:
        p = build_fasttree4zorro_commands(con)
        run_script(p)
    if jump <= 2.45 and stop > 2.45:
        analyze_zorro_fasttrees(con)
    if jump <= 2.46 and stop > 2.46:
        measure_fasttree_distances(con)
    if jump <= 2.47 and stop > 2.47:
        compare_fasttrees(con)
    if jump <= 2.48 and stop > 2.48:
        cleanup_zorro_analysis(con)
    if jump <= 2.49 and stop > 2.49:
        write_alignment_for_raxml(con)

if jump <= 2.7 and stop > 2.7:
    write_log(con, "Checkpoint: converting all FASTA to PHYLIP")
    if enable_aws:
        aws_update_status("Converting all FASTA to PHYLIP", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(2.7, S3_BUCKET, S3_KEYBASE)
    convert_all_fasta_to_phylip(con)

""" ML Trees """
if jump <= 3 and stop > 3:
    write_log(con, "Checkpoint: finding ML trees with RAxML (be patient)")
    if enable_aws:
        aws_update_status("Finding ML Phylogenies with RAxML (this can take a while)", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(3, S3_BUCKET, S3_KEYBASE)
    ap.params["checkpoint"] = 2
    ap.params["pending_checkpoint"] = 3
    p = write_raxml_commands(con)
    run_script(p)

if jump <= 3.1 and stop > 3.1:
    write_log(con, "Checkpoint: checking RAxML output")
    if enable_aws:
        aws_update_status("Checking RAxML Output", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(3.1, S3_BUCKET, S3_KEYBASE)
    """ML trees, part 2"""
    check_raxml_output(con)
    get_mlalpha_pp(con)

""" Branch Support """
if jump <= 4 and stop > 4:
    write_log(con, "Checkpoint: calculating branch support")
    if enable_aws:
        aws_update_status("Calculating Phylogenetic Branch Support", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(4, S3_BUCKET, S3_KEYBASE)
    ap.params["checkpoint"] = 3
    ap.params["pending_checkpoint"] = 4
    x = calc_alrt(con)
    run_script(x)
    calc_alr(con)

if jump <= 4.1 and stop > 4.1:
    import_supported_trees(con)

""" A.S.R. """
if jump <= 5 and stop > 5:
    write_log(con, "Checkpoint: reconstructing ancestral sequences")
    if enable_aws:
        aws_update_status("Reconstructing Ancestral Sequences", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(5, S3_BUCKET, S3_KEYBASE)
    ap.params["checkpoint"] = 4
    ap.params["pending_checkpoint"] = 5
    x = get_asr_commands(con)
    run_script(x)
    
if jump <= 5.1 and stop > 5.1:
    check_asr_output(con)

if jump <= 5.11 and stop > 5.11:
    if enable_aws:
        aws_update_status("Matching Ancestors Across Models", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(5.11, S3_BUCKET, S3_KEYBASE)
    match_ancestors_across_models(con)
    
if jump <= 5.2 and stop > 5.3:
    write_log(con, "Checkpoint: extracting relevant ancestors")
    if enable_aws:
        aws_update_status("Extracting Relevant Ancestors", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(5.2, S3_BUCKET, S3_KEYBASE)
    ap.params["checkpoint"] = 5
    ap.params["pending_checkpoint"] = 5.1
    x = get_getanc_commands(con)
    run_script(x)
    check_getanc_output(con)


""" Predict sites of functional evolution """
if jump <= 6 and stop > 6:
    if "compareanc" in ap.params:
        write_log(con, "Checkpoint: performing Df ranking for functional loci")
        if enable_aws:
            aws_update_status("Calculating Df Ranks for Functional Loci", S3_BUCKET, S3_KEYBASE)
            push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
            aws_checkpoint(6, S3_BUCKET, S3_KEYBASE)
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

#if jump <= 6.2 and stop > 6.2:
#    index_mutations(con)

"""
dN/dS Tests
"""
x = get_setting_values(con, "ergntpath")
if x != None:
    if jump <= 6.5 and stop > 6.5:
        """Do dn/ds test."""
        write_log(con, "Checkpoint: performing dN/dS tests")
        if enable_aws:
            aws_update_status("Calculating dN/dS Test Metrics", S3_BUCKET, S3_KEYBASE)
            push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
            aws_checkpoint(6.5, S3_BUCKET, S3_KEYBASE)
        x = get_dnds_commands(con)
        run_script(x)
    if jump <= 6.6 and stop > 6.6:
        write_log(con, "Checkpoint: checking dN/dS output")
        parse_dnds_results(con)
        
    if jump <= 6.8 and stop > 6.8:
        write_log(con, "Checkpoint: comparing Df to dN/dS")
        if enable_aws:
            aws_update_status("Comparing Df Ranks to dN/dS Test Scores", S3_BUCKET, S3_KEYBASE)
            push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
            aws_checkpoint(6.8, S3_BUCKET, S3_KEYBASE)
        setup_compare_functional_loci(con)
        compare_functional_loci(con)
else:
    write_log(con, "Codon sequences were not given by the user; I will skip dN/dS analysis.")



"""
    December 2014: The new Django-version of this code should stop here.
    Rather than building static HTML pages (code below), we'll use Django
    to produce dynamic HTML content on the fly.
"""


if jump <= 7 and stop > 7:
    write_log(con, "Checkpoint: cleaning-up residual files")
    if enable_aws:
        aws_update_status("Almost Done: Cleaning Residual Files", S3_BUCKET, S3_KEYBASE)
        push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
        aws_checkpoint(7, S3_BUCKET, S3_KEYBASE)
    cleanup(con) 

write_log(con, "Checkpoint: Analysis is complete.")
if enable_aws:
    push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE)
    aws_update_status("Finished", S3_BUCKET, S3_KEYBASE)
    aws_checkpoint(8, S3_BUCKET, S3_KEYBASE)
    
    """The job sends the stop signal to the SQS queue, essentially committing suicide."""
    sqs_stop(S3_KEYBASE, attempts=0)

exit()

"""
    December 2014: The following HTML-generation methods are disabled.
    The new architecture will dynamically generate HTML views, using
    Django, instead of pre-computing all the static HTML pages.
    This new dynamic approach will consume less disk space, and be
    more flexible, than the static approach.
"""
 
if jump <= 8.2 and stop > 8.3:
    ap.params["checkpoint"] = 7.1
    ap.params["pending_checkpoint"] = 7.2
    write_ancseq_fasta(con, ap)
 
if stop >= 9:
    ap.params["checkpoint"] = 100
    write_log(con, "Done")