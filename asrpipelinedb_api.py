"""
asrpipelinedb_api.py
This is the API for the asrpipeline Sqlite3 database.

Author: Victor Hanson-Smith
Email: victorhansonsmith@gmail.com
"""

import sqlite3 as lite
import os, re, sys
from log import *

from asrpipelinedb import *

def get_setting_values(con, keyword):
    """Returns a list of one, or more, values for the setting keyword"""
    cur = con.cursor()
    sql = "select value from Settings where keyword='" + keyword + "'"
    cur.execute(sql)
    x = cur.fetchall()
    values = []
    for ii in x:
        values.append( ii[0] )
    if values.__len__() == 0:
        return None
    return values

def add_setting_value(con, keyword, value, unique=False):
    """Adds a keyword/value pair to the Settings table.
    If unique==True, then it will overwrite any existing entries for keyword."""
    cur = con.cursor()
    if unique == True:
        sql = "delete from Settings where keyword='" + keyword + "'"
        cur.execute(sql)
        con.commit()
    sql = "insert into Settings (keyword, value) VALUES('" + keyword + "','" + value.__str__() + "')"
    cur.execute(sql)
    con.commit()

def build_aligned_codon_seq(con, taxonid):
    pass
    

def import_original_seq(con, shortname, sequence, datatype=1):
    """Returns the taxonID of the newly imported sequence."""
    cur = con.cursor()
    
    """Have we seen this sequence name before?"""
    sql = "select count(*) from Taxa where shortname='" + shortname + "'"
    cur.execute(sql)
    count = cur.fetchone()[0]
    if count == 0:
        sql = "insert into Taxa (fullname, shortname) VALUES('" + shortname + "','" + shortname + "')"
        cur.execute(sql)
        con.commit()
    
    taxonid = get_taxonid(con, shortname)
    
    sql = "select count(*) from OriginalSequences where taxonid=" + taxonid.__str__() + " and datatype=" + datatype.__str__()
    cur.execute(sql)
    count = cur.fetchone()[0]
    if count == 0:
        sql = "insert into OriginalSequences (taxonid, sequence, datatype) VALUES(" + taxonid.__str__() + ",'" + sequence + "'," + datatype.__str__() + ")"
        cur.execute(sql)
        con.commit()
    
    return taxonid

def import_aligned_seq(con, taxonid, almethodid, seq, datatype=1):
    cur = con.cursor()    
    sql = "select count(*) from AlignedSequences where taxonid=" + taxonid.__str__() + " and almethod=" + almethodid.__str__() + " and datatype=" + datatype.__str__()
    cur.execute(sql)
    count = cur.fetchone()[0]
    if count == 0:      
        sql = "insert or replace into AlignedSequences (taxonid, alsequence, almethod, datatype) VALUES("
        sql += taxonid.__str__() + ",'" + seq + "'," + almethodid.__str__() + "," + datatype.__str__() + ")"
        cur.execute(sql)
        con.commit()
    
def import_alignment_method(con, name, exe):
    cur = con.cursor()
    sql = "select count(*) from AlignmentMethods where name='" + name + "'"
    cur.execute(sql)
    count = cur.fetchone()[0]
    if count == 0:  
        sql = " insert into AlignmentMethods(name, exe_path) VALUES("
        sql += "'" + name + "','" + exe + "')"
        cur.execute(sql)
        con.commit()

def get_taxonid(con, shortname):
    cur = con.cursor()
    sql = "select id from Taxa where shortname='" + shortname + "'"
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    return x[0][0]

def get_taxon_name(con, taxonid):
    cur = con.cursor()
    sql = "select shortname from Taxa where id=" + taxonid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    return x[0][0]

def get_ancestorid(con, alias, almethod, phylomodel):
    cur = con.cursor()
    sql = "select ancid from AncestorsAlias where alias='" + alias + "' and ancid in (select id from Ancestors where almethod=" + almethod.__str__() + " and phylomodel=" + phylomodel.__str__() + ")"
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    return x[0][0]

def get_ancestorname(con, ancid):
    cur = con.cursor()
    sql = "select name from Ancestors where ancid=" + ancid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    return x[0][0]

def get_ancestral_alias(con, ancid):
    cur = con.cursor()
    sql = "select alias from AncestorsAlias where ancid=" + ancid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    return x[0][0]

def get_ancestral_comparison_pairs(con):
    cur = con.cursor()
    sql = "select alias1, alias2 from CompareAncestors"
    cur.execute(sql)
    x = cur.fetchall()
    pairs = []
    for ii in x:
        pairs.append(  (ii[0],ii[1])  )
    return pairs

def get_ingroup_ids(con):
    cur = con.cursor()
    sql = "select id from TaxaGroups where name <> 'outgroup'"
    cur.execute(sql)
    ids = []
    for ii in cur.fetchall():
        ids.append(ii[0])
    return ids

def get_original_seq(con, taxonid, datatype=1):
    cur = con.cursor()
    sql = "select sequence from OriginalSequences where taxonid=" + taxonid.__str__() + " and datatype=" + datatype.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    return x[0][0] 

def get_aligned_seq(con, taxonid, almethodid,datatype=1):
    cur = con.cursor()
    sql = "select alsequence from AlignedSequences where taxonid=" + taxonid.__str__() + " and almethod=" + almethodid.__str__() + " and datatype=" + datatype.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    return x[0][0]    

def get_sequences(con, almethod = None, sitesets = [], sites = [], taxagroups = [], datatype=1):
    """
    Returns a hash of taxa:sequences.
    con = the sqlite3 db
    almethod: if None, then returns all OriginalSequences, ignoring the value of siteset.
            : if a number N, then it returns only aligned sequences generated by the alignment method N.
    sitesets = the IDs of a SiteSet. Only sites corresponding to that set will be
        included in the output sequences. If [], then it will return all sites
    sites = a list of sites, of which only these sites should be included. You can't specify
        both sites and sitesets
    taxagroups: if [], then all possible sequences are returned,
                else, only sequences in the groups listed will be returned.
    relative_to_trimmed: If True, then the site numbers provided in sites should be interpreted
        as being relative to the aligned sequences trimmed to the seed sequences.
    """
    cur = con.cursor()
    
    sql = "select "
    
    """Different queries for original sequences, versus aligned sequences."""
    if almethod == None:
        sql += " taxonid, sequence from OriginalSequences"
    else:
        if False == is_valid_almethod( con, almethod ):
            print "\n. Error (31): the alignment method", almethod, "is invalid."
            exit()
        sql += " taxonid, alsequence from AlignedSequences where almethod=" + almethod.__str__()

    if taxagroups != []:
        sql += " where taxonid in (select taxonid from GroupsTaxa "
        taxaids = []
        pieces = []
        for tgid in taxagroups:
            pieces.append( " where groupid=" + tgid )
        sql += " or ".join( pieces )
        sql += ")"
    sql += " and datatype=" + datatype.__str__()
    
    cur.execute(sql)
    x = cur.fetchall()
    
    taxa_seqs = {} # Key = taxonid, value = aligned or original sequence, untrimmed.
    for ii in x:
        taxonid = ii[0]
        seq = ii[1]
        taxa_seqs[ taxonid ] = seq
    
    if sitesets != []:
        siteranges = []
        for id in sitesets:
            siteranges += get_siteranges_in_set(con, id, almethod)
        
        for taxon in taxa_seqs:
            seq = taxa_seqs[taxon]
            trimseq = ""
            for r in siteranges:
                """Notice the -1 here; its to translate sites (1-based) to list indices (0-based)"""
                trimseq += seq[ r[0]-1 : r[1] ]
            taxa_seqs[taxon] = trimseq
    
    elif sites != [] and almethod != None:
        
        for taxon in taxa_seqs:
            seq = taxa_seqs[taxon]
            trimseq = ""
            for s in sites:
                """Note, the -1 in the following array references allow us to translate a site
                number (which are 1-based) to a list index (which are 0-based)"""
                print "asrpipelinedb.py 137:", almethod, s-1, seq.__len__()
                trimseq += seq[s-1]
            taxa_seqs[taxon] = trimseq        
    
    return taxa_seqs
        
    
def is_valid_almethod(con, methodid):
    cur = con.cursor()
    sql = "select count(*) from AlignmentMethods where id=" + methodid.__str__()
    cur.execute(sql)
    x = cur.fetchone()[0]
    if x > 0:
        return True
    else:
        return False

def get_taxaid_in_group(con, groupid):
    cur = con.cursor()
    sql = "select taxonid from GroupsTaxa where groupid=" + groupid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    taxaids = []
    for ii in x:
        taxaids.append( ii[0] )
    return taxaids

def get_taxaid_from_ancestor(con, ancid):
    cur = con.cursor()
    sql = "select taxonid from GroupsTaxa where groupid in (select ingroupid from AncestorsGroups where ancid=" + ancid.__str__() + " )"
    cur.execute(sql)
    x = cur.fetchall()
    taxonids = []
    for ii in x:
        taxonids.append( ii[0] )
    return taxonids
    

def get_siteranges_in_set(con, sitesetid, almethodid):
    """Returns a list of tuples [ (from,to), (from,to), etc. ]"""
    cur = con.cursor()
    sql = "select fromsite, tosite from SiteSetsAlignment where setid=" + sitesetid.__str__() + " and almethod=" + almethodid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    siteranges = []
    for ii in x:
        siteranges.append( (ii[0], ii[1]) )
    return siteranges

def get_lower_bound_in_siteset(con, sitesetid, almethodid):
    cur = con.cursor()
    sql = "select fromsite, tosite from SiteSetsAlignment where setid=" + sitesetid.__str__() + " and almethod=" + almethodid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    minsite = None
    for ii in x:
        if minsite == None:
            minsite = ii[0]
        elif minsite > ii[0]:
            minsite = ii[0]
    return minsite

def write_fasta(seqs, fpath):
    fout = open(fpath, "w")
    for s in seqs:
        fout.write(">" + s + "\n")
        fout.write(seqs[s] + "\n")
    fout.close()
    
def write_phylip(seqs, ppath):
    """Sanity check: are all the sequences the same length?"""
    l = None
    for s in seqs:
        if l == None:
            l = seqs[s].__len__()
        elif seqs[s].__len__() != l:
            print "\n. Error. (115) You are trying to write a phylip-formatted alignment using a set of sequences that are different lengths."
            exit()
    """Okay, write the phylip file."""
    fout = open(ppath, "w")
    fout.write( seqs.__len__().__str__() + "   " + l.__str__() + "\n")
    for s in seqs:
        fout.write(s + "   " + seqs[s] + "\n")
    fout.close()

def get_sitesetid(con, sitesetname):
    cur = con.cursor()
    sql = "select id from SiteSets where setname='" + sitesetname + "'"
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    else:
        return x[0][0]

def get_alignmentsitescores(con, almethod, scoringmethoid):
    cur = con.cursor()
    sql = "select site, score from AlignmentSiteScores where almethodid=" + almethod.__str__() + " and scoringmethodid=" + scoringmethoid.__str__()
    cur.execute(sql)
    y = cur.fetchall()
    site_scores = {}
    for jj in y:
        site_scores[ jj[0] ] = jj[1]
    #print "235: almethod", almethod, ", N sitescores=", site_scores.__len__().__str__() + ", from " + min( site_scores.keys()).__str__() + " to " + max( site_scores.keys()).__str__()
    return site_scores


def get_phylo_modelids(con):
    cur = con.cursor()
    sql = "select modelid from PhyloModels"
    cur.execute(sql)
    modelids = []
    x = cur.fetchall()
    for ii in x:
        modelids.append( ii[0] )
    return modelids

def get_phylo_modelnames(con):
    cur = con.cursor()
    sql = "select name from PhyloModels"
    cur.execute(sql)
    modelnames = []
    x = cur.fetchall()
    for ii in x:
        modelnames.append( ii[0] )
    return modelnames

def get_alignment_method_names(con):
    cur = con.cursor()
    sql = "select name from AlignmentMethods"
    cur.execute(sql)
    names = []
    x = cur.fetchall()
    for ii in x:
        names.append( ii[0] )
    return names

def get_alignment_method_ids(con):
    cur = con.cursor()
    sql = "select id from AlignmentMethods"
    cur.execute(sql)
    ids = []
    x = cur.fetchall()
    for ii in x:
        ids.append( ii[0] )
    return ids