"""
This is the Sqlite3 database engine.

Author: Victor Hanson-Smith
Email: victorhansonsmith@gmail.com
"""
import sqlite3 as lite
import os, sys
from splash import *

def build_db(dbpath = None):
    """Initializes all the tables. Returns the DB connection object.
    If tables already exist, they will NOT be overwritten."""
        
    if dbpath == None or dbpath == False:
        dbpath = "asr.db"
        print "\n. Creating a new database at", dbpath
    else:
        print "\n. Restoring the existing database at", dbpath

    con = lite.connect(dbpath)

    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS About(version TEXT)")
    cur.execute("INSERT OR REPLACE INTO About(version) VALUES('" + VERSION.__str__() + "')")
    
    cur.execute("create table if not exists Settings(keyword TEXT, value TEXT)")
    cur.execute("create table if not exists Log(id INTEGER primary key, time DATETIME DEFAULT CURRENT_TIMESTAMP,  message TEXT, code INT)")
    cur.execute("create table if not exists ErrorLog(id INTEGER primary key, time DATETIME DEFAULT CURRENT_TIMESTAMP,  message TEXT, code INT)")
    
    # Sequence Alignment
    """Note: each Taxa can have only one sequence; see the 'unique' modified on taxonid table fields."""
    cur.execute("create table if not exists Taxa(id INTEGER primary key autoincrement, fullname TEXT unique, shortname TEXT unique)")
    cur.execute("create table if not exists OriginalSequences(id INTEGER primary key autoincrement, taxonid INT unique, sequence TEXT, datatype INT)") # taxonid is the ID of a row in Taxa; datatypes: 0=nucleotide, 1=amino acid
    cur.execute("create table if not exists AlignmentMethods(id INTEGER primary key autoincrement, name TEXT unique, exe_path TEXT)")
    cur.execute("create table if not exists AlignedSequences(id INTEGER primary key autoincrement, taxonid INT, alsequence TEXT, almethod INT)") # these sequences contain indels


    """A SiteSet is a collection of sites in the alignment. For example, after trimming the sequences, instead of storing a new
        instance of the sequences, a SiteSet can be used to define the trimmed boundaries.
        SiteSets are declared in the table SiteSets, and then defined for each alignment in the table SiteSetsAlignment."""
    cur.execute("create table if not exists SiteSets(id INTEGER primary key autoincrement, setname TEXT unique)") # sets of sequence sites
    
    """Defines a set of sites, 'from' site to the 'to' site, inclusive. A set/sequence pair can have
        multiple from/to pairs in this table, for example if the SiteSet is not contiguous."""
    cur.execute("create table if not exists SiteSetsAlignment(setid INTEGER, almethod INT, fromsite INT, tosite INT)") # alseqid is the ID of an aligned sequence
    
    """These tables are for ZORRO and FastTree analysis."""
    cur.execute("create table if not exists AlignmentSiteScoringMethods(id INTEGER primary key autoincrement, name TEXT unique)")
    cur.execute("create table if not exists AlignmentSiteScores(almethodid INTEGER, scoringmethodid INTEGER, site INT, score FLOAT)")
    cur.execute("create table if not exists ZorroThreshStats(almethod INTEGER, thresh FLOAT, min_score FLOAT, nsites INT)")
    cur.execute("create table if not exists ZorroThreshFasttree(treeid INTEGER primary key autoincrement, almethod INTEGER, thresh FLOAT, newick TEXT)")
    cur.execute("create table if not exists FasttreeStats(treeid INTEGER, mean_bootstrap FLOAT, sum_of_branches FLOAT)")
    cur.execute("create table if not exists FasttreeSymmetricDistances(treeid1 INTEGER, treeid2 INTEGER, distance FLOAT)")
    cur.execute("create table if not exists FasttreeRFDistances(treeid1 INTEGER, treeid2 INTEGER, distance FLOAT)")

    
    # a list of phylo softwares, such as PhyML and RAxML
    cur.execute("create table if not exists PhyloSoftwares(id INTEGER primary key autoincrement, name TEXT unique, exe_path TEXT)") 
    # a list of phylogenetic models, such as LG+I+G
    cur.execute("create table if not exists PhyloModels(modelid INTEGER primary key autoincrement, name TEXT unique)") 
    # command-line strings used to invoke phylo models in different softwares
    cur.execute("create table if not exists PhyloModelsSoftware(modelid INTEGER, softwareid INTEGER, software_model_string TEXT)") 
    cur.execute("create table if not exists ParsimonyPhylogenies(id INTEGER primary key autoincrement, seqsetid INTEGER, newick TEXT)") # seqsetid is the site set used to make this tree
    cur.execute("create table if not exists UnsupportedMlPhylogenies(id INTEGER primary key autoincrement, almethod INT, phylomodelid INTEGER, newick TEXT)")

    cur.execute("create table if not exists TreeMl(mltreeid INTEGER, likelihood FLOAT)") # mltreeid is an ID for an UnsupportedMlPhylogeny
    cur.execute("create table if not exists TreeAlpha(mltreeid INTEGER, alpha FLOAT)")
    cur.execute("create table if not exists TreePP(mltreeid INTEGER, pp FLOAT)")
    
    # a list of branch support methods, such as aLRT and PP
    cur.execute("create table if not exists BranchSupportMethods(id INTEGER primary key autoincrement, name TEXT unique)")
    
    # maps unsupported trees to a support method
    cur.execute("create table if not exists SupportedMlPhylogenies(id INTEGER primary key autoincrement, unsupportedmltreeid INTEGER, newick TEXT, supportmethodid INT)")
    
    # A TaxaGroup is an "ingroup" of related taxa.
    # A special group named 'outgroup' is also created.
    cur.execute("create table if not exists TaxaGroups(id INTEGER primary key autoincrement, name TEXT unique)") # the name of the TaxaGroup will be the name of the group's ancestor.
    cur.execute("create table if not exists GroupsTaxa(groupid INTEGER, taxonid INTEGER)") # groupid is the ID of either an ingroup or outgroup
    cur.execute("create table if not exists GroupSeedTaxa(groupid, INTEGER, seed_taxonid INTEGER)") # each group can have one, or multiple, seed taxon
    
    cur.execute("create table if not exists Ancestors(id INTEGER primary key autoincrement, almethod INT, phylomodel INT, name TEXT)")
    cur.execute("create table if not exists AncestralStates(ancid INTEGER, site INT, state CHAR, pp FLOAT)") # site is specific to the almethod for this Ancestor
    cur.execute("create table if not exists AncestralCladogram(id INTEGER primary key autoincrement, unsupportedmltreeid INTEGER, newick TEXT)")
    
    """Some ancestors have alias names, in addition to their default number-based name."""
    cur.execute("create table if not exists AncestorsAlias(ancid INTEGER, alias TEXT)") # alias names for this ancestor
    """Some ancestors are special, with pre-defined mappings to known ingroups and outgroups."""
    cur.execute("create table if not exists AncestorsGroups(ancid INTEGER unique, ingroupid INTEGER, outgroupid INTEGER)") # some ancestors, but not all, will have a mapping to known taxa groups.
    
    cur.execute("create table if not exists CompareAncestors(ancid1 INTEGER, ancid2 INTEGER)")
    
    
    """For HTML visualization"""
    cur.execute("create table if not exists HtmlPieces(id INTEGER primary key autoincrement, keyword TEXT unique, htmlstring TEXT)")
    con.commit()
    
    return con