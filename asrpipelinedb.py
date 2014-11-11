"""
This is the Sqlite3 database engine.

Author: Victor Hanson-Smith
Email: victorhansonsmith@gmail.com
"""
import sqlite3 as lite
import os, sys
from version import *

def build_db(dbpath = None):
    """Initializes all the tables. Returns the DB connection object.
    If tables already exist, they will NOT be overwritten."""
        
    if dbpath == None or dbpath == False:
        dbpath = "test.db"
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
    cur.execute("create table if not exists Taxa(id INTEGER primary key autoincrement, fullname TEXT, shortname TEXT unique)")
    cur.execute("create table if not exists OriginalSequences(id INTEGER primary key autoincrement, taxonid INT, sequence TEXT, datatype INT)") # taxonid is the ID of a row in Taxa; datatypes: 0=nucleotide, 1=amino acid
    cur.execute("create table if not exists AlignmentMethods(id INTEGER primary key autoincrement, name TEXT, exe_path TEXT)")
    cur.execute("create table if not exists AlignedSequences(id INTEGER primary key autoincrement, taxonid INT, alsequence TEXT, almethod INT)") # these sequences contain indels

    # When sequences get trimmed, rather than store a new sequence, define a SiteSet mapping from->to site
    # regions.
    cur.execute("create table if not exists SiteSets(id INTEGER primary key autoincrement, setname TEXT)") # sets of sequence sites
    
    # Defines a set of sites, 'from' site to the 'to' site, inclusive. A set/sequence pair can have
    # multiple listings in this table, for example, if the desired sites are non-contiguous
    cur.execute("create table if not exists SiteSetsAlignment(setid INTEGER, almethod INT, fromsite INT, tosite INT)") # alseqid is the ID of an aligned sequence
    
    # for ZORRO:
    cur.execute("create table if not exists AlignmentSiteScoringMethods(id INTEGER primary key autoincrement, name TEXT)")
    cur.execute("create table if not exists AlignmentSiteScores(almethodid INTEGER, scoringmethodid INTEGER, site INT, score FLOAT)")
    cur.execute("create table if not exists SampledSequences(id INTEGER primary key autoincrement, taxonid INT, samsequence TEXT, datatype INT, almethod INT)")
    
    # a list of phylo softwares, such as PhyML and RAxML
    cur.execute("create table if not exists PhyloSoftwares(id INTEGER primary key autoincrement, name TEXT, exe_path TEXT)") 
    # a list of phylogenetic models, such as LG+I+G
    cur.execute("create table if not exists PhyloModels(modelid INTEGER primary key autoincrement, name TEXT)") 
    # command-line strings used to invoke phylo models in different softwares
    cur.execute("create table if not exists PhyloModelsSoftware(modelid INTEGER, softwareid INTEGER, software_model_string TEXT)") 
    cur.execute("create table if not exists ParsimonyPhylogenies(id INTEGER primary key autoincrement, seqsetid INTEGER, newick TEXT)") # seqsetid is the site set used to make this tree
    cur.execute("create table if not exists UnsupportedMlPhylogenies(id INTEGER primary key autoincrement, phylomodelid INTEGER, seqsetid INTEGER, newick TEXT)")

    # Phyml and Mr. Bayes:
    
    # a list of branch support methods, such as aLRT and PP
    cur.execute("create table if not exists BranchSupportMethods(id INTEGER primary key autoincrement, name TEXT)")
    # maps unsupported trees to a support method
    cur.execute("create table if not exists SupportedMlPhylogenies(id INTEGER primary key autoincrement, unsupportedmltreeid INTEGER, newick TEXT, supportmethodid INT)")
    
    # define groups of taxa
    cur.execute("create table if not exists Ingroups(id INTEGER primary key autoincrement, name TEXT)")
    cur.execute("create table if not exists Outgroups(id INTEGER primary key autoincrement, name TEXT)")
    cur.execute("create table if not exists GroupsTaxa(groupid INTEGER, taxonid INTEGER)") # groupid is the ID of either an ingroup or outgroup
    
    cur.execute("create table if not exists Ancestors(id INTEGER primary key autoincrement, ingroupid INTEGER, outgroupid INTEGER, name TEXT)")
    cur.execute("create table if not exists AncestorsTreeNumbers(ancestorid INTEGER, unsupportedmltreeid INTEGER, number INT)")
    cur.execute("create table if not exists AncestralCladogram(id INTEGER primary key autoincrement, unsupportedmltreeid INTEGER, newick TEXT)")
    
    """For HTML visualization"""
    cur.execute("create table if not exists HtmlPieces(id INTEGER primary key autoincrement, keyword TEXT, htmlstring TEXT)")
    con.commit()
    
    return con