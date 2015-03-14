from dendropy import Tree

from asrpipelinedb_api import *
from tools import *

def index_mutations(con):
    """Builds an index of all mutations"""
    cur = con.cursor()
    
    for msaid in get_alignment_method_ids(con):
        for modelid in get_phylo_modelids(con):
            newick = get_anc_cladogram(con, msaid, modelid)
            t = Tree()
            t.read_from_string(newick, "newick")
            for edge in t.preorder_edge_iter():
                if edge.head_node == None or edge.tail_node == None:
                    continue
                if edge.head_node.label == None or edge.tail_node.label == None:
                    continue 
                print msaid, modelid, edge.head_node.label, edge.tail_node.label
                anc1name = "Node" + edge.head_node.label.__str__()
                anc2name = "Node" + edge.tail_node.label.__str__()
                index_mutations_helper(con, msaid, modelid, anc1name, anc2name)
              
                    
def index_mutations_helper(con, msaid, phylomodelid, anc1name, anc2name):
    
    """To-do continue here -- there may be some redundancy here.
        I.e., this method is being called for each msa/model, but then later in the method
        we iterate over msa/models AGAIN."""
    
    cur = con.cursor()
    
    sql = "select name from AlignmentMethods where id=" + msaid.__str__()
    cur.execute(sql)
    msaname = cur.fetchone()[0]
    
    sql = "select name from PhyloModels where modelid=" + phylomodelid.__str__()
    cur.execute(sql)
    phylomodelname = cur.fetchone()[0]

    """Get the seedsequence with indels"""
    seedsequence = get_seed_sequence(con, msaid)
    nsites = seedsequence.__len__()
    sql = "select id, shortname from Taxa where shortname in (select value from Settings where keyword='seedtaxa')"
    cur.execute(sql)
    x = cur.fetchone()
    seedtaxonid = x[0]
    seedtaxonname = x[1]
    sl = seedtaxonname.__len__()
    if sl > 10:
        seedtaxonnameshort = seedtaxonname[0:4] + "..." + seedtaxonname[sl-7:sl]
    else:
        seedtaxonnameshort = seedtaxonname
    
    """Which ancestors are selected (used to define the branch)"""
    ancid1 = None
    ancid2 = None
    sql = "select name, id from Ancestors where name='" + anc1name + "' and almethod=" + msaid.__str__() + " and phylomodel=" + phylomodelid.__str__()
    cur.execute(sql)
    x = cur.fetchone()
    if x != None:
        ancid1 = x[1]
    sql = "select name, id from Ancestors where name='" + anc2name + "' and almethod=" + msaid.__str__() + " and phylomodel=" + phylomodelid.__str__()
    cur.execute(sql)
    x = cur.fetchone()
    if x != None:
        ancid2 = x[1]
        
    sql = "select name from Ancestors where id=" + ancid1.__str__()
    cur.execute(sql)
    ancname1 = cur.fetchone()[0]

    sql = "select name from Ancestors where id=" + ancid2.__str__()
    cur.execute(sql)
    ancname2 = cur.fetchone()[0]
    
    """Get the site map between different sequence alignments"""
    msa_site1_site2 = {} # key = msaid, value = hash; key = site in user-specified msa, value = site in msaid
    sql = "select id, name from AlignmentMethods"
    cur.execute(sql)
    msaid_name = {}
    for ii in cur.fetchall():
        this_msaid = int(ii[0])
        this_msaname = ii[1].__str__()
        msaid_name[ this_msaid ] = this_msaname
        msa_site1_site2[this_msaid] = {}
        
        """For each alignment method, get the mapped-to sites from the SiteMap table."""
        sql = "select site1, site2 from SiteMap where"
        sql += " almethod1=" + msaid.__str__() # the user-specified alignment method
        sql += " and almethod2=" + ii[0].__str__()
        sql += " and taxonid=" + seedtaxonid.__str__()
        cur.execute(sql)
        query = cur.fetchall()
        for qq in query:
            site1 = qq[0]
            site2 = qq[1]
            msa_site1_site2[this_msaid][site1] = site2  
                
    """Note: at this point,  msa_site1_site2 contains data about all alignments EXCEPT
        for the user-specified msaid"""

    """ Get a list of these ancestors in other alignments and models """ 
    matched_ancestors = get_ancestral_matches(con, ancid1, ancid2)
    matched_ancestors = [ (ancid1,ancid2) ] + matched_ancestors
     
    ancid_msaid = {}
    ancid_model = {}
    ancid_site_state_pp = {}
    ancid_name = {}    
    ancid_site_statepp = {}
    for match in matched_ancestors:
        for id in [ match[0], match[1] ]:
            sql = "Select almethod, phylomodel, name from Ancestors where id=" + id.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            ancid_msaid[id] = int(x[0])
            ancid_model[id] = int(x[1])
            ancname = x[2]
            ancid_name[id] = ancname
            #msaname = msaid_name[ ancid_msaid[id] ]
            ancid_site_statepp[id] = get_site_ml(con, id, skip_indels = False)   

    """We need model names and alignment names for printing in the mutation table header rows."""
    phylomodelid_name = {}
    sql = "select modelid, name from PhyloModels"
    cur.execute(sql)
    query = cur.fetchall()
    for ii in query:
        phylomodelid_name[   ii[0]  ] = ii[1]
               
    mutation_header = []
    mutation_rows = []
    mutation_shortlist = [] # each item is a tuple: (site in msa, site in seed, ancestral state 1, anc. state 2, weight)
    seed_site = 0 # the index into the seed sequence sans indels
    """For each site in the seed sequence"""
    for site in range(1, nsites+1):
        found_content = False # did we find any non-indel characters?
        this_row = []
        
        if site == 1:
            """Add header information"""
            mutation_header = ["Site in " + msaname]
            mutation_header.append( "Site in\n" + seedtaxonnameshort )

        if seedsequence[site-1] != "-":
            seed_site += 1
        
        count_replicates = 0 # how many of the columns support this mutation?
         
        """For each ancestor"""
        for match in matched_ancestors:
            #print "988:", match
            this_anc1 = match[0]
            this_anc2 = match[1]
            
            this_msaid1 = ancid_msaid[this_anc1]
            this_msaid2 = ancid_msaid[this_anc2]
            
            this_msaname1 = msaid_name[this_msaid1]
            this_msaname2 = msaid_name[this_msaid2]
            
            this_modelid1 = ancid_model[this_anc1]
            this_modelid2 = ancid_model[this_anc2]
            
            this_modelname1 = phylomodelid_name[this_modelid1]
            this_modelname2 = phylomodelid_name[this_modelid2]
            
            if this_modelid1 != this_modelid2:
                print "954: error"
                continue
            if this_msaid1 != this_msaid2:
                print "957: error"
                continue

            if site == 1:
                """Add header information"""
                mutation_header.append( this_modelname1 + "\n" + this_msaname1 )
  
            site2 = None
            if int(this_msaid1) == int(msaid):
                site2 = site
            elif seedsequence[site-1] != "-":
                """Get a mapped site"""
                if site in msa_site1_site2[this_msaid1]:
                    site2 = msa_site1_site2[this_msaid1][site]

            this_column = ("","","","","", "") # Six values for the Django template to use. see the line down below starting with this_column = (
            if site2 != None:
                (anc1state, anc1pp) = ancid_site_statepp[this_anc1][site2]
                (anc2state, anc2pp) = ancid_site_statepp[this_anc2][site2]
                mutation_flag = ""
                if anc1state != "-" or anc2state != "-":
                    found_content = True
                if anc1state == "-" and anc2state != "-":
                    mutation_flag = "i"
                elif anc2state == "-" and anc1state != "-":
                    mutation_flag = "d"
                elif anc1state != anc2state and anc1pp > 0.7 and anc2pp > 0.7:
                    mutation_flag = "1"
                elif anc1state != anc2state and anc1pp > 0.5 and anc2pp > 0.5:
                    mutation_flag = "2"
                elif anc1state != anc2state and anc1pp > 0.4 and anc2pp > 0.4:
                    mutation_flag = "3"
                    
                if mutation_flag in ["1", "2", "3"]:
                    count_replicates += 1
                
                sql = "insert into MutationIndex (ancid1, ancid2, mlstate1, mlstate2, state1pp, state2pp)"
                sql += " values("
                sql += this_anc1.__str__() + ","
                sql += this_anc2.__str__() + ","
                sql += "'" + anc1state.__str__() + "',"
                sql += "'" + anc2state.__str__() + "',"
                sql += anc1pp.__str__() + ","
                sql += anc2pp.__str__()
                sql += ")"
                cur.execute(sql)
            con.commit()