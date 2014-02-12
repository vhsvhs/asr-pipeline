import os, re, sys, dendropy
from Bio import Phylo
from tools import *

def annotate_phyloxml(line, ancdir):
    """Input a line of PhyloXML containing a cladogram with ancestral node numbers, outputs that line with additional markup (possibly) added.
    Specifically, changes <confidence> to <name>, add adds <annotations> with a clickable URL
    to the ancestral PP file."""
    line = re.sub("confidence", "name", line)
    line = re.sub("type=\"unknown\"", "", line)

    # The following line breaks jsphylosvg, for unknown reasons.
    # apparently, the whitespace is important in their javascript.
    #line = re.sub(" ", "", line)

    if line.__contains__("name"):
        name = line
        name = re.sub("\n", "", name)
        name = re.sub(" ", "", name)
        name = re.sub("<name>", "", name)
        name = re.sub("<\/name>", "", name)
        if name.isdigit(): # this taxa name is all numbers, and therefore, an ancestor
            line += "<annotation>"
            line += "<desc>Ancestor #" + name + "</desc>"
            url = ancdir + "/node" + name + ".html"
            line += "<uri>" + url + "</uri>"
            line += "</annotation>"
    #print "Retline:", line
    return line


def newick_to_xml(dir, model):
    npath = get_cladogram_path(dir, model)
    ancdir = "asr." + get_runid(dir, model)

    """Builds PhyloXML version of newick tree at npath."""
    tree = Phylo.read(npath,'newick')
    xmlpath = re.sub(".tre", ".xml", npath)
    Phylo.write(tree, xmlpath, 'phyloxml')
    fin = open(xmlpath, "r")
    o = ""
    for l in fin.readlines():
        l = annotate_phyloxml(l, ancdir)
        o += l
    o = re.sub("\n", "", o)
    #o = re.sub(" ", "", o)
    return (xmlpath, o)