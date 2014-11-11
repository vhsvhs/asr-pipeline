from tools import *


def clean_seqs( apath ):
    fin = open(apath, "r")
    
    taxa = []
    
    for l in fin.readlines():
        if l.startswith(">"):
            tokens = l.split()
            # The name is the first token on the line.
            name = tokens[0]
            name = re.sub("_", ".", name)
            if name.__contains__("|"):
                name = name.split("|")[0] + "." + name.split("|")[2]
            print name
            if name not in taxa:
                taxa.append( name )
            else:
                print "Your sequence database", sys.argv[1], "contains two (or more) sequences with the name", name
                print "Please fix this problem before proceeding."
                exit()
        else:
            l = l.strip()
            l = re.sub("\*", "", l)
            print l
    fin.close()

#print "I found a total of", taxa.__len__(), "sequences."

#clean_seqs.py