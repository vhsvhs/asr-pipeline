import math,re,sys,os

def alrt_to_alr( treepath, outpath ):
    fin = open(treepath, "r")
    newline = ""
    for line in fin.readlines():
        #print line
        tokens = line.split(":")
        for t in tokens:
            if t.__contains__(")") and False == t.__contains__(";"):
                #print t
                ts = t.split(")")
                newline += ts[0] + ")"
                alrt = float(ts[1])
                if alrt > 1418:
                    alrt = 1418.0
                alr = math.exp(alrt/2.0)

                """Wrap the ALR value in single quotes '...'
                    because the scientific notation causes problems
                    for Newick viewers later in this pipeline."""
                newline += "'%.4e'"%alr
                newline += ":"
            else:
                newline += t
                newline += ":"
    fin.close()
    
    fout = open(outpath, "w")
    fout.write( newline + "\n")
    fout.close()
