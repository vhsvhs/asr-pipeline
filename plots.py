import math, os, re, sys

from scipy import stats as scipystats

def calculate_bins(values):    
    maxval = max(values)
    minval = min(values)
    
    MAX_BIN_LENGTH = float("%.5f"%maxval)
    MIN_BIN_LENGTH = float("%.5f"%minval)
    BIN_SIZE = (MAX_BIN_LENGTH - MIN_BIN_LENGTH) / 20
    MAX_BIN_LENGTH += BIN_SIZE
    MIN_BIN_LENGTH -= BIN_SIZE
    
    bins = {} # key = floor, value = count
    bins[MAX_BIN_LENGTH] = 0.0    
    
    for b in values:        
        if b > MAX_BIN_LENGTH:
            bins[MAX_BIN_LENGTH] += 1
        else:    
            bin = MIN_BIN_LENGTH
            while (bin <= b):
                if bin not in bins:
                    bins[bin] = 0.0
                if bin+BIN_SIZE > b and bin <= b:
                    bins[bin] += 1
                bin += BIN_SIZE
                
    normalized_bins = {}
    for b in bins.keys():
        normalized_bins[b] = float(bins[b]) / values.__len__()
    return normalized_bins


def histogram(data, filekeyword, xlab="score", ylab="proportion"):
    """Plots a histogram of the values in the data
    data = a list of values, unsorted
    """
    
    bins = calculate_bins(data)
    
    pointsets = bins.keys()
    pointsets.sort()
        
    finalset = pointsets[ pointsets.__len__()-1 ]
    tablepath = filekeyword + ".tmp"
    fout = open(tablepath, "w")
    for p in pointsets:
        if p != finalset:
            fout.write(p.__str__() + "\t")
        else:
            fout.write(p.__str__() )
    fout.write("\n")
    for p in pointsets:
        if p != finalset:
            fout.write( bins[p].__str__() + "\t")
        else:
            fout.write( bins[p].__str__() )            
    fout.write("\n")
    fout.close()
    
    pdfpath = filekeyword + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=6, height=4);\n"    
    cranstr += "bars <- read.table(\"" + tablepath + "\", header=T, sep=\"\\t\")\n"

    cranstr += "pointsets <- c("
    for p in pointsets:
        cranstr += (p).__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "barx = barplot(as.matrix(bars), beside=TRUE, ylim=range(0,1.0), names.arg=pointsets, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\");\n"
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)



def scatter1(values_ii, values_jj, filekeyword, xlab="", ylab="", force_square=False):    
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    
    pdfpath = filekeyword + ".pdf"
    cranstr += "pdf(\"" + pdfpath + "\", width=6, height=6);\n"    

    print "\n. Writing a scatterplot to", pdfpath

    # X values
    cranstr += "x<-c("
    for v in values_ii:
        cranstr += v.__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"

    # Y values
    cranstr += "y<-c("
    for v in values_jj:
        cranstr += v.__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "plot(x, y, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\""

    maxx = max(values_ii)
    maxy = max(values_jj)
    lim = max( [maxx, maxy] )    
    if force_square:
        cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
    
    cranstr += ");\n"
    
    if force_square:
        cranstr += "abline(0,1)\n"
    
    """Pearson's linear value correlation."""
    corr_valsa = []
    corr_valsb = []
    for ww in range(0, values_ii.__len__()):
        if values_ii[ww] != 0 and values_jj[ww] != 0:
            corr_valsa.append( values_ii[ww] )
            corr_valsb.append( values_jj[ww] )
    (rho, pvalue) = scipystats.pearsonr( corr_valsa, corr_valsb )
    cranstr += "text(" + ((lim-min(values_ii))/2).__str__() + ", " + lim.__str__() + ", \"Prs R=%.2f"%rho + ", P=%.2f"%pvalue + "\");\n"
    
    """Spearman's non-linear non-parametric rank correlation."""
    (rho, pvalue) = scipystats.spearmanr( corr_valsa, corr_valsb )
    cranstr += "text(" + ((lim-min(values_ii))/2).__str__() + ", " + (0.92*lim).__str__() + ", \"Spr R=%.2f"%rho + ", P=%.2f"%pvalue + "\");\n"
    
    cranstr += "dev.off();\n"
    
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save --slave < " + cranpath)
    
    return cranpath