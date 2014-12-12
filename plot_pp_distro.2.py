import os, sys, re

from tools import *

#dirs = ["MUSCLE/asr.muscle.PROTGAMMALG", "MUSCLE/asr.muscle.PROTCATLG", "MUSCLE/asr.muscle.PROTGAMMAWAG", "MUSCLE/asr.muscle.PROTCATWAG", "MSAPROBS/asr.msaprobs.PROTGAMMALG", "MSAPROBS/asr.msaprobs.PROTGAMMAWAG", "MSAPROBS/asr.msaprobs.PROTCATLG", "MSAPROBS/asr.msaprobs.PROTCATWAG"]

#
# These are the sites in the known catalyic domain
# . . . i.e., the sites we really care about.
#
#msaprobs_struct_sites = "1352 3316"
#muscle_struct_sites = "749 2082"

path_mean = {}
path_seq = {}
path_sd = {}

#
# depricated
#
# def parse_summary(path, keyword):
#     fin = open(path, "r")
#     lines = fin.readlines()
#     last_was_cran = False
#     for l in lines:
#         l = l.strip()
#         if l.startswith(">"):
#             last_was_cran = True
#         elif last_was_cran:
#             path_seq[ keyword ] = l
#             last_was_cran = False
# 
#         if l.__contains__("mean PP ="):
#             path_mean[ keyword ] = l.split()[3]
#         if l.__contains__("standard deviation ="):
#             path_sd[ keyword ] = l.split()[3]
#     fin.close()
# 
# 
# for d in dirs:
#     msapath = get_fullphylippath(d)
#     [start,stop] = get_boundary_sites(msapath,"Saccharomyces.cerevisiae.IME2")
# 
#     for f in os.listdir(d + "/tree1"):
#         if f.startswith("node") and f.endswith(".dat"):
#             short = re.sub(".dat", "", f)
#             print dir, f
# 
#             #
#             # Analyze full length protein sequences
#             #
#             c = "python ~/Documents/SourceCode/Lazarus/plot_pp_distribution.py " + d + "/tree1/" + f + " > " + d + "/tree1/" + short + ".summary.txt"
#             print c
#             os.system(c)
# 
#             parse_summary(d + "/tree1/" + short + ".summary.txt", d + "/" + f)
#             os.system("mv barplot.anc.pdf " + d + "/tree1/" + short + ".barplot.pdf")
#             os.system("mv barplot.anc.cran " + d + "/tree1/" + short + ".barplot.cran")
#             os.system("mv barplot.table.anc.txt " + d + "/tree1/" + short + ".pptable.txt")
# 
#             #
#             # Kinase domain structural sites:
#             #
# #             os.system("rm -rf " + d + "/tree1/*.32-387.*")
# #             if d.__contains__("msaprobs"):
# #                 os.system("python ~/Documents/SourceCode/Lazarus/plot_pp_distribution.py " + d + "/tree1/" + f + " " + msaprobs_struct_sites + " > " + d + "/tree1/" + short + ".struct_only.summary.txt")
# #             else:
# #                 os.system("python ~/Documents/SourceCode/Lazarus/plot_pp_distribution.py " + d + "/tree1/" + f + " " + muscle_struct_sites + " > " + d + "/tree1/" + short + ".struct_only.summary.txt")
# #             parse_summary(d + "/tree1/" + short + ".struct_only.summary.txt", d + "/" + f + "-struct_only")
# #             os.system("mv barplot.anc.pdf " + d + "/tree1/" + short + ".struct_only.barplot.pdf")
# #             os.system("mv barplot.anc.cran " + d + "/tree1/" + short + ".struct_only.barplot.cran")
# #             os.system("mv barplot.table.anc.txt " + d + "/tree1/" + short + ".struct_only.pptable.txt")
# 
#             """
#             #
#             # Only those sites in alpha helices and beta sheets. . .
#             #
#             if d.__contains__("msaprobs"):
#                 os.system("python ~/Documents/SourceCode/Lazarus/plot_pp_distribution.py " + d + "/tree1/" + f + " " + msaprobs_ab_sites + " > " + d + "/tree1/" + short + ".ab-only.summary.txt")
#             else:
#                 os.system("python ~/Documents/SourceCode/Lazarus/plot_pp_distribution.py " + d + "/tree1/" + f + " " + muscle_ab_sites + " > " + d + "/tree1/" + short + ".ab-only.summary.txt")
#             parse_summary(d + "/tree1/" + short + ".ab-only.summary.txt", d + "/" + f + "-ab-only")
#             os.system("mv barplot.anc.pdf " + d + "/tree1/" + short + ".ab-only.barplot.pdf")
#             os.system("mv barplot.anc.cran " + d + "/tree1/" + short + ".ab-only.barplot.cran")
#             os.system("mv barplot.table.anc.txt " + d + "/tree1/" + short + ".ab-only.pptable.txt")
#             """

# hash: key = dat file path, value = tuple
node_full = {}
node_struct = {}
node_ab = {}

for path in path_mean:
    if path.endswith(".dat"):
        node_full[ path ] = ( path_mean[path], path_sd[path], path_seq[path] )
    elif path.endswith("ab-only"):
        node = re.sub("-ab-only", "", path)
        node_ab[ node ] = ( path_mean[path], path_sd[path], path_seq[path] )
    #elif path.endswith("-struct_only"):
    #    node = re.sub("-struct_only", "", path)
    #    node_struct[ node ] =( ( path_mean[path], path_sd[path], path_seq[path] ) )

for node in node_full:
    line = node
    line += "\t" + node_full[node][0]
    line += "\t" + node_full[node][1]

    line += "\t" + node_struct[node][0]
    line += "\t" + node_struct[node][1]

    #line += "\t" + node_ab[node][0]
    #line += "\t" + node_ab[node][1]

    line += "\t" + path_seq[node]
    line += "\t" + path_seq[node + "-struct_only"]
    line += "\t" + path_seq[node + "-ab-only"]
    print line

for path in path_mean:
    print path, path_mean[path], path_sd[path], path_seq[path]
