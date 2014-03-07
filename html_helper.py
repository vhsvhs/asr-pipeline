import os, sys, time
from configuration import *
from tools import *
from phyloxml_helper import *

HTMLDIR = "HTML"

def get_header(head = "", urlpre = ""):
    """Writes the top required HTML lines for a report webpage.  head = extra lines to stick in the <head>,
    urlpre = the directory containing HTML, or by default './' """
    h = ""
    h += "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n"
    h += "<html>\n"  
    h += "<head>\n"
    h += "<title>" + ap.params["project_title"] + " Ancestral Library</title>\n"
    h += "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n"
    h += "<link rel=\"stylesheet\" href=\"" + urlpre + "asrpipeline.css\">\n"
    h += head
    h += "</head>\n"    
    h += "<body><br>\n"
    if "HTML_SPECIAL1" in ap.params:
        h += ap.params["HTML_SPECIAL1"] + "\n"
    h += "<h1>" + ap.params["project_title"] + " Ancestral Library</h1>\n"    
    h += "<hr>\n"
    h += "<p><a href='" + urlpre + "index.html'>Overview</a>"
    h += " | <a href='" + urlpre + "alignments.html'>Alignments</a>" 
    h += " | <a href='" + urlpre + "trees.html'>Trees & Ancestors</a>"
    if os.path.exists("SIMULATION"):
        h += " | <a href='" + urlpre + "errorsimulation.html'>Accuracy Assessment</a>"
    if os.path.exists("SCRIPTS/compareanc_commands.sh"):
        h += " | <a href='" + urlpre + "anccomp.html'>Mutations</a>"    
    #
    # to-do: bayesian samples, with an interactive widget
    #
    h += "</p>\n"
    h += "<hr>\n"
    return h

def get_footer():
    f = "\n"
    f += "<hr>\n"
    f += "<p class='smalltext' align='right'>Last Updated "
    now = time.strftime("%c")
    f += "%s"%now 
    f += "<br>\n"
    f += "Questions/Comments? Email: victorhansonsmith@gmail.com"
    f += "</p>"
    f += "</body>\n"
    f += "</html>\n"
    return f

def write_index():
    fout = open(HTMLDIR + "/index.html", "w") 
    fout.write( get_header() )
    
    fout.write("<p>This library was calculated with the following set of curated protein sequences:</p>\n")
    fout.write("<p>Download: <a href=\"../" + ap.params["ergseqpath"] + "\">fasta</a></p>\n")
    
    fout.write( get_footer() )

    fout.close()

def write_alignments():
    out = ""
    out += get_header()

    out += "<table width=100%>\n"
    out += "<tr>"
    out += "<th align='left'>Alignment Method</th>"
    out += "<th align='center'>N taxa</th>"
    out += "<th align='center'>N sites</th>"
    out += "<th align='center'>Download</th>"
    out += "</tr>\n"
    
    for d in ap.params["msa_algorithms"]:
        out += "<tr>"
        out += "<td aling='left'><p>" + d 
        if d == "msaprobs":
            out += " <a class='smalltext' href=\"http://msaprobs.sourceforge.net/\">[ref]</a>"
        if d == "muscle":
            out += " <a class='smalltext' href=\"https://www.ebi.ac.uk/Tools/msa/muscle/\">[ref]</a>"
        out += "</p></td>"
        fpath = get_fastapath(d)
        ppath = get_phylippath(d)
        (ntaxa, nsites) = get_phylipstats(ppath)
        out += "<td align='center'>" + ntaxa.__str__() + "</td>"
        out += "<td align='center'>" + nsites.__str__() + "</td>"
        out += "<td align='center'><a href='../" + fpath + "'>fasta</a> | <a href='../" + ppath + "'>phylip</a></td>"
        out += "</tr>\n"
    out += "</table>\n"

    out += get_footer()

    fout = open(HTMLDIR + "/alignments.html", "w") 
    fout.write(out)
    fout.close()


def read_lnllog(dir):
    fin = open(dir + "/raxml.lnl.summary.txt", "r")
    lines = fin.readlines()
    fin.close()
    model_data = {}
    for l in lines:
        if l.__len__() > 3:
            tokens = l.split()
            model = tokens[0].split(".")[1]
            lnl = "%.3f"%float(tokens[1])
            pp = "%.3f"%float(tokens[2])
            model_data[model] = (lnl,pp)
    return model_data

def write_treesancs():
    out = ""
    out += get_header()

    for d in ap.params["msa_algorithms"]:
        out += "<table width=100%>\n"
        out += "<tr>"
        out += "<th align='left'>Alignment</th>"
        out += "<th align='left'>Model</th>"
        out += "<th align='right'>log(Likelihood)</th>"
        out += "<th>Relative Probability</th>"
        out += "<th>&sum; Branch Lengths</th>"
        out += "<th>Download</th>"
        out += "<th>Ancestors</th>"
        out += "</tr>\n"

        model_data = read_lnllog(d)
        for model in ap.params["raxml_models"]:
            tpath = get_alrt_treepath(d, model)
            tlength = get_tree_length( tpath )
            out += "<tr>"
            out += "<td>" + d + "</td>"
            out += "<td>" + model + "</td>"
            out += "<td align='right'>" + model_data[model][0].__str__() + "</td>"
            out += "<td align='center'>" + model_data[model][1].__str__() + "</td>"
            out += "<td align='center'>" + "%.3f"%tlength + "</td>"
            out += "<td align='center'><a href='../" + tpath + "'>newick</a></td>"
            out += "<td align='center'><a href='tree." + d + "." + model + ".html'>view</a></td>"
            out += "</tr>\n"
            write_anctree(d, model)
        out += "</table>\n"

    out += "<br>"
    out += "<p class='smalltext'>***More information about these phylogenetic models can found in the RAxML Software Manual <a href='http://www.trex.uqam.ca/documents/RAxML-Manual.7.0.3.pdf'>[ref]</a>.</p>"

    out += get_footer()

    fout = open(HTMLDIR + "/trees.html", "w")
    fout.write( out )
    fout.close()


def write_anctree(d, model):
    rid = get_runid(d, model)
    tpath = get_alrt_treepath(d, model)

    #
    # to-do: copy the contents of HTML_SUPPORT to HTMLDIR. Otherwise, the necessary Javascript will be missing.
    #

    js = "<script type=\"text/javascript\" src=\"HTML_SUPPORT/raphael-min.js\" ></script>\n"
    js += "<script type=\"text/javascript\" src=\"HTML_SUPPORT/jsphylosvg-min.js\"></script>\n"

    js += "<script type=\"text/javascript\">"
    js += "window.onload = function(){\n"

    npath = get_cladogram_path(d, model)
    (xmlpath, xmlstring) = newick_to_xml(d, model)

    js += "       var dataObject = { phyloxml: "
    js += "'" + xmlstring + "' ,"
    js += "fileSource:false};\n"

    width = 800
    (ntaxa,nsites) = get_phylipstats( get_phylippath(d) )
    height = ntaxa * 16
    js += "\nphylocanvas = new Smits.PhyloCanvas(dataObject,\"svgCanvas\","
    js += "800," + height.__str__() + ");\n"
    js += "};\n" # closes onload

    """
    js += "Smits.PhyloCanvas.Render.Style = {\n"
    js += "   line: {\n"
    js += "       \"stroke\":'rgb(0,0,0)',\n"
    js += "       \"stroke-width\": 1.5\n"
    js += "   },\n"
    js += "   text: {\n"
    js += "       \"font-family\":'Helvetica',\n"
    js += "       \"font-size\": 10\n"
    js += "   },\n"
    js += "}\n"
    """
    js += "</script>"

    out = ""
    out += get_header(head=js)
    out += "<h2>Cladogram of Ancestors, Alignment: " + DIR_nick[d] + ", Model: " + model + "</h2>\n"
    out += "<h3>Click an ancestral node for details.</h3>\n"
    out += "<div id=\"svgCanvas\"> </div>\n"
    out += get_footer()
    fout = open(HTMLDIR + "/tree." + d + "." + model + ".html", "w")
    fout.write( out )
    fout.close()


def get_color_for_pp(pp):
    style = ""
    if pp >= 1.0:
        style = "#1E90FF"
    elif pp >= 0.9:
        style = "#00BFFF"
    elif pp >= 0.8:
        style = "#90EE90"
    elif pp >= 0.7:
        style = "#ADFF2F"
    elif pp >= 0.6:
        style = "#FFD700"
    elif pp >= 0.5:
        style = "#FFA07A"
    else:
        style = "#FFB6C1"
    return style

def get_style_for_pp(pp):
    style = ""
    if pp >= 1.0:
        style = "ppfull"
    elif pp >= 0.9:
        style = "pp9"
    elif pp >= 0.8:
        style = "pp8"
    elif pp >= 0.7:
        style = "pp7"
    elif pp >= 0.6:
        style = "pp6"
    elif pp >= 0.5:
        style = "pp5"
    else:
        style = "pplow"
    return style

def pp_distro_to_html(data):
    """Input: a hashtable of site-PP values, gotten from get_pp_distro().
    Output: a nice HTML table."""
    out = "<table>\n"
    out += "<tr class='ppheader'>"
    out += "<td align='center'><strong>Site</strong></td>"
    out += "<td colspan=20 align='left'><strong>Support</strong></td>"
    out += "</tr>\n"
    for site in data:
        if False == data[site][0].__contains__("-"):
            """Only write rows for sites with non-indels."""
            out += "<tr>"
            out += "<td>" + site.__str__() + "</td>"
            for ss in data[site]:
                style = get_style_for_pp( ss[1] )
                out += "<td class='" + style + "'>" + ss[0] + " " + ss[1].__str__() + "</td>"
            for i in range(data[site].__len__(),21):
                out += "<td></td>"
            out += "</tr>"
    
    out += "</table>\n"
    return out

def get_ppdistro_key():
    out = "<table>"
    out += "<tr><td colspan=8 align='left' class='ppheader'><strong>Legend</strong></td></tr>\n"
    out += "<tr>"
    out += "<td cellpadding=4>Probability:</td>"
    out += "<td class='ppfull'>1.0</td>"
    out += "<td class='pp9'>(1.0, 0.9]</td>"
    out += "<td class='pp8'>(0.9,0.8]</td>"
    out += "<td class='pp7'>(0.8,0.7]</td>"
    out += "<td class='pp6'>(0.7,0.6]</td>"
    out += "<td class='pp5'>(0.6,0.5]</td>"
    out += "<td class='pplow'>(0.5,0.0]</td>"
    out += "</tr>\n"
    return out

def get_ppdistro_summary(data):
    pps = []
    for site in data:
        if data[site][0][0] != '-' and data[site][0][0] != "-":
            pps.append( data[site][0][1] )
    mean = get_mean(pps)
    sd = get_sd(pps)
    return "<p>&#x3BC; PP = %.3f"%get_mean(pps) + ", &#x3C3; = %.3f"%get_sd(pps) + "</p>\n"


def write_ppdistro_plot(data):
    """Writes a JS script tag for a Google barplot"""
    keys = ["[0.0,0.1)","[0.1,0.2)","[0.2,0.3)","[0.3,0.4)","[0.4,0.5)","[0.5,0.6)","[0.6,0.7)","[0.7,0.8)","[0.8,0.9)","[0.9,1.0)", "1.0"]
    key_prop = {}
    for k in keys:
        key_prop[k] = 0.0
    nsites = 0
    for site in data:
        if data[site][0][0] != '-' and data[site][0][0] != "-":
            nsites += 1
            p = data[site][0][1]
            if p < 0.1:
                key_prop["[0.0,0.1)"] += 1
            if p < 0.2:
                key_prop["[0.1,0.2)"] += 1
            elif p < 0.3:
                key_prop["[0.2,0.3)"] += 1
            elif p < 0.4:
                key_prop["[0.3,0.4)"] += 1
            elif p < 0.5:
                key_prop["[0.4,0.5)"] += 1
            elif p < 0.6:
                key_prop["[0.5,0.6)"] += 1
            elif p < 0.7:
                key_prop["[0.6,0.7)"] += 1
            elif p < 0.8:
                key_prop["[0.7,0.8)"] += 1
            elif p < 0.9:
                key_prop["[0.8,0.9)"] += 1
            elif p < 1.0:
                key_prop["[0.9,1.0)"] += 1
            else:
                key_prop["1.0"] += 1
    for k in key_prop:
        key_prop[k] = float(key_prop[k]) / nsites

    out = ""
    out += "<script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>"
    out += "<script type=\"text/javascript\">\n"
    out += "google.load(\"visualization\", \"1\", {packages:[\"corechart\"]});\n"
    out += "google.setOnLoadCallback(drawChart);\n"
    out += "function drawChart() {\n"
    out += "var data = google.visualization.arrayToDataTable([\n"
    out += "['Probability', 'Proportion of Sites', {role:'style'}],\n"
    for k in keys:
        style = "'black'"
        pp = 0.0
        if k.__contains__("0.1)"):        
            pp = 0.05
        elif k.__contains__("0.2)"):        
            pp = 0.15
        elif k.__contains__("0.3)"):        
            pp = 0.25
        elif k.__contains__("0.4)"):        
            pp = 0.35
        elif k.__contains__("0.5)"):        
            pp = 0.45
        elif k.__contains__("0.6)"):        
            pp = 0.55
        elif k.__contains__("0.7)"):        
            pp = 0.65
        elif k.__contains__("0.8)"):        
            pp = 0.75
        elif k.__contains__("0.9)"):        
            pp = 0.85
        elif k.__contains__("1.0)"):        
            pp = 0.95
        elif k.__contains__("1.0"):        
            pp = 1.0
        style = get_color_for_pp( pp )

        out += "['" + k + "', " + key_prop[k].__str__() + ",'" + style + "']"
        if k != "1.0":
            out += ",\n"
        else:
            out += "\n"
    out += "]);\n"
    out += "var options = {\n"
    out += "  title: '',\n"
    out += "  vAxis: {title: 'Proportion of Sites',  titleTextStyle: {color: 'black', fontName:'Helvetica',fontSize:'12pt'}},\n"
    out += "  chartArea:{left:60,top:10,width:\"100%\",height:\"75%\"},\n"
    out += "  orientation: 'horizontal',\n"
    out += "  backgroundColor: '#f0f0f0',\n"
    out += " legend: { position: \"none\" }\n"
    out += "};\n"
    out += "var chart = new google.visualization.BarChart(document.getElementById('chart_div'));\n"
    out += "chart.draw(data, options);\n"
    out += "}\n"
    out += "</script>\n"
    return out


def write_ancestors_indi():
    """Writes on HTML page for each ancestor."""
    for d in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            outdir = HTMLDIR + "/asr." + get_runid(d, model)
            os.system("mkdir " + outdir)
            ancdir = d + "/asr." + model + "/tree1"
            for f in os.listdir(ancdir):
                if f.__contains__(".dat"):            
                    data = get_pp_distro( ancdir + "/" + f )

                    outpath = outdir + "/" + f.split(".")[0] + ".html" 
                    out = ""

                    # writes the script for a Google barplot
                    google_header = write_ppdistro_plot(data)
                    out += get_header(urlpre="../",head=google_header)
                    nodenum = f.split(".")[0]
                    nodenum = re.sub("node", "", nodenum)
                    
                    
                    out += "<h2>Ancestral Node " + nodenum + ", Alignment: " + DIR_nick[d] + ", Model: " + model + "</h2>\n"
                    #out += "<hr>\n"
                    out += "<div>\n"
                    out += "<h3>Maximum Likelihood Sequence:</h3>\n" 
                    out += "<p>" + get_ml_sequence_from_file( ancdir+"/"+f ) + "</p>\n"
                    out += "</div>"
        
                    out += "<hr>"

                    out += "<h3><strong>Support Summary:</strong></h3>\n"
                    out += get_ppdistro_summary(data) 
                    # this div is required for the Google barplot...
                    out += "<div id=\"chart_div\" style=\"width: 100%; height: 200px;\"></div>\n"
                    out += "<hr>"
                    out += "<h3>Support by Site:</h3>\n"
                    out += "<p>Site numbers correspond to the alignment.  Sites with insertions/deletions are not shown.</p>\n"
                    #out += "<div align='center'>\n"
                    #out += get_ppdistro_key()
                    #out += "</div>"
                    out += "<div align='left'>"
                    out += pp_distro_to_html( data )
                    out += "</div>\n"

                    out += get_footer()
                    fout = open(outpath, "w")
                    fout.write(out)


def write_anccomp(ap):
    """Writes the header for the anccomp page, plus calls the method write_anccomp_indi."""
    outpath = HTMLDIR + "/anccomp.html"
    fout = open( outpath, "w")
    fout.write( get_header() )
    #fout.write("<h2>Predicting Functional Loci, Using &Delta;F</h2>\n")
    #fout.write("<p>The following ancestors were compared using the &Delta;F metric (Hanson-Smith and Baker, 2014).</p>")

    fout.write("<h3>Select a Phylogenetic Branch:</h3>\n")
    fout.write("<ul>\n")
    for pair in ap.params["compareanc"]:
        write_anccomp_indi(pair, ap)
        indi_path = pair[0] + "to" + pair[1] + ".html"
        fout.write("<li><p>")
        fout.write("<a href='" + indi_path + "'>")
        fout.write(pair[0] + " to " + pair[1])
        fout.write("</a>")
        fout.write("</p></li>\n")
    fout.write("</ul>\n")
    fout.write( get_footer() )
    fout.close()
    
    
def write_anccomp_indi(pair, ap):
    colwidth = {}
    colwidth[0] = "10%"
    colwidth[1] = "16%"
    colwidth[2] = "10%"
    colwidth[3] = "12%"
    colwidth[4] = "12%"
    colwidth[5] = "12%"
    colwidth[6] = "12%"
    colwidth[7] = "16%"
    
    """Writes on HTML page for each ancestral comparison"""
    outpath = HTMLDIR + "/" + pair[0] + "to" + pair[1] + ".html"    
    fout = open( outpath, "w" )
    plotstring = write_anccomp_plot(pair)    
    fout.write( get_header(head=plotstring) )
    fout.write("<h2>" + pair[0] + " &#10142; " + pair[1] + "</h2>\n")
    fout.write("<p>On the phylogenetic branch(es) leading from " + pair[0] + " to " + pair[1] + ", several amino acid substitutions occurred. These substitutions can be classified into the following categories:</p>\n")
    
    fout.write("<ul>\n")    
    fout.write("<li><strong>Type 1:</strong> the maximum likelihood (ML) amino acid changed between " + pair[0] + " and " + pair[1] + ", and there is poor support for " + pair[0] + "'s ML amino acid in " + pair[1] + ", and vice versa.</li>\n")
    fout.write("<li><strong>Type 2:</strong> the ML amino acid changed, but the ML residue in " + pair[0] + " is supported as an alternate residue in " + pair[1] + ", or vice versa.</li>\n")
    fout.write("<li><strong>Type 3:</strong> the ML state did not change, but either " + pair[0] + " or " + pair[1] + " has uncertainty about this state.</li>\n")
    fout.write("<li><strong>Indel Events:</strong> a residue was inserted or deleted in the sequence.</li>\n")
    fout.write("</ul>\n")
        
    #get_ml_sequence_from_file(path, getindels=False)
    #for msa in ap.params["msa_algorithms"]:
    #        for model in ap.params["raxml_models"]:
    #            asrpath1 = msa + "/asr." + model + "/" + pair[0] + ".dat"
    #            asrpath2 = msa + "/asr." + model + "/" + pair[1] + ".dat"
    #            fout.write("<p>" + msa + ":" + model + ":" + get_ml_sequence_from_file(asrpath1, getindels=True) + "<\p>\n")
    #            fout.write("<p>" + msa + ":" + model + ":" + get_ml_sequence_from_file(asrpath1, getindels=True) + "<\p>\n")
    
    
    fin = open(pair[0] + "to" + pair[1] + "/ancestral_changes.txt", "r")
    for l in fin.xreadlines(): 
        fout.write("<table width=\"100%\">\n")
        tokens = l.split("\t")        

        # Write a row for this msa and model
        fout.write("<tr>\n")
        if False == l.__contains__("Alignment") and l.__len__() > 2:
            msa = tokens[0]
            model = tokens[1]
        
        if l.__contains__("Alignment"):
            #fout.write("<th>&nbsp;</th>")
            for tt in range(0, tokens.__len__()):
                align = "left"
                if tt > 1:
                    align = "center"
                t = tokens[tt]
                fout.write("<th width=\"" + colwidth[tt] + "\" align='" + align + "'>" + t + "</th>")
        else:
            #fout.write("<td><p style='tinytext'><a onclick=\"toggle_visibility('tr" + msa + "." + model + "'); toggle_visibility('hide" + msa + "." + model + "'); toggle_visibility('show" + msa + "." + model + "');\">(X)</a></p></td>")
            for tt in range(0, tokens.__len__()):
                align = "left"
                if tt > 1:
                    align = "center"
                t = tokens[tt]
                if tt == tokens.__len__()-1: 

                    fout.write("<td style='tinytext' width=\"" + colwidth[tt] + "\" align='" + align + "'><p>" + t + "</p></td>")
                else:
                    fout.write("<td width=\"" + colwidth[tt] + "\" align='" + align + "'><p>" + t + "</p></td>")

        if False == l.__contains__("Alignment") and l.__len__() > 2:
            fout.write("<td align='center'>")
            fout.write("<div id='show" + msa + "." + model + "' style='display:block; background:#ffff99;'>")
            fout.write("<p><a onclick=\"toggle_visibility('tr" + msa + "." + model + "'); ")
            fout.write("    toggle_visibility('hide" + msa + "." + model + "'); ")
            fout.write("    toggle_visibility('show" + msa + "." + model + "');\">(+) show details</a></p></div>")
            fout.write("<div id='hide" + msa + "." + model + "' style='display:none; background:#ffff99;'>")
            fout.write("<p><a onclick=\"toggle_visibility('tr" + msa + "." + model + "'); ")
            fout.write("    toggle_visibility('hide" + msa + "." + model + "'); ")
            fout.write("    toggle_visibility('show" + msa + "." + model + "');\">(-) hide details</a></p></div>")
            fout.write("</td>")
        else:
            fout.write("<th>&nbsp;</th>")
        fout.write("</tr>\n")
        fout.write("</table>\n")
        #fout.write("<br>\n")

        if False == l.__contains__("Alignment") and l.__len__() > 2:
            fout.write("<div align='left' style='display:none;' id='tr" + msa + "." + model + "'>\n")
            fout.write("<hr class='thinhr'>")
            
            # Key:
            fout.write("<table width='40%' align='center'><tr><td align='center'><p class='smalltext'>Key:</p></td><td class='redrow' align='center'><p class='smalltext'>Type 1</p></td><td class='orangerow' align='center'><p class='smalltext'>Type 2</p></td><td class='bluerow' align='center'><p class='smalltext'>Type 3</p></td><td class='indelrow' align='center'><p class='smalltext'>Indel Event</p></td></tr></table>")
            
            # HTML table, pre-built by anccomp_tools.py
            fout.write("<hr class='thinhr'>")
            fout.write("<div align='left'>\n")
            fin = open(pair[0] + "to" + pair[1] + "/ancestral_changes." + msa + "." + model + ".html", "r")
            for seql in fin.xreadlines():
                fout.write(seql)
            fin.close()
            fout.write("</div>")
            fout.write("<hr class='thinhr'>")
            fout.write("</div>\n")
    fin.close() 
   
    
    fout.write("\n<hr>\n")
    
    #
    # Score across sites
    #
    #fout.write("<div id=\"chart_div\" style=\"width: 100%; height: 300px;\"></div>")
    fout.write("<h2>&Delta;F Scores</h2>\n")
    fout.write("<p>The &Delta;F metric ranks the shift in entropy between two ancestral amino acid probability distributions, ")
    fout.write("from an ancient ancestor and a more-recently derived ancestor.")
    fout.write("Sites with extreme &Delta;F scores (either positive or negative) should be considered strong hypotheses ")
    fout.write(" for genetic loci that caused changes to the overall protein function. Sites with positive &Delta;F scores indicate that ")
    fout.write("the derived ancestor gained conservation at that site, which may have played a role in constructing or shifting ")
    fout.write(" the protein's function. Sites with negative $Delta;F scores indicate that ")
    fout.write("the derived ancestor lost conservation at the site, which may have played a role in degenerating or relaxing ")
    fout.write(" the protein's function.</p>\n")

    fout.write("<div align='center'>")
    fout.write("<table width='50%'>\n")
    fout.write("<tr>")
    fout.write("<td>")
    fout.write("<h3>&Delta;F Scores, Ranked:</h3>\n")
    fout.write("<p>Download: <a href='../" + pair[0] + "to" + pair[1] + "/Df.details.txt'>spreadsheet</a></p>")
    fout.write("</td>")
    fout.write("<td>")
    fout.write("<h3>Full Details of &Delta;F Analaysis:</strong></h3>")
    fout.write("<p>Download: <a href='../" + pair[0] + "to" + pair[1] + "/summary.txt'>spreadsheet</a></p>")
    fout.write("<td>")
    fout.write("</tr>\n")
    fout.write("</table>\n")
    fout.write("</div>\n")
    
    fout.write("<hr class='thinhr'>")
    
    #fout.write("<h3>&Delta;F Scores Across Sequence Sites:</h3>")
    fout.write("<p align='left'>Download: <a href='../" + pair[0] + "to" + pair[1] + "/Df-by-site.w=1.pdf'>pdf</a> | <a href='../" + pair[0] + "to" + pair[1] + "/Df-by-site.w=1.png'>png</a> | <a href='../" + pair[0] + "to" + pair[1] + "/Df-by-site.w=1.pdf.rscript'>R script</a> | <a href='../" + pair[0] + "to" + pair[1] + "/Df.ranked.txt'>spreadsheet</a></p>\n")
    fout.write("<a href='../" + pair[0] + "to" + pair[1] + "/Df-by-site.w=1.pdf'><img src='../" + pair[0] + "to" + pair[1] + "/Df-by-site.w=1.png'></a>\n")
    fout.write("<br>\n")
    
    fout.write("<p align='left'>Download <a href='../" + pair[0] + "to" + pair[1] + "/Df-histogram.pdf'>pdf</a> | <a href='../" + pair[0] + "to" + pair[1] + "/Df-histogram.png'>png</a> | <a href='../" + pair[0] + "to" + pair[1] + "/Df-histogram.rscript'>R script</a></p>")
    fout.write("<a href='../" + pair[0] + "to" + pair[1] + "/Df-histogram.pdf'><img src='../" + pair[0] + "to" + pair[1] + "/Df-histogram.png'></a>\n")
    fout.write("<br>\n")
    
    #
    #
    #
    if ap.params["do_pdb_analysis"]:
        fout.write("<hr>\n")
        fout.write("<h3>&Delta;F Scores on Structural Homology Models</h3>\n")
        fout.write("<table>\n")
        fout.write("<tr>\n")
        for metric in ["Df", "k", "p"]:
            fout.write("<td>")
            fout.write("<a href=\"../" + pair[0] + "to" + pair[1] + "/pymol_script." + metric + "." + pair[1] + ".pse\">")
            fout.write("<img width=\"200px\" height=\"200px\" src=\"../" + pair[0] + "to" + pair[1] + "/pymol_ray." + metric + "." + pair[1] + ".png\">\n")
            fout.write("</a>")
            fout.write("</td>")
        fout.write("</tr>\n")

        fout.write("<tr>\n")
        for metric in ["Df", "k", "p"]:
            fout.write("<td align='center'><p>" + metric + " scores</p></td>")
        fout.write("</tr>")
        
        fout.write("<tr>\n")
        for metric in ["Df", "k", "p"]:
            fout.write("<td align='center'><p>Download: <a href=\"../" + pair[0] + "to" + pair[1] + "/pymol_script." + metric + "." + pair[1] + ".pse\"><p>pymol</a></p></td>")
        fout.write("</tr>")
        fout.write("</table>\n")
    
    #
    # Print the detailed sites summary (pre-computed in a text file).
    #
    
    #fout.write("<br>\n")
    fin = open(pair[0] + "to" + pair[1] + "/Df.details.txt", "r")
    msa_site_df = {}
    msa_site_rank = {}
    msa_site_pp = {}
    for msa in ap.params["msa_algorithms"]:
        msa_site_df[msa] = {}
        msa_site_rank[msa] = {}
        msa_site_pp[msa] = {}
    msa_site_df["averaged"] = {}
    msa_site_rank["averaged"] = {}
    msa_site_pp["averaged"] = {}
        
    last_msa = "averaged"
    for l in fin.xreadlines():
        if l.startswith("-->"):
            last_msa = "averaged"
            tokens = l.split()
            df = float(tokens[3])
            rank = int(tokens[5])
            site = int(tokens[ tokens.__len__()-2 ])
            msa_site_df[last_msa][site] = df
            msa_site_rank[last_msa][site] = rank
        elif l.__len__() > 2:
            tokens = l.split()
            if tokens[0] in ap.params["msa_algorithms"]:
                last_msa = tokens[0]
                df = float(tokens[5])
                site = int(tokens[2])
                msa_site_df[last_msa][site] = df
                msa_site_pp[last_msa][site] = ""
            elif tokens.__len__() > 3: # skip the seed context lines
                site = int(tokens[2])
                if site in msa_site_pp[last_msa]:
                    msa_site_pp[last_msa][site] += l 
                else:
                    msa_site_pp[last_msa][site] = l
    
    """
    fout.write("<table width='100%'>\n")
    fout.write("<tr>")
    fout.write("<th>Alignment</th><th>Model</th><th>Rel. Probability</th><th>show Df scores</th>\n")
    fout.write("</tr>\n")
    for msa in msa_site_df:
        msanick = msa.split(".")[0]
        model
        if msa != "averaged":
            print msa
            #model = msa.split(".")[1]
            model_data = read_lnllog(msanick)

        fout.write("<tr>")
        if msa != "averaged":
            fout.write("<td>" + msa + "</td><td>" + model + "</td><td>" + "%.3f"%model_data[model][1] + "</td><td onclick=\"toggle_visibility('tr" + msa + "');>show Df scores</td>\n")
        else:
            fout.write("<td>averaged</td><td>averaged</td><td>n/a</td><td onclick=\"toggle_visibility('tr" + msa + "');>show Df scores</td>\n")
        fout.write("</tr>\n")    
        fout.write("<tr id='tr" + msa + "'><td  colspan=\"4\"><div width='100%'>")
        fout.write("<table>")
        fout.write("<tr>\n")
        fout.write("<th>Df Score</th><th>Site #</th><th>Site in Seed</th><th>Probabilities:</th>\n")
        fout.write("</tr>\n")
        for site in msa_site_df[msa]:
            fout.write("<tr>\n")
            fout.write("<td><p class='smalltext'>" + msa_site_df[msa][site] + "</p></td>")
            fout.write("<td><p class='smalltext'>" + site + "</p></td>")
            fout.write("<td><p class='smalltext'>???</p></td>")
            fout.write("<td><p class='smalltext'>" + msa_site_pp[msa][site] + "</p></td>")
            fout.write("</tr>\n")
        fout.write("</table>")
        fout.write("</div></td></tr>\n")
    fout.write("</table>\n")
    """
    
    fin.close()
    fout.write("<hr>\n")
    fout.write("")

    
    """
    fout.write("<hr>\n")
    fout.write("<h3>Prediction Summary For Each Site:</h3>\n")
    
    for msa in ap.params["msa_algorithms"]:
        for model in ap.params["raxml_models"]:
            runid = get_runid(msa, model)
            htmlfrag_path = pair[0] + "to" + pair[1] + "/" + runid + ".html"
            if os.path.exists(htmlfrag_path):
                #fout.write("<a onclick=\"ToggleList(" + runid + ")\">toggle</a>")
                fout.write("<div class=\"divInfo\" id=\"" + runid + "\">\n")
                fout.write("<h3>" + runid + "</h3>")
                fin = open(pair[0] + "to" + pair[1] + "/" + runid + ".html", "r")
                for l in fin.xreadlines():
                    fout.write( l )
                fin.close()

                fout.write("</div>\n")
    """
    
    fout.write( get_footer() )
    fout.close()
    
def write_anccomp_plot(pair):
    #fout = open(HTMLDIR + "/sorttable.js", "w")
    #fout.write( get_sorttable() )
    #fout.close()
    
    out = ""
    #
    # Sortable table stuff
    #
    #out += "<script src=\"sorttable.js\"></script>\n"
    
    #
    # Google Chart Stuff
    #
    out += "<script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>"
    out += "<script type=\"text/javascript\">"
    out += "google.load(\"visualization\", \"1\", {packages:[\"corechart\"]});"
    out += "google.setOnLoadCallback(drawChart);"
    out += "function drawChart() {"
    out += "var data = google.visualization.arrayToDataTable(["
    out += "['Site', 'Score'],"
    fin = open(pair[0] + "to" + pair[1] + "/summary.txt")
    site_score = {}
    for l in fin.xreadlines():
        if l.__len__() > 5 and False == l.startswith("site"):
            tokens = l.split()
            site_score[ int(tokens[0]) ] = float(tokens[3])
            out += "[" + tokens[0] + "," + tokens[3] + "],"
    out += "]);"

    out += "var options = {"
    out += "title: 'Delta-f Score',"
    out += "hAxis: {title: 'Site', titleTextStyle: {color: 'red'}}"
    out += "};"

    out += "var chart = new google.visualization.ColumnChart(document.getElementById('chart_div'));"
    out += "chart.draw(data, options);"
    out += "}"
    
    out += "\nfunction ToggleList(IDS) {\n"
    out += "var CState = document.getElementById(IDS);\n"
    out += "if (CState.style.display != \"none\") { CState.style.display = \"none\"; }\n"
    out += "else { CState.style.display = \"block\"; }\n"
    out += "}\n"
    
    out += "</script>"
    
    #
    # Collapsible div stuff. . .
    #
    out += "<script type=\"text/javascript\">\n"
    out += "   function toggle_visibility(id) {\n"
    out += "       var e = document.getElementById(id);\n"
    out += "       if(e.style.display == 'block')\n"
    out += "          e.style.display = 'none';\n"
    out += "       else\n"
    out += "          e.style.display = 'block';\n"
    out += "   }\n"
    out += "</script>\n"

    return out



def write_css():
    if False == os.path.exists(HTMLDIR):
        os.system("mkdir " + HTMLDIR)
    fout = open(HTMLDIR + "/asrpipeline.css", "w")
    fout.write("body {\n")
    fout.write("background-color: #f0f0f0;\n")
    fout.write("/*background-color: #e5e5e5;*/\n")
    fout.write("width: 75%;\n")
    fout.write("margin: auto;\n")
    fout.write("font-family: sans-serif;\n")
    fout.write("line-height: 140%;\n")
    fout.write("}\n")
    fout.write("a:link {\n")
    fout.write("text-decoration: none;\n")
    fout.write("color: #0066ff;\n")
    fout.write("}\n")
    fout.write("a:visited {\n")
    fout.write("color: #0066ff;\n")
    fout.write("text-decoration: none;\n")
    fout.write("}\n")
    fout.write("a:mouseover {\n")
    fout.write("color: #ff9900;\n")
    fout.write("text-decoration: underline;\n")
    fout.write("}\n")
    fout.write("a:hover {\n")
    fout.write("color: #ff9900;\n")
    fout.write("text-decoration: underline;\n")
    fout.write("}\n")
    fout.write("a:active {\n")
    fout.write("text-decoration: none;\n")
    fout.write("color: #0066ff;\n")
    fout.write("}\n")
    fout.write("\n")
    fout.write("h1{ /* Big titles */\n")
    fout.write("color: #444444;\n")
    fout.write("font-size: 25pt;\n")
    fout.write("font-weight: bold;\n")
    fout.write("font-family: Helvetica, Optima, sans-serif;\n")
    fout.write("}\n")
    fout.write("h2{ /* sub-titles, such as \"Section 1\", \"Section 2\"*/\n")
    fout.write("font-size: 16pt;\n")
    fout.write("color: #444444;\n")
    fout.write("font-family: Helvetica, Optima, sans-serif;\n")
    fout.write("}\n")
    fout.write("h3{ /* navigation text, on the top*/\n")
    fout.write("color: #444444;\n")
    fout.write("font-size: 12pt;\n")
    fout.write("font-family: Helvetica, sans-serif;\n")
    fout.write("}\n")
    fout.write("h4{ /* clickable links, usually underneath clickable icons*/\n")
    fout.write("font-size: 10pt;\n")
    fout.write("padding: 0;\n")
    fout.write("color: #444444;\n")
    fout.write("}\n")
    fout.write("h5{ /* misc sub-titles, smaller than h2 and usually in-line */\n")
    fout.write("font-size: 12pt;\n")
    fout.write("color: #444444;\n")
    fout.write("}\n")
    fout.write("hr{\n")
    fout.write("height: 5px;\n")
    fout.write("color: #444444;\n")
    fout.write("background-color: #444444;\n")
    fout.write("}\n")

    fout.write(".thinhr{\n")
    fout.write("height: 1px;\n")
    fout.write("color: #444444;\n")
    fout.write("background-color: #444444;\n")
    fout.write("}\n")
    
    fout.write(".verythinhr{\n")
    fout.write("height: 0.5px;\n")
    fout.write("color: #444444;\n")
    fout.write("background-color: #444444;\n")
    fout.write("}\n")
    
    fout.write("th { font-size: 12pt; color: #444444;}")
    fout.write("td, p, ul{ /* body text */\n")
    fout.write("font-size: 11pt;\n")
    fout.write("color: #444444;\n")
    fout.write("font-family: Helvetica, serif;\n")
    fout.write("line-height: 160%;\n")
    fout.write("}\n")
    fout.write(".bigp{ /* Nav links */\n")
    fout.write("font-size: 12pt;\n")
    fout.write("font-family: Helvetica;\n")
    fout.write("}\n")
    fout.write(".biggerp{\n")
    fout.write("font-size: 16pt;\n")
    fout.write("font-family: Helvetica;\n")
    fout.write("}\n")
    fout.write(".longp{\n")
    fout.write("font-family:Georgia;\n")
    fout.write("}\n")
    fout.write("\n")
    fout.write("img{\n")
    fout.write("border: 0;\n")
    fout.write("}\n")
    fout.write("\n")
    fout.write("\n")
    fout.write(".aqua{color: #3399ff;}\n")
    fout.write(".yellow{color: #ffcc33;}\n")
    fout.write(".white{color: #ffffff;}\n")
    
    fout.write(".red{color:#ff6666;}\n")
    fout.write(".orange{color:#ff6666;}\n")
    fout.write(".blue{color:#ff6666;}\n")
    
    fout.write(".redrow{background: #ffcccc;}\n")
    
    fout.write(".orangerow{background: #ffffcc;}\n")
    fout.write(".whiterow{background: #ffffff;}\n")
    fout.write(".bluerow{background: #C2E0FF;}\n")
    fout.write(".indelrow{background: #ffcc00;}\n")
    
    fout.write(".redtd{background: #ffcccc; font-size: 7pt;}\n")
    fout.write(".orangetd{background: #ffffcc; font-size: 7pt;}\n")
    fout.write(".whitetd{font-size: 7pt;}\n")
    fout.write(".bluetd{background: #C2E0FF; font-size: 7pt;}\n")
    fout.write(".indeltd{background: #ffcc00; font-size: 7pt;}\n")

    fout.write(".smalltext{font-size: 9pt;}\n")
    fout.write(".tinytext{font-size: 6pt;}\n")
    fout.write(".verytinytext{font-size: 3pt;}\n")
    fout.write("div{word-wrap: break-word;}\n")
    fout.write(".ppheader{background: #C0C0C0;}\n")
    fout.write(".headerrow{background: #C0C0C0; font-weight: bold; font-size: 12pt;}\n")
    fout.write(".ppfull{background: #1E90FF;}\n")
    fout.write(".pp9{background: #00BFFF;}\n")
    fout.write(".pp8{background: #90EE90;}\n")
    fout.write(".pp7{background: #ADFF2F;}\n")
    fout.write(".pp6{background: #FFD700;}\n")
    fout.write(".pp5{background: #FFA07A;}\n")
    fout.write(".pplow{background: #FFB6C1;}\n")
    fout.write(".divInfo { display:block; }\n")
    
    
    #fout.write(".#blanket {\n")
    #fout.write(".background-color:#111;\n")
    #fout.write(".opacity: 0.65;\n")
    #fout.write(".filter:alpha(opacity=65);\n")
    #fout.write("position:absolute;\n")
    #fout.write("z-index: 9001;\n")
    #fout.write("top:0px;\n")
    #fout.write("left:0px;\n")
    #fout.write("width:100%;\n")
    #fout.write("}\n")
    #fout.write("#popUpDiv {\n")
    #fout.write("position:absolute;\n")
    #fout.write("background-color:#eeeeee;\n")
    #fout.write("width:300px;\n")
    #fout.write("height:300px;\n")
    #fout.write("z-index: 9002;}\n")
    
    fout.write("\n")
    fout.close()
