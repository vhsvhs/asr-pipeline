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
    h += "<h1>" + ap.params["project_title"] + " Ancestral Library</h1>\n"    
    h += "<hr>\n"
    h += "<p><a href='" + urlpre + "index.html'>Overview</a>"
    h += " | <a href='" + urlpre + "alignments.html'>Alignments</a>" 
    h += " | <a href='" + urlpre + "trees.html'>Trees & Ancestors</a>"
    if os.path.exists("SIMULATION"):
        h += " | <a href='" + urlpre + "errorsimulation.html'>Accuracy Assessment</a>"
    if os.path.exists("SCRIPTS/compareanc_commands.sh"):
        h += " | <a href='" + urlpre + "anccomp.html'>Functional Loci Prediction</a>"    
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
        out += "<td aling='left'>" + d + "</td>"
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
        out += "<h2>" + d + "</h2>"
        out += "<table width=100%>\n"
        out += "<tr>"
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
            out += "<td>" + model + "</td>"
            out += "<td align='right'>" + model_data[model][0].__str__() + "</td>"
            out += "<td align='center'>" + model_data[model][1].__str__() + "</td>"
            out += "<td align='center'>" + "%.3f"%tlength + "</td>"
            out += "<td align='center'><a href='../" + tpath + "'>newick</a></td>"
            out += "<td align='center'><a href='tree." + d + "." + model + ".html'>view</a></td>"
            out += "</tr>\n"
            write_anctree(d, model)
        out += "</table>\n"

    out += get_footer()

    fout = open(HTMLDIR + "/trees.html", "w")
    fout.write( out )
    fout.close()


def write_anctree(d, model):
    rid = get_runid(d, model)
    tpath = get_alrt_treepath(d, model)

    js = "<script type=\"text/javascript\" src=\"../SCRIPTS/HTML_SUPPORT/raphael-min.js\" ></script>\n"
    js += "<script type=\"text/javascript\" src=\"../SCRIPTS/HTML_SUPPORT/jsphylosvg-min.js\"></script>\n"

    js += "<script type=\"text/javascript\">"
    js += "window.onload = function(){\n"

    npath = get_cladogram_path(d, model)
    (xmlpath, xmlstring) = newick_to_xml(d, model)

    js += "       var dataObject = { phyloxml: "
    js += "'" + xmlstring + "' ,"
    #js += "'" + xmlpath + "' ,"
    js += "fileSource:false};\n"

    #js += "YUI().use('oop', 'json-stringify', 'io-base', 'event', 'event-delegate', function(Y){"
    width = 800
    (ntaxa,nsites) = get_phylipstats( get_phylippath(d) )
    height = ntaxa * 16
    #print height
    js += "\nphylocanvas = new Smits.PhyloCanvas(dataObject,\"svgCanvas\","
    js += "800," + height.__str__() + ");\n"
    #js += "});\n" # closes use
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


def write_anccomp():
    """Writes the header for the anccomp page, plus calls the method write_anccomp_indi."""
    outpath = HTMLDIR + "/anccomp.html"
    fout = open( outpath, "w")
    fout.write( get_header() )
    fout.write("<h2>&Delta-F Comparisons;</h2>\n")
    for pair in ap.params["compareanc"]:
        write_anccomp_indi(pair)
        indi_path = pair[0] + "to" + pair[1] + ".html"
        fout.write("<p>")
        fout.write("<a href='" + indi_path + "'>")
        fout.write(pair[0] + " to " + pair[1])
        fout.write("</a>")
        fout.write("</p>\n")
    fout.write( get_footer() )
    fout.close()
    
    
def write_anccomp_indi(pair):
    """Writes on HTML page for each ancestral comparison"""
    outpath = HTMLDIR + "/" + pair[0] + "to" + pair[1] + ".html"    
    fout = open( outpath, "w" )
    plotstring = write_anccomp_plot(pair)
    #print "442:", plotstring
    fout.write( get_header(head=plotstring) )
    fout.write("<h2>Comparison of" + pair[0] + " to " + pair[1] + "</h2>\n")
    
    #
    # Score across sites
    #
    fout.write("<div id=\"chart_div\" style=\"width: 100%; height: 300px;\"></div>")
    fout.write("<img src='../" + pair[0] + "to" + pair[1] + "/hb-by-site.w=1.pdf'>")
    fout.write("<img src='../" + pair[0] + "to" + pair[1] + "/hb-histogram.pdf'>")
    fout.write( get_footer() )
    fout.close()
    
def write_anccomp_plot(pair):
    out = "<script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>"
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
    out += "</script>"
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
    fout.write("\n")
    fout.write(".smalltext{font-size: 9pt;}\n")
    fout.write("div{word-wrap: break-word;}\n")
    fout.write(".ppheader{background: #C0C0C0;}\n")
    fout.write(".ppfull{background: #1E90FF;}\n")
    fout.write(".pp9{background: #00BFFF;}\n")
    fout.write(".pp8{background: #90EE90;}\n")
    fout.write(".pp7{background: #ADFF2F;}\n")
    fout.write(".pp6{background: #FFD700;}\n")
    fout.write(".pp5{background: #FFA07A;}\n")
    fout.write(".pplow{background: #FFB6C1;}\n")
    fout.write("\n")
    fout.close()
