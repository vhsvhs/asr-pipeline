def run_similarity_network_analysis(con):

    #Import Pythoscape
    from pythoscape.auxiliary.re_patterns import RE_PATTERNS
    import pythoscape.main.environments as env
    import pythoscape.interface.local_interface as l_i
    
    #Import plug-ins
    import pythoscape.plugin.input.import_sequences as i_s
    
    #Create interface and environment
    my_interface = l_i.LocalInterface("pythoscape")
    my_pytho = env.PythoscapeEnvironment(my_interface)
    
    #Create plug-ins
    
    #Import sequences
    plugin_1 = i_s.ImportFromFastaFile('sigma_gst.fasta', id_re=RE_PATTERNS['UniprotKB'], id_name='Uniprot')
    
    #Execute plug-ins
    my_pytho.execute_plugin(plugin_1)
    
    #Input edges to Pythoscape
    blastp_exe = get_setting_values(con, "blastp_exe")
    if blastp_exe == None:
        write_log(con, "I can't find your BLAST executable. I'm skipping the similarity network analysis.")
        return
    blastp_exe = blastp_exe[0]
    plugin_2 = a_l_b.AddBLASTEdgesFromLocalBLAST(blastp_exe, include_atts=['-log10(E)'])
    
    #Execute plug-ins
    my_pytho.execute_plugin(plugin_2)
    
    # Skipping plugin 3 (UniProt annotations)
    
    #Create representative network
    cdhit_exe = get_setting_values(con, "cdhit_exe")
    if chdit_exe == None:
        write_log(con, "I canot find your CD-HIT executable. I'm skipping the similarity network analysis.")
        return
    cdhit_exe = cdhit_exe[0]
    plugin_4 = m_c_r.CreateCDHITRepnodes('rep-net', cdhit_exe, c=0.7, n=5)
    
    #Execute plug-ins
    my_pytho.execute_plugin(plugin_4)

    #Create representative network
    plugin_5 = r_s.CalcNodeSize('rep-net','cdhit 0.70','node size')
    
    #Execute plug-ins
    my_pytho.execute_plugin(plugin_5)
    
    #Create representative network
    plugin_6 = a_r_e.AddEdgesToRepnodeNetwork('rep-net','cdhit 0.70','-log10(E)',filt_dir='>',filt_value=0)
    
    #Execute plug-ins
    my_pytho.execute_plugin(plugin_6)

    #Create representative node attributes
    plugin_7 = a_r_a.AddAttributesByIfAny('rep-net','cdhit 0.70','family','repnode family')
    
    #Execute plug-ins
    my_pytho.execute_plugin(plugin_7)
    
    #Output network
    out_plugin_1 = o_x.OutputXGMML('gst-sigma-rep-net-20.xgmml','sigma gst representative node network @ 20',filt_name='rep-net mean',filt_dir='>',filt_value=20)
    
    #Execute plug-ins
    my_net.execute_plugin(out_plugin_1)