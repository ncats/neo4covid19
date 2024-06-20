# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
# Script: sg.py
#
# Aim: Convert the Cytoscape JSON file exported by SmartGraph (https://smartgraph.ncats.io) to dataframe.
#
# References: 
#
# [1] Zahoránszky-Kőhalmi, G., Sheils, T. & Oprea, T.I. SmartGraph: a network pharmacology investigation platform.
#     J Cheminform 12, 5 (2020). DOI: https://doi.org/10.1186/s13321-020-0409-9
#
# [2] SmartGraph [https://smartgraph.ncats.io/]
# 
# [3] Cytoscape.js [https://js.cytoscape.org/]
#
# [4] P. Shannondoi et al. Cytoscape: A Software Environment for Integrated Models of Biomolecular Interaction Networks. Genome Res. 2003. 13: 2498-2504
#     DOI: 10.1101/gr.1239303
#
# Ref: https://www.programiz.com/python-programming/json
# Ref: https://networkx.github.io/documentation/networkx-1.10/tutorial/tutorial.html#adding-attributes-to-graphs-nodes-and-edges
# Ref: http://json.parser.online.fr/
# Ref: https://www.geeksforgeeks.org/reading-and-writing-json-to-a-file-in-python/
# Ref: https://reqbin.com/code/python/poyzx88o/python-requests-get-example#:~:text=To%20send%20an%20HTTP%20GET,method%20with%20the%20headers%3D%20parameter.
# Ref: https://www.w3schools.com/python/ref_requests_get.asp
#


import json
import sys
import pandas as pd
import requests
from n4c_commons import parse_urls
 

FILE_URLS = '../cfg/urls.cfg'







def sg_analysis (sources, targets):
    urls = parse_urls (FILE_URLS)
    url_sg = urls['url_sg'] #Ensure this is env variable is changed to remove '/{source_uniprot_ids}/{target_uniprot_ids}'

    headers = {'Content-Type': 'application/json', 'accept': 'application/json'}

    cargo = {
        'source_uniprot_ids': sources,
        'target_uniprot_ids': targets,
        'shortest_paths': True,
        'max_length': 3,
        'confidence_cutoff': 0.0,
        'directed': True
    }

    print ()
    print ("[->] SmartGraph analysis started. This may take a while ...")
    print ()

    try:
        api_response = requests.post(url = url_sg, data=json.dumps(cargo), headers = headers)

        result = api_response.json()
    except:
        raise Exception ("[ERROR] Something went wrong when calling SmartGraph endpoint.")

    print ("[*] SmartGraph analysis done.")
    print()
    print()

    return (result)




def list2string (l):
    s = ''
    first = True
    
    for i in range(len(l)):
        if first:
            s = l[i]
            first = False
        else:
            s += ';' + l[i]

    return (s)
        

def process_sg (json_sg_api):

    #print (json_sg_api)

    data = json_sg_api

    edges_json = data['edges']

    nodes_json = data['nodes']


    l_nodes_uuid = []
    l_nodes_uniprot = []
    l_nodes_genesymbol = []
    l_nodes_name = []
    l_nodes_synonyms = []
    l_nodes_node_id = []


    for n in nodes_json:

        #print (n)

    
        node_id = n['node_id']
        uuid = n['uuid']
        uniprot_id = n['uniprot_id']
        name = n['fullname']
        synonyms = n['synonyms']
        synonyms = list2string (synonyms)
        gene_symbols = n['gene_symbols']
        gene_symbols = list2string (gene_symbols)

        l_nodes_node_id.append(node_id)
        l_nodes_uuid.append(uuid)
        l_nodes_uniprot.append(uniprot_id)
        l_nodes_name.append(name)
        l_nodes_synonyms.append(synonyms)
        l_nodes_genesymbol.append(gene_symbols)


    df_nodes = pd.DataFrame ({
            'node_id': l_nodes_node_id,
            'uuid': l_nodes_uuid,
            'gene_symbol': l_nodes_genesymbol,
            'uniprot': l_nodes_uniprot,
            'name': l_nodes_name,
            'synonyms': l_nodes_synonyms
        })


    #print (df_nodes)

    l_edges_uuid = []
    l_edges_sourcedb = []
    l_edges_conf = []
    l_edges_ppiuid = []
    l_edges_mechanism = []
    l_edges_source = []
    l_edges_target = []
    l_edges_modulation = []

    for e in edges_json:
        uuid = e['uuid']
    
        source_node = e['start_node']
        target_node = e['end_node']

        edge_type = e['edge_type']
    

    
        mechanism = e['mechanism_details']
        mechanism = list2string (mechanism)
        source_db = e['sourceDB']
        confidence_score = e['max_confidence_value']
        ppiuid = e['edge_label']
        modulation = e['action_type']

        l_edges_uuid.append(uuid)
        l_edges_sourcedb.append(source_db)
        l_edges_conf.append(confidence_score)
        l_edges_ppiuid.append (ppiuid)
        l_edges_mechanism.append(mechanism)
        l_edges_source.append(source_node)
        l_edges_target.append(target_node)
        l_edges_modulation.append(modulation)    

    df_edges = pd.DataFrame ({

        'uuid': l_edges_uuid,
        'sourcedb': l_edges_sourcedb,
        'source_specific_score': l_edges_conf,
        'ppiuid': l_edges_ppiuid,
        'mechanism': l_edges_mechanism,
        'source': l_edges_source,
        'target': l_edges_target,
        'interaction': l_edges_modulation
    })

    #print (df_edges)


    df = df_edges.merge (df_nodes, left_on = 'source', right_on = 'node_id', how = 'inner')
    df = df.rename (columns = {'gene_symbol': 'source_node'})
    df = df.merge (df_nodes, left_on = 'target', right_on = 'node_id', how = 'inner')
    df = df.rename (columns = {'gene_symbol': 'target_node'})

    #print (df)

    #print ('[*] Processing input SmartGraph subnetwork done.')

    return (df)

#process_sg (FILE_sg)
