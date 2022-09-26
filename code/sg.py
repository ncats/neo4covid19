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


import json
import sys
import pandas as pd

FILE_sg = '../data/input/SG_HATs_dist_3_conf_0.15.json'


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
		

def process_sg (fname):

	with open(fname) as f:
		data = json.load(f)

	edges_json = data['elements']['edges']

	nodes_json = data['elements']['nodes']


	l_nodes_uuid = []
	l_nodes_uniprot = []
	l_nodes_genesymbol = []
	l_nodes_name = []
	l_nodes_synonyms = []


	for n in nodes_json:
		node_record = n['data']
		node_attributes = node_record['node']
	
		uuid = node_attributes['uuid']
		uniprot_id = node_attributes['uniprot_id']
		name = node_attributes['fullname']
		synonyms = node_attributes['synonyms']
		synonyms = list2string (synonyms)
		gene_symbol = node_attributes['genes']

		l_nodes_uuid.append(uuid)
		l_nodes_uniprot.append(uniprot_id)
		l_nodes_name.append(name)
		l_nodes_synonyms.append(synonyms)
		l_nodes_genesymbol.append(gene_symbol)


	df_nodes = pd.DataFrame ({
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
	l_edges_source_uuid = []
	l_edges_target_uuid = []
	l_edges_modulation = []

	for e in edges_json:
		edge_record = e['data']
		uuid = edge_record['id']
	
		source_node = edge_record['source']
		target_node = edge_record['target']





		edge_attributes = edge_record['properties']['properties']
		uuid = node_attributes['uuid']
		edge_type = edge_record['properties']['type']
	

	
		mechanism = edge_record['properties']['properties']['edgeInfo']
		mechanism = list2string (mechanism)
		source_db = edge_record['properties']['properties']['sourceDB']
		confidence_score = edge_record['properties']['max_confidence_value']
		ppiuid = edge_record['properties']['properties']['ppi_uid']
		modulation = edge_record['properties']['edgeType']

		l_edges_uuid.append(uuid)
		l_edges_sourcedb.append(source_db)
		l_edges_conf.append(confidence_score)
		l_edges_ppiuid.append (ppiuid)
		l_edges_mechanism.append(mechanism)
		l_edges_source_uuid.append(source_node)
		l_edges_target_uuid.append(target_node)
		l_edges_modulation.append(modulation)	

	df_edges = pd.DataFrame ({

		'uuid': l_edges_uuid,
		'sourcedb': l_edges_sourcedb,
		'source_specific_score': l_edges_conf,
		'ppiuid': l_edges_ppiuid,
		'mechanism': l_edges_mechanism,
		'source_uuid': l_edges_source_uuid,
		'target_uuid': l_edges_target_uuid,
		'interaction': l_edges_modulation
	})

	#print (df_edges)


	df = df_edges.merge (df_nodes, left_on = 'source_uuid', right_on = 'uuid', how = 'inner')
	df = df.rename (columns = {'gene_symbol': 'source_node'})
	df = df.merge (df_nodes, left_on = 'target_uuid', right_on = 'uuid', how = 'inner')
	df = df.rename (columns = {'gene_symbol': 'target_node'})

	#print (df)

	#print ('[*] Processing input SmartGraph subnetwork done.')

	return (df)

#process_sg (FILE_sg)
