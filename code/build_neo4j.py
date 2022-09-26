# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advancing Translational Sciences (NCATS/NIH)
#
# Ref: https://groups.google.com/forum/#!topic/cytoscape-discuss/Nmm1zhVJpd4
# Ref: https://nbviewer.jupyter.org/github/idekerlab/py2cytoscape/blob/develop/examples/New_wrapper_api_sample.ipynb
# Ref: https://py2neo.org/v4/database.html
# Ref: https://stackoverflow.com/questions/37530309/attributeerror-graph-object-has-no-attribute-cypher-in-migration-of-data-fr
# Ref: https://apassionatechie.wordpress.com/2017/12/27/create-multiple-pandas-dataframe-columns-from-applying-a-function-with-multiple-returns/
#
#
# Ref: https://docs.python.org/3/library/uuid.html



import sys

import uuid

import dbutils as dbu
from dbutils import *
#open_neo4j_connection
# get_all_target_node_metadata, get_all_ppi_edge_metadata

import pandas as pd

from py2neo import Graph

FILE_nodes_targets = '../data/output/nodes_proteins.tab'
FILE_nodes_drugs = '../data/output/nodes_drugs.tab'
FILE_edges_ppi = '../data/output/edges_ppi.tab'
FILE_edges_dtis = '../data/output/edges_dtis.tab'




protein_node_uuids = {}
compound_node_uuids = {}

uuids = []


def extract_data_sources (df):
	unique_abbr_data_sources = list(set(list(df['abbreviated_data_source'])))
	all_unique_abbr_data_sources = []
	
	for u in unique_abbr_data_sources:
		ds = u.split ('||')
		for d in ds:
			all_unique_abbr_data_sources.append (d)
	all_unique_abbr_data_sources = list (set (all_unique_abbr_data_sources))
	
	return (all_unique_abbr_data_sources)



def insert_host_protein_nodes (conn, df):

	#def create_target_node(name, target_type, tdl):
	
	def generate_uuid (node_name):
		ok = False
		u = None
		while not ok:
			u = str(uuid.uuid1())
	

		
			if u not in uuids:
				protein_node_uuids[node_name] = u
				uuids.append (u)
				ok = True
		return (u)

		
	df['uuid'] = df.apply (lambda x: generate_uuid (x['host_protein']), axis = 1)
	

	all_unique_abbr_data_sources = extract_data_sources (df)
	
	
	df.apply (lambda x: create_host_protein_node (conn, x['host_protein'], x['abbreviated_data_source'], all_unique_abbr_data_sources, x['is_experimental'], x['data_source'], x['acquisition_method'], x['source_specific_score'], x['source_specific_score_type'], x['activation'], x['activation_type'], x['metadata'], x['uniprot'], x['tdl'], x['uuid'], ), axis = 1)
	#print (df)

	dbu.create_host_protein_indices(conn)
	

def insert_pathogen_protein_nodes (conn, df):

	#def create_target_node(name, target_type, tdl):
	
	def generate_uuid (node_name):
		ok = False
		u = None
		while not ok:
			u = str(uuid.uuid1())
	

		
			if u not in uuids:
				protein_node_uuids[node_name] = u
				uuids.append (u)
				ok = True
		return (u)

		
	df['uuid'] = df.apply (lambda x: generate_uuid (x['pathogen_protein']), axis = 1)
	
	all_unique_abbr_data_sources = extract_data_sources (df)
	
	
	df.apply (lambda x: create_pathogen_protein_node (conn, x['pathogen_protein'], x['abbreviated_data_source'], all_unique_abbr_data_sources, x['is_experimental'], x['data_source'], x['acquisition_method'], x['source_specific_score'], x['source_specific_score_type'], x['activation'], x['activation_type'], x['metadata'],  x['uuid'], ), axis = 1)
	#print (df)

	dbu.create_pathogen_protein_indices(conn)	
	
	
def insert_drug_nodes (conn, df):

	#def create_target_node(name, target_type, tdl):
	
	def generate_uuid (node_name):
		ok = False
		u = None
		while not ok:
			u = str(uuid.uuid1())
	

		
			if u not in uuids:
				compound_node_uuids[node_name] = u
				uuids.append (u)
				ok = True
		return (u)

		
	df['uuid'] = df.apply (lambda x: generate_uuid (x['drug_name']), axis = 1)
	
	
	all_unique_abbr_data_sources = extract_data_sources (df)
	
		
	df.apply (lambda x: create_drug_node (conn, x['drug_name'], x['abbreviated_data_source'], all_unique_abbr_data_sources, x['data_source'], x['acquisition_method'], x['is_experimental'],  x['smiles'], x['inchi'], x['inchi_key'], x['ns_inchi_key'], x['CAS_RN'], x['metadata'], x['uuid']), axis = 1)
	#print (df)

	dbu.create_compound_indices(conn)	

def insert_ppi_edges (conn, df):
	#df = pd.read_csv (filename, sep = '\t')
	
	### !!!!! ##### This step should be remove, for testing only
	#df = df[df['start_node'].isin(list(node_uuids.keys()))].copy()
	#df = df[df['end_node'].isin(list(node_uuids.keys()))].copy()
	### !!!!! ####



	def get_node_uuid (node_name):
		return (protein_node_uuids[node_name])

		
	def generate_uuid ():
		ok = False
		u = None
		while not ok:
			u = str(uuid.uuid1())
	

		
			if u not in uuids:
				uuids.append (u)
				ok = True
		return (u)

	all_unique_abbr_data_sources = extract_data_sources (df)

	#def cleanup_ds (ds):
	#	return (ds.replace("'", ""))

	#df = df[['start_node', 'end_node', 'DataSource', 'edgeInfo']].copy()
	#df['t2t_uid'] = df.apply (lambda x: generate_t2t_uid (x['start_node'], x['end_node']), axis = 1)
	
	df['source_node_uuid'] = df.apply (lambda x: get_node_uuid (x['host_protein_a']), axis = 1)
	df['target_node_uuid'] = df.apply (lambda x: get_node_uuid (x['host_protein_b']), axis = 1)
	df['uuid'] = df.apply (lambda x: generate_uuid (), axis = 1)
	
	#df = df.fillna ('unkown')

	#df['data_source_mod'] = df.apply (lambda x: cleanup_ds (x['sourceDataSource']), axis = 1) 
	#df = df.drop (columns = ['DataSource'])
	#df = df.rename (columns = {'data_source_mod': 'DataSource'})
	#print (df['start_node'].values)


	
	df.apply (lambda x: dbu.create_ppi_edge (conn, x['host_protein_a'], x['host_protein_b'], x['interaction'], x['mechanism'], x['source_specific_score'], x['is_experimental'], x['data_source'], x['abbreviated_data_source'], all_unique_abbr_data_sources, x['acquisition_method'], x['source_specific_score_type'], x['directed'], x['edge_label'], x['ref_annotation'], x['ref_direction'], x['ref_score'], x['ref_interaction'], x['metadata'], x['source_node_uuid'], x['target_node_uuid'], x['uuid']), axis = 1)






def insert_hpi_edges (conn, df):
	#df = pd.read_csv (filename, sep = '\t')
	
	### !!!!! ##### This step should be remove, for testing only
	#df = df[df['start_node'].isin(list(node_uuids.keys()))].copy()
	#df = df[df['end_node'].isin(list(node_uuids.keys()))].copy()
	### !!!!! ####



	def get_node_uuid (node_name):
		return (protein_node_uuids[node_name])

		
	def generate_uuid ():
		ok = False
		u = None
		while not ok:
			u = str(uuid.uuid1())
	

		
			if u not in uuids:
				uuids.append (u)
				ok = True
		return (u)

	all_unique_abbr_data_sources = extract_data_sources (df)

	#def cleanup_ds (ds):
	#	return (ds.replace("'", ""))

	#df = df[['start_node', 'end_node', 'DataSource', 'edgeInfo']].copy()
	#df['t2t_uid'] = df.apply (lambda x: generate_t2t_uid (x['start_node'], x['end_node']), axis = 1)
	
	df['source_node_uuid'] = df.apply (lambda x: get_node_uuid (x['pathogen_protein']), axis = 1)
	df['target_node_uuid'] = df.apply (lambda x: get_node_uuid (x['host_protein']), axis = 1)
	df['uuid'] = df.apply (lambda x: generate_uuid (), axis = 1)
	
	#df = df.fillna ('unkown')

	#df['data_source_mod'] = df.apply (lambda x: cleanup_ds (x['sourceDataSource']), axis = 1) 
	#df = df.drop (columns = ['DataSource'])
	#df = df.rename (columns = {'data_source_mod': 'DataSource'})
	#print (df['start_node'].values)



	df.apply (lambda x: dbu.create_hpi_edge (conn, x['pathogen_protein'], x['host_protein'], x['interaction'], x['mechanism'], x['source_specific_score'], x['is_experimental'], x['data_source'], x['abbreviated_data_source'], all_unique_abbr_data_sources, x['acquisition_method'], x['source_specific_score_type'], x['directed'], x['edge_label'], x['metadata'], x['source_node_uuid'], x['target_node_uuid'], x['uuid']), axis = 1)




def insert_dti_edges (conn, df):
	#df = pd.read_csv (filename, sep = '\t')
	
	### !!!!! ##### This step should be remove, for testing only
	#df = df[df['start_node'].isin(list(node_uuids.keys()))].copy()
	#df = df[df['end_node'].isin(list(node_uuids.keys()))].copy()
	### !!!!! ####



	def get_compound_node_uuid (node_name):
		return (compound_node_uuids[node_name])

	def get_protein_node_uuid (node_name):
		return (protein_node_uuids[node_name])
		
	def generate_uuid ():
		ok = False
		u = None
		while not ok:
			u = str(uuid.uuid1())
	

		
			if u not in uuids:
				uuids.append (u)
				ok = True
		return (u)
		
		
	all_unique_abbr_data_sources = extract_data_sources (df)


	#def cleanup_ds (ds):
	#	return (ds.replace("'", ""))

	#df = df[['start_node', 'end_node', 'DataSource', 'edgeInfo']].copy()
	#df['t2t_uid'] = df.apply (lambda x: generate_t2t_uid (x['start_node'], x['end_node']), axis = 1)
	
	df['source_node_uuid'] = df.apply (lambda x: get_compound_node_uuid (x['drug_name']), axis = 1)
	df['target_node_uuid'] = df.apply (lambda x: get_protein_node_uuid (x['host_protein']), axis = 1)
	df['uuid'] = df.apply (lambda x: generate_uuid (), axis = 1)
	
	#df = df.fillna ('unkown')

	#df['data_source_mod'] = df.apply (lambda x: cleanup_ds (x['sourceDataSource']), axis = 1) 
	#df = df.drop (columns = ['DataSource'])
	#df = df.rename (columns = {'data_source_mod': 'DataSource'})
	#print (df['start_node'].values)


	
	df.apply (lambda x: dbu.create_dti_edge (conn, x['drug_name'], x['host_protein'], x['action_type'], x['p_chembl'], x['is_experimental'], x['data_source'], x['abbreviated_data_source'], all_unique_abbr_data_sources, x['acquisition_method'], x['source_specific_score_type'], x['directed'], x['source_specific_score'], x['edge_label'], x['metadata'], x['source_node_uuid'], x['target_node_uuid'], x['uuid']), axis = 1)



	
### WORKFLOW ####


def build_neo4j (dbcfg, df_hprot, df_pprot, df_drug, df_hpi, df_ppi, df_dti):

	#df_targets = pd.read_csv (FILE_nodes_targets, sep = '\t')
	#df_drugs = pd.read_csv (FILE_nodes_drugs, sep = '\t')
	#df_ppi = pd.read_csv (FILE_edges_ppi, sep = '\t')
	#df_dti = pd.read_csv (FILE_edges_dtis, sep = '\t')

	print ('')


	print ('\t[->] Connecting to Neo4j db ...')

	## Connect
	conn = dbu.open_neo4j_connection (dbcfg)

	print ('\t[**] Neo4j connection ESTABLISHED.')


	## Reset DB
	print ('\t[->] Deleting all nodes and edges ...')

	#df = test_query(graph)

	dbu.erase_db (conn)

	print ('\t[**] All nodes and edges DELETED.')

	## Insert protein nodes

	print ('\t[->] Create host protein nodes ...')

	insert_host_protein_nodes (conn, df_hprot)

	
	print ('\t[**] .. done')



	print ('\t[->] Create pathogen protein nodes ...')

	insert_pathogen_protein_nodes (conn, df_pprot)
		
	print ('\t[**] .. done')





	## Insert drug nodes

	print ('\t[->] Create drug nodes ...')

	insert_drug_nodes (conn, df_drug)

	print ('\t[**] .. done')



	## Insert HPI edges

	print ('\t[->] Create HPI edges ...')


	insert_hpi_edges (conn, df_hpi)


	print ('\t[**] .. done')
		


	## Insert DTI edges

	print ('\t[->] Create DTI edges ...')


	insert_dti_edges (conn, df_dti)


	print ('\t[**] .. done')
	
	
	

	## Insert PPI edges

	print ('\t[->] Create PPI edges ...')


	insert_ppi_edges (conn, df_ppi)


	print ('\t[**] .. done')
	


	






