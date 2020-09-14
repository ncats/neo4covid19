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

#import networkx as nx
#from networkx import *
#import networkxgmml
#from networkxgmml import XGMMLWriter

FILE_nodes_targets = '../data/output/nodes_proteins.tab'
FILE_nodes_drugs = '../data/output/nodes_drugs.tab'
FILE_edges_ppi = '../data/output/edges_ppi.tab'
FILE_edges_dtis = '../data/output/edges_dtis.tab'




protein_node_uuids = {}
compound_node_uuids = {}

uuids = []

def insert_protein_nodes (conn, df):

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

		
	df['uuid'] = df.apply (lambda x: generate_uuid (x['gene']), axis = 1)
	df.apply (lambda x: create_target_node (conn, x['gene'], x['is_in_hats'], x['is_in_jdti'], x['is_in_phipster'], x['is_in_preprint'], x['is_in_crispr'], x['is_in_string'], x['is_in_taiml'], x['is_in_natdt'], x['is_in_drugcentral'], x['node_type'], x['target_type'], x['tdl'], x['uniprot'], x['uuid'], x['metadata']), axis = 1)
	#print (df)

	dbu.create_target_indices(conn)
	
	
	
	
def insert_compound_nodes (conn, df):

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
	
	df['CAS_RN'] = df['CAS_RN'].fillna('')
	df['SMILES'] = df['SMILES'].fillna('')
	df['InChI'] = df['InChI'].fillna('')
	df['InChIKey'] = df['InChIKey'].fillna('')
	
	
		
	df.apply (lambda x: create_compound_node (conn, x['struct_id'], x['drug_name'], x['SMILES'], x['InChI'], x['InChIKey'], x['CAS_RN'], x['is_in_drugcentral'], x['is_in_jdti'], x['is_in_hcq'], x['is_in_nhc'], x['is_in_cam'], x['is_in_drugs'], x['node_type'], x['uuid'], x['metadata']), axis = 1)
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

	#def cleanup_ds (ds):
	#	return (ds.replace("'", ""))

	#df = df[['start_node', 'end_node', 'DataSource', 'edgeInfo']].copy()
	#df['t2t_uid'] = df.apply (lambda x: generate_t2t_uid (x['start_node'], x['end_node']), axis = 1)
	
	df['source_node_uuid'] = df.apply (lambda x: get_node_uuid (x['source_node']), axis = 1)
	df['target_node_uuid'] = df.apply (lambda x: get_node_uuid (x['target_node']), axis = 1)
	df['uuid'] = df.apply (lambda x: generate_uuid (), axis = 1)
	
	#df = df.fillna ('unkown')

	#df['data_source_mod'] = df.apply (lambda x: cleanup_ds (x['sourceDataSource']), axis = 1) 
	#df = df.drop (columns = ['DataSource'])
	#df = df.rename (columns = {'data_source_mod': 'DataSource'})
	#print (df['start_node'].values)

	df['comment'] = df['comment'].fillna('')
	df['pmids'] = df['pmids'].fillna('')
	df['metadata'] = df['metadata'].fillna('')
	
	
	df.apply (lambda x: dbu.create_ppi_edge (conn, x['source_node'], x['target_node'], x['interaction'], x['mechanism'], x['reactome_mechanism'], x['reactome_regdir'], x['reactome_score'], x['metadata'], x['comment'], x['pmids'], x['priority'], x['source'], x['data_origin'], x['relationship'], x['source_specific_score'], x['is_in_preprint'], x['is_in_phipster'], x['is_in_string'], x['is_in_hats'], x['is_in_reactome'], x['edge_label'], x['source_node_uuid'], x['target_node_uuid'], x['uuid'], x['inxtype']), axis = 1)


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

	#def cleanup_ds (ds):
	#	return (ds.replace("'", ""))

	#df = df[['start_node', 'end_node', 'DataSource', 'edgeInfo']].copy()
	#df['t2t_uid'] = df.apply (lambda x: generate_t2t_uid (x['start_node'], x['end_node']), axis = 1)
	
	df['source_node_uuid'] = df.apply (lambda x: get_compound_node_uuid (x['drug_name']), axis = 1)
	df['target_node_uuid'] = df.apply (lambda x: get_protein_node_uuid (x['target_node']), axis = 1)
	df['uuid'] = df.apply (lambda x: generate_uuid (), axis = 1)
	
	#df = df.fillna ('unkown')

	#df['data_source_mod'] = df.apply (lambda x: cleanup_ds (x['sourceDataSource']), axis = 1) 
	#df = df.drop (columns = ['DataSource'])
	#df = df.rename (columns = {'data_source_mod': 'DataSource'})
	#print (df['start_node'].values)

	df['comment'] = df['comment'].fillna('')
	df['pmids'] = df['pmids'].fillna('')
	df['metadata'] = df['metadata'].fillna('')
	
	
	df.apply (lambda x: dbu.create_dti_edge (conn, x['action_type'], x['data_origin'], x['drug_name'], x['target_node'], x['is_activity_known'], x['p_chembl'], x['priority'], x['source'], x['source_node'], x['edge_label'], x['relationship'], x['is_in_drugcentral'], x['is_in_jdti'], x['comment'], x['pmids'], x['metadata'], x['source_node_uuid'], x['target_node_uuid'], x['uuid']), axis = 1)



	
### WORKFLOW ####


def build_neo4j (dbcfg, df_targets, df_drugs, df_ppi, df_dti):

	#df_targets = pd.read_csv (FILE_nodes_targets, sep = '\t')
	#df_drugs = pd.read_csv (FILE_nodes_drugs, sep = '\t')
	#df_ppi = pd.read_csv (FILE_edges_ppi, sep = '\t')
	#df_dti = pd.read_csv (FILE_edges_dtis, sep = '\t')

	print ('')


	print ('[->] Connecting to Neo4j db ...')

	## Connect
	conn = dbu.open_neo4j_connection (dbcfg)

	print ('[**] Neo4j connection ESTABLISHED.')


	## Reset DB
	print ('[->] Deleting all nodes and edges ...')

	#df = test_query(graph)

	dbu.erase_db (conn)

	print ('[**] All nodes and edges DELETED.')

	## Insert protein nodes

	print ('[->] Create protein nodes ...')

	insert_protein_nodes (conn, df_targets)

	print ('[**] Protein nodes were CREATED.')


	## Insert PPI edges

	print ('[->] Create PPI edges ...')


	insert_ppi_edges (conn, df_ppi)


	print ('[**] PPI edges were CREATED.')
	
	
	## Insert compound nodes

	print ('[->] Create compound nodes ...')

	insert_compound_nodes (conn, df_drugs)

	print ('[**] Compound nodes were CREATED.')
	
	## Insert DTI edges

	print ('[->] Create DTI edges ...')


	insert_dti_edges (conn, df_dti)


	print ('[**] DTI edges were CREATED.')
	
	print ('[Done.]')





