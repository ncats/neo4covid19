# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advanving Translational Sciences (NCATS/NIH)
#
# References
#
# Ref: https://stackoverflow.com/questions/23418600/python-convert-a-list-of-float-to-string
# Ref: https://stackoverflow.com/questions/6475314/python-for-in-loop-preceded-by-a-variable
# Ref: https://stackoverflow.com/questions/3887381/typeerror-nonetype-object-is-not-iterable-in-python

import pandas as pd

from n4c_commons import*
from drugcentral import *
from ppireactome import *
from idgtdl import *
from build_neo4j import *


# Input

FILE_host_proteins = '../data/input/standardized/STD_host_protein.tsv'
FILE_pathogen_proteins = '../data/input/standardized/STD_pathogen_protein.tsv'
FILE_ppi = '../data/input/standardized/STD_ppi.tsv'
FILE_hpi = '../data/input/standardized/STD_hpi.tsv'
FILE_dti = '../data/input/standardized/STD_dti.tsv'

# Output

FILE_final_hprot = '../data/output/FINAL_host_proteins.tsv'
FILE_final_pprot = '../data/output/FINAL_pathogen_proteins.tsv'
FILE_final_hpi = '../data/output/FINAL_hpi.tsv'
FILE_final_ppi = '../data/output/FINAL_ppi.tsv'
FILE_final_dti = '../data/output/FINAL_dti.tsv'
FILE_final_drug = '../data/output/FINAL_drugs.tsv'

# DB Config
FILE_db_config = '../cfg/db.cfg'

def read_input():
	df_hprot = pd.read_csv (FILE_host_proteins, sep = '\t')
	df_pprot = pd.read_csv (FILE_pathogen_proteins, sep = '\t')
	df_ppi = pd.read_csv (FILE_ppi, sep = '\t')
	df_hpi = pd.read_csv (FILE_hpi, sep = '\t')
	df_dti = pd.read_csv (FILE_dti, sep = '\t')
	
	return (df_hprot, df_pprot, df_ppi, df_hpi, df_dti)
	



def aggregate_input (df_hprot, df_pprot, df_ppi, df_hpi):
	
	
	# Unique host proteins


	agg_cols = ['host_protein']
	df_hprot = aggregate_by_concat (df_hprot, agg_cols, sep = '||')

	
	# Unique pathogen proteins
	agg_cols = ['pathogen_protein']
	df_pprot = aggregate_by_concat (df_pprot, agg_cols, sep = '||')
	
	# Unique PPIs
	
	agg_cols = ['host_protein_a', 'host_protein_b']
	df_ppi = aggregate_by_concat (df_ppi, agg_cols, sep = '||')
	
	
	# Unique HPIs
	
	agg_cols = ['pathogen_protein', 'host_protein']
	df_hpi = aggregate_by_concat (df_hpi, agg_cols, sep = '||')
	
	return (df_hprot, df_pprot, df_ppi, df_hpi)




def aggregate_dtis (df_dc, df_dti):
	df = df_dc.append (df_dti, ignore_index = True)
	
	agg_cols = ['drug_name', 'host_protein']
	
	df = aggregate_by_concat (df, agg_cols, sep = '||')
	
	return (df)
	
	
	
	

def write_output (df_hprot, df_pprot, df_ppi, df_hpi, df_dti, df_drugs):
	df_hprot.to_csv (FILE_final_hprot, sep = '\t', index = False)
	df_pprot.to_csv (FILE_final_pprot, sep = '\t', index = False)
	df_hpi.to_csv (FILE_final_hpi, sep = '\t', index = False)
	df_ppi.to_csv (FILE_final_ppi, sep = '\t', index = False)
	df_dti.to_csv (FILE_final_dti, sep = '\t', index = False)
	df_drugs.to_csv (FILE_final_drug, sep = '\t', index = False)
	

def build ():
	print ('\n\n[-->] Building Neo4COVID19 database ...\n\n')
	
	
	# Process input
	
	print ('[->] Processing input ...')

	df_hprot, df_pprot, df_ppi, df_hpi, df_dti = read_input ()
	
	

	
	# Extract host proteins from SmartGraph & STRING expansion results
	
	df_new_host_proteins = extract_host_proteins_from_ppi (df_ppi)
	agg_cols = ['host_protein',  'activation', 'activation_type', 'metadata', 'abbreviated_data_source']
	

	
	df_new_host_proteins = aggregate_by_first (df_new_host_proteins, agg_cols)
	
	

	df_hprot = df_hprot.append (df_new_host_proteins, ignore_index = True)

	
	
	
	df_hprot, df_pprot, df_ppi, df_hpi = aggregate_input (df_hprot, df_pprot, df_ppi, df_hpi)
	
	df_hpi['edge_label'] = df_hpi.apply (lambda x: '_'.join([x['pathogen_protein'], x['host_protein']]), axis = 1)

	df_ppi['edge_label'] = df_ppi.apply (lambda x: '_'.join([x['host_protein_a'], x['host_protein_b']]), axis = 1)


		
	print ('[**] .. done\n')


	# Identifying DTIs based on host proteins using DrugCentral
	
	print ('[->] Identifying drugs from DrugCentral that act on host proteins ...')

	df_dc = get_dc_dti ()
	df_dc = standardize_dc (df_dc)
	df_dc = df_dc[df_dc['host_protein'].isin (list(df_hprot['host_protein']))].copy()
	df_dti = aggregate_dtis (df_dc, df_dti)
	
	df_dti['edge_label'] = df_dti.apply (lambda x: '_'.join([x['drug_name'], x['host_protein']]), axis = 1)
	
	print ('[**] .. done\n')
		

	# Extract drugs from DTIs
	
	df_drugs = extract_drugs_from_dti (df_dti)
	df_drug_ref = get_drug_ref()
	df_drugs = df_drugs.merge (df_drug_ref, on = 'drug_name', how = 'left')
	
	agg_cols = ['drug_name', 'abbreviated_data_source']
	df_drugs = aggregate_by_first (df_drugs, agg_cols)
	
	agg_cols = ['drug_name']
	df_drugs = aggregate_by_concat (df_drugs, agg_cols)
	#df_drugs = df_drugs.groupby (['drug_name'], as_index = False).agg ('first')


	# Overlay Reactome to PPI
	
	print ('[->] Getting reference PPI data ...')
	df_ppi_ref = get_reactome_ppi ()
	df_ppi_ref = standardize_reactome (df_ppi_ref)
	print ('[**] .. done\n')
	
		
	df_ppi = df_ppi.merge(df_ppi_ref, on = 'edge_label', how= 'left')
	df_ppi = df_ppi.drop(columns = ['prioritized_for_pathway_analysis', 'do_ppi_expansion', 'host_protein_a_y', 'host_protein_b_y']).copy()
	df_ppi = df_ppi.rename (columns = {
		'host_protein_a_x': 'host_protein_a',
		'host_protein_b_x': 'host_protein_b'
	})
	
	# Annotate TDL

	print ('[->] Annotating IDG Target Level Development (TDL) categories of host proteins ...')
	df_hprot = annotate_tdl (df_hprot, 'host_protein')
	print ('[**] .. done\n')

	
	
	df_hprot = df_hprot.fillna('na')
	df_pprot = df_pprot.fillna('na')
	df_ppi = df_ppi.fillna('na')
	df_hpi = df_hpi.fillna ('na')
	df_dti = df_dti.fillna ('na')
	df_drugs = df_drugs.fillna('na')
	
	
	
	# Write output files
	
	
	print ('[->] Writing output files ...')
	write_output (df_hprot, df_pprot, df_ppi, df_hpi, df_dti, df_drugs)
	print ('[**] .. done\n')
	
	
	
	# Build db
	
	print ('[->] Create and populate Neo4COVID19 knowledgebase ...\n')
	
	
	build_neo4j (FILE_db_config, df_hprot, df_pprot, df_drugs, df_hpi, df_ppi, df_dti)
	
	print ('[**] .. done\n')

	
	print ('[Done.]')


	
build ()
