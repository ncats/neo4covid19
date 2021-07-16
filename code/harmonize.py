# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advancing Translational Sciences (NCATS/NIH)
#
# References
#
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.aggregate.html
# Ref: https://stackoverflow.com/questions/27298178/concatenate-strings-from-several-rows-using-pandas-groupby
# Ref: https://stackoverflow.com/questions/32117848/pandas-groupby-concatenate-strings-in-multiple-columns


import pandas as pd
import sys

from standardize import *
from drugcentral import *
from uniprot_map import *
from stringify import stringify
from n4c_commons import *


# Input

FILE_data_registry = '../data/input/data_source_registry.txt'


# Output

STD_hpi = '../data/input/standardized/STD_hpi.tsv'
STD_ppi = '../data/input/standardized/STD_ppi.tsv'
STD_dti = '../data/input/standardized/STD_dti.tsv'
STD_host_protein = '../data/input/standardized/STD_host_protein.tsv'
STD_pathogen_protein = '../data/input/standardized/STD_pathogen_protein.tsv'

FILE_SG_A = '../data/output/sg_proteins_a.tsv'
FILE_SG_B = '../data/output/sg_proteins_b.tsv'
FILE_GENE_TO_UNIPROT = '../data/output/sg_gene2uniprot_mapping.tsv'



def standardize_data (fname_in, schema, params):
	


	if schema == 'dh_krogan':
		df = dh_krogan (fname_in, params)
	elif schema == 'dh_phipster':
		df = dh_phipster (fname_in, params)
	elif schema == 'dh_taiml':
		df = dh_taiml (fname_in, params)
	elif schema == 'dh_nat':
		df = dh_nat (fname_in, params)
	elif schema == 'dh_crispr':
		df = dh_crispr (fname_in, params)
	elif schema == 'dh_hats':
		df = dh_hats (fname_in, params)
	elif schema == 'dh_dc':
		df = dh_dc (params)
	elif schema == 'dh_jm_dti':
		df = dh_jm_dti (fname_in, params)

				
	return (df)
		




def generate_sg_inputs (df_host_proteins):
	df = df_host_proteins
	df_a = df[df['prioritized_for_pathway_analysis'] == True].copy()
	df_b = df[df['prioritized_for_pathway_analysis'] == False].copy()
	
	df_a = df_a[['host_protein', 'abbreviated_data_source']].copy()
	df_a = df_a.groupby(['host_protein'], as_index = False).agg ({'abbreviated_data_source':'first'})
	df_a = df_a[['host_protein']].copy()
	
	df_b = df_b[['host_protein', 'abbreviated_data_source']].copy()
	df_b = df_b.groupby(['host_protein'], as_index = False).agg ({'abbreviated_data_source':'first'})
	df_b = df_b[['host_protein']].copy()
	
	df_a = gene2uniprot(list(df_a['host_protein']))
	df_b = gene2uniprot(list(df_b['host_protein']))
	
	
	df_a = df_a.groupby(['gene_symbol','uniprot'], as_index=False).agg('first')
	df_b = df_b.groupby(['gene_symbol','uniprot'], as_index=False).agg('first')

	
	df_map = df_a.append(df_b, ignore_index = True)
	
	df_map = df_map.groupby(['gene_symbol','uniprot'], as_index=False).agg('first')


	df_a = df_a[['uniprot']].copy()
	df_b = df_b[['uniprot']].copy()
	df_b = df_b[~df_b['uniprot'].isin(list(df_a['uniprot']))].copy()


	df_a.to_csv (FILE_SG_A, sep = '\t', index = False)
	df_b.to_csv (FILE_SG_B, sep = '\t', index = False)


	df_map.to_csv (FILE_GENE_TO_UNIPROT, sep = '\t', index = False)
	



def expand_by_string (df_host_proteins):

	
	df_host_proteins = df_host_proteins[df_host_proteins['do_ppi_expansion'] == True].copy()
	
	df_host_proteins = df_host_proteins.groupby(['host_protein'], as_index = False).agg('first')
	
	unique_human_proteins = list(df_host_proteins['host_protein'])
	
	df = stringify (unique_human_proteins, species_ncbi = 9606, limit_of_mapped_genes = 1, max_interactor = 100, score_cutoff = 0, alpha = 0.5, request_id = 'Neo4COVID19', comment = '', priority = 0)



	df = df.rename (columns = {
					'source_node': 'host_protein_a',
					'target_node': 'host_protein_b',
					'source': 'data_source',
					'data_origin': 'acquisition_method'
				})



	df['interaction'] = 'undefined'
	df['mechanism'] = 'unknown'
	
	
	#print (df.columns)
	#print (df.head())
	
	df = df[['host_protein_a', 'host_protein_b', 'interaction', 'mechanism', 'source_specific_score']].copy()







	#agg_cols = ['host_protein_a', 'host_protein_b']
	#df = aggregate_by_concat (df, agg_cols)

	#df = merge_metadata ('hpi', df, md)
	
	# Generic to all datasets
	df['is_experimental'] = False
	df['data_source'] = 'STRING'
	df['abbreviated_data_source'] = 'ppi_string'
	df['acquisition_method'] = 'STRING API'
	df['source_specific_score_type'] = 'STRING DB score'
	df['prioritized_for_pathway_analysis'] = False
	df['do_ppi_expansion'] = False
	df['directed'] = False
	df['metadata'] = 'na'
	
	
	
	agg_cols = ['host_protein_a', 'host_protein_b', 'interaction', 'mechanism', 'metadata']
	df = aggregate_by_first (df, agg_cols)
	
	#print (df_string.columns)
	#print (df_string.head())

	#print (df_string.head())

#exit()

	print ('[*] STRINGifying done.')

	return (df)
	



def harmonize ():
	
	print ('\n\n[*] Process started ...')
	
	df_hpi = pd.DataFrame ()
	df_ppi = pd.DataFrame ()
	df_dti = pd.DataFrame ()
	df_host_proteins = pd.DataFrame ()
	df_pathogen_proteins = pd.DataFrame ()
	df_drugs = pd.DataFrame ()
	

	
	df_registry = pd.read_csv (FILE_data_registry, sep = '\t')
	
	
	
	#df_drug_ref = get_drug_ref()
	#df_drug_ref.to_csv ('../data/output/tmp_drug_ref.csv', sep = '\t', index = False)


	
	data_sets = list(df_registry['input'])
	harmonization_schemas = list(df_registry['harmonization_schema'])
	data_types = list(df_registry['data_type'])
	is_experimental_data = list(df_registry['is_experimental'])
	data_sources = list(df_registry['data_source'])
	abbreviations = list(df_registry['abbreviated_data_source'])
	original_scores = list(df_registry['original_score'])
	original_score_types = list(df_registry['original_score_type'])
	proteins_for_sg = list(df_registry['prioritized_for_pathway_analysis'])
	proteins_for_string = list(df_registry['do_ppi_expansion'])
	acquisition_methods = list(df_registry['acquisition_method'])
	metadata = list(df_registry['metadata_columns'])
	directed = list(df_registry['is_directed'])

	
	
	for i in range (len(data_sets)):
		
		data_set = data_sets[i]
		
		print ('[*] Processing data set: %s ...' % (data_set))

		
		harmonization_schema = harmonization_schemas[i]
		data_type = data_types[i]

		params = {}
		
		params['data_type'] = data_types[i]
		params['is_experimental'] = is_experimental_data[i]
		params['data_source'] = data_sources[i]
		params['abbreviated_data_source'] = abbreviations[i]
		params['original_score'] = original_scores[i]
		params['original_score_type'] = original_score_types[i]
		params['prioritized_for_pathway_analysis'] = proteins_for_sg[i]
		params['do_ppi_expansion'] = proteins_for_string[i]
		params['acquisition_method'] = acquisition_methods[i]
		params['metadata_columns'] = metadata[i]
		params['directed'] = directed[i]
		

		# Standardize data structure

		df = standardize_data (data_set, harmonization_schema, params)
	
	
			
		# Merge data of same type
	
		if data_type == 'hpi':


			agg_cols = ['host_protein', 'pathogen_protein', 'interaction', 'mechanism', 'metadata']
			if not uniqueness_check (df, agg_cols):
				
				df = aggregate_by_first (df, agg_cols)
				
				
			#agg_cols = ['host_protein', 'pathogen_protein']
				
			#if not uniqueness_check (df, agg_cols):
			#	df = aggregate_by_concat (df, agg_cols)



			df_hpi = df_hpi.append (df)

			df_hprot = extract_host_proteins_from_hpi (df)
			df_host_proteins = df_host_proteins.append (df_hprot, ignore_index = True)

			df_pprot = extract_pathogen_proteins_from_hpi (df)
			df_pathogen_proteins = df_pathogen_proteins.append (df_pprot, ignore_index = True)






		elif data_type == 'dti':
		
			
			agg_cols = ['drug_name', 'host_protein', 'action_type', 'metadata']
			
			if not uniqueness_check (df, agg_cols):
				
				df = aggregate_by_first (df, agg_cols)
	
	
			#agg_cols = ['drug_name', 'host_protein']
			
			#if not uniqueness_check (df, agg_cols):
			#	df = aggregate_by_concat (df, agg_cols)
		
		
			df_dti = df_dti.append (df)
			
			df_hprot = extract_host_proteins_from_dti (df)
			df_host_proteins = df_host_proteins.append (df_hprot, ignore_index = True)



		elif data_type == 'host_protein':
		

			agg_cols = ['host_protein',  'activation', 'activation_type', 'metadata']

			if not uniqueness_check (df, agg_cols):

				df = aggregate_by_first (df, agg_cols)
				
				
			#agg_cols = ['host_protein']

			#if not uniqueness_check (df, agg_cols):

			#	df = aggregate_by_concat (df, agg_cols)
		
		
			df_host_proteins = df_host_proteins.append (df, ignore_index = True)
			
			


		
		else:
			print (data_set)
			print ('[E] Invalid data type defined for data set %s in data source registry. Possible values: [hpi/dti/host_protein]. Terminating ...' % (data_set))
			sys.exit(-1)
	
		print ('[*] .. done')



	return (df_hpi, df_ppi, df_dti, df_host_proteins, df_pathogen_proteins, df_drugs)



def write_standardized_data (df_hpi, df_ppi, df_dti, df_host_proteins, df_pathogen_proteins, df_drugs):
	df_host_proteins.to_csv (STD_host_protein, sep = '\t', index = False)
	df_pathogen_proteins.to_csv (STD_pathogen_protein, sep = '\t', index = False)
	df_hpi.to_csv (STD_hpi, sep = '\t', index = False)
	df_ppi.to_csv (STD_ppi, sep = '\t', index = False)
	df_dti.to_csv (STD_dti, sep = '\t', index = False)
	


df_hpi, df_ppi, df_dti, df_host_proteins, df_pathogen_proteins, df_drugs = harmonize ()


#print (df_hpi, df_ppi, df_dti, df_host_proteins, df_pathogen_proteins, df_drugs)




generate_sg_inputs (df_host_proteins)

df_ppi = expand_by_string (df_host_proteins)


write_standardized_data (df_hpi, df_ppi, df_dti, df_host_proteins, df_pathogen_proteins, df_drugs)

#separate_stringify_inputs ()





print ('\n[Done.]\n')
