# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advancing Translational Sciences (NCATS/NIH)
#
# References
#
# Ref: https://stackoverflow.com/questions/23418600/python-convert-a-list-of-float-to-string
# Ref: https://stackoverflow.com/questions/6475314/python-for-in-loop-preceded-by-a-variable


import pandas as pd

from drugcentral import *



FILE_virus_id_map_source_file = '../data/input/Merged.xlsx'
FILE_virus_id_map_source_sheet = 'ID_Mapping'



### Standardization templates

# Implement a harmonization scheme for each dataset that contains either
#	- host proteins
#	- pathogen - host protein interactions
#	- drug - target interactions
#
# according to the mandatory data field listed below (fields marked with asterisk can take 'na' literal string as value is the information is not available):
#
#
#	Host proteins:
#
#		- host_protein  					(string, gene_symbol, key)
#		- abbreviated_data_source			(string)
#		- is_experimental					(boolean)
#		- data_source						(string)
#		- acquisition_method				(string)
#		- prioritized_for_pathway_analysis	(boolean)
#		- do_ppi_expansion					(boolean)
#		- *source_specific_score			(string)
#		- *source_specific_score_type		(string)
#		- *activation						(string)
#		- *activation_type					(string)
#		- *metadata							(string, format:  key1:value1;key2:value2)
#
#
#
#	Host-Pathogen Interactions (HPIs)
#
#		- pathogen_protein  				(string, name)
#		- host_protein  					(string, gene symbol)
#		- interaction  						(string)
#		- mechanism  						(string)
#		- abbreviated_data_source			(string)
#		- is_experimental					(boolean)
#		- data_source						(string)
#		- acquisition_method				(string)
#		- prioritized_for_pathway_analysis	(boolean)
#		- do_ppi_expansion					(boolean)
#		- directed							(boolean)
#		- *source_specific_score			(string)
#		- *source_specific_score_type		(string)
#		- *metadata							(string, format:  key1:value1;key2:value2)
#
#
#
#	Drug-Target Interactions (DTIs)
#
#		
#		- drug_name  						(string, name, key)
#		- host_protein  					(string, gene symbol)
#		- action_type  						(string)
#		- p_chembl  						(string)
#		- abbreviated_data_source			(string)
#		- is_experimental					(boolean)
#		- data_source						(string)
#		- acquisition_method				(string)
#		- prioritized_for_pathway_analysis	(boolean)
#		- do_ppi_expansion					(boolean)
#		- directed							(boolean)
#		- *source_specific_score			(string)
#		- *source_specific_score_type		(string)
#		- *metadata							(string, format:  key1:value1;key2:value2)




# In the future support for additional data types might be possible, such as
#	- protein-protein interactions (between host proteins)
#	- pathogen proteins
#	- drugs.
#
# Please note, that these type of inputs are currently not yet supported, however, PPI, DTI and drug information are pulled in from resources, which logic is implemented in the code, therefore the data format standardization is already taken care of in those cases.
#
#
#	Protein-Protein Interactions (PPIs)
#
#		- host_protein_a  					(string, gene symbol)
#		- host_protein_b  					(string, gene symbol)
#		- interaction  						(string)
#		- mechanism  						(string)
#		- abbreviated_data_source			(string)
#		- is_experimental					(boolean)
#		- data_source						(string)
#		- acquisition_method				(string)
#		- prioritized_for_pathway_analysis	(boolean)
#		- do_ppi_expansion					(boolean)
#		- directed							(boolean)
#		- *source_specific_score			(string)
#		- *source_specific_score_type		(string)
#		- *metadata							(string, format:  key1:value1;key2:value2)
#
#
#
#
#	Pathogen proteins
#
#		- pathogen_protein  				(string, name, key)
#		- abbreviated_data_source			(string)
#		- is_experimental					(boolean)
#		- data_source						(string)
#		- acquisition_method				(string)
#		- prioritized_for_pathway_analysis	(boolean)
#		- do_ppi_expansion					(boolean)
#		- *source_specific_score			(string)
#		- *source_specific_score_type		(string)
#		- *activation						(string)
#		- *activation_type					(string)
#		- *metadata							(string, format:  key1:value1;key2:value2)
#
#	Fields marked with asterisk can take 'na' literal string as value is the information is not available.
#
#
#	Drugs
#
#		- drug_name  						(string, drug name, key)
#		- abbreviated_data_source			(string)
#		- is_experimental					(boolean)
#		- data_source						(string)
#		- acquisition_method				(string)
#		- *smiles							(boolean)
#		- *inchi							(boolean)
#		- *inchi_key						(string)
#		- *ns_inchi_key						(string)
#		- *CAS_RN							(string)
#		- *metadata							(string, format:  key1:value1;key2:value2)



	
	


df_virus_id_map = pd.read_excel (FILE_virus_id_map_source_file, sheet_name = FILE_virus_id_map_source_sheet)	






def standardize_pathogen_protein_names (df):
	def modify_name (name):
		new_name = ''
		name = name.replace ('SARS-CoV2 ', '')
		new_name = 'SARS-CoV2 ' + name
			
		return (new_name)
	
	df['new_name'] = df.apply (lambda x: modify_name (x['pathogen_protein']), axis = 1)
	
	df = df.drop(columns = ['pathogen_protein']).copy()
	df = df.rename (columns = {'new_name': 'pathogen_protein'})

	return (df)



def add_general_fields (df, params):
	# Generic to all datasets
	df['is_experimental'] = params['is_experimental']
	df['data_source'] = params['data_source']
	df['abbreviated_data_source'] = params['abbreviated_data_source']
	df['acquisition_method'] = params['acquisition_method']
	df['prioritized_for_pathway_analysis'] = params['prioritized_for_pathway_analysis']
	df['do_ppi_expansion'] = params['do_ppi_expansion']
	
	if params['original_score'] != 'na':
		df = df.rename (columns = {
			params['original_score']: 'source_specific_score'
		})
	else:
		df['source_specific_score'] = 'na'
	
	df['source_specific_score_type'] = params['original_score_type']
	
	return (df)
	




def dh_krogan (fname_in, params):
	## Experimental host-pathogen interactions (HPIs) by Krogan et al. (https://www.biorxiv.org/content/10.1101/2020.03.22.002386v1)
	
	df = pd.read_csv (fname_in, sep = '\t')

	df = df.rename (columns = {
					'Bait': 'pathogen_protein',
					'PreyGene': 'host_protein'
				})



	df = standardize_pathogen_protein_names (df)


	md = extract_metadata (df, params['metadata_columns'])



	# Interaction specific data
	df['interaction'] = 'undefined'
	df['mechanism'] = 'unknown'
	df = df[['host_protein', 'pathogen_protein', 'interaction', 'mechanism']].copy()

	df = merge_metadata ('hpi', df, md)

	df['directed'] = params['directed']
	
	df = add_general_fields (df, params)


	return (df)

"""
	if params['metadata_columns'] == 'na':
		df['metadata'] = ''
	else:
		df = record_metadata (df, params['metadata_columns'])
"""



def dh_phipster(fname_in, params):
	# Predicted host-pathogen interactions (HPIs) by P-HIPSTer (http://phipster.org/)
	
	
	file_in = fname_in.split (';')[0].strip()
	sn = fname_in.split (';')[1].strip()
	
	df = pd.read_excel (file_in, sheet_name = sn)
	
	
	df = df.merge (df_virus_id_map, left_on = 'Virus Protein', right_on = 'orig_id', how = 'inner')


	df = df.rename (columns = {
					'Virus Protein': 'pathogen_protein',
					'Human Protein': 'host_protein'
				})


	df = standardize_pathogen_protein_names (df)


	md = extract_metadata (df, params['metadata_columns'])

	# Interaction specific data
	df['interaction'] = 'undefined'
	df['mechanism'] = 'unknown'
	df = df[['host_protein', 'pathogen_protein', 'interaction', 'mechanism']].copy()

	df = merge_metadata ('hpi', df, md)
	
	df['directed'] = params['directed']
	
	df = add_general_fields (df, params)

	
	return (df)




def dh_taiml (fname_in, params):
	# Predicted host-pathogen interactions (HPIs) by P-HIPSTer (http://phipster.org/)

	
	file_in = fname_in.split (';')[0].strip()
	sn = fname_in.split (';')[1].strip()
	
	df = pd.read_excel (file_in, sheet_name = sn)
	
	df = df.rename (columns = {
					'Symbol': 'host_protein'
				})


	md = extract_metadata (df, params['metadata_columns'])

	# Specific to host roteins
	df['activation'] = 'na'
	df['activation_type'] = 'na'
	df = df[['host_protein', 'activation', 'activation_type']].copy()	
	
	
	df = merge_metadata ('host_protein', df, md)
	
	df = add_general_fields (df, params)

	
	return (df)
	

def dh_nat (fname_in, params):
	
	
	df = pd.read_csv (fname_in, sep = '\t')

	df = df.rename (columns = {
					'Gene Symbol01': 'host_protein',
					'Quantity_Change': 'activation'
				})



	md = extract_metadata (df, params['metadata_columns'])

	# Specific to host roteins
	df['activation_type'] = 'quantity change'
	df = df[['host_protein', 'activation', 'activation_type']].copy()	
	
	
	df = merge_metadata ('host_protein', df, md)
	
	df = add_general_fields (df, params)

	
	return (df)
	
	
	
def dh_crispr (fname_in, params):
	
	
	
	file_in = fname_in.split (';')[0].strip()
	sn = fname_in.split (';')[1].strip()
	
	df = pd.read_excel (file_in, sheet_name = sn)



	df = df.rename (columns = {
					'gene_symbol': 'host_protein',
				})


	# Specific to host roteins


	md = extract_metadata (df, params['metadata_columns'])

	df = df[['host_protein']].copy()	


	
	df['activation'] = 'na'
	df['activation_type'] = 'na'

	df = df[['host_protein', 'activation', 'activation_type']].copy()	
	
	df = merge_metadata ('host_protein', df, md)

	df = add_general_fields (df, params)

	
	return (df)


	
	
def dh_hats (fname_in, params):
	
	
	df = pd.read_csv (fname_in, sep = '\t')



	df = df.rename (columns = {
					'gene_symbol': 'host_protein',
				})


	# Specific to host roteins


	md = extract_metadata (df, params['metadata_columns'])

	df = df[['host_protein']].copy()	


	
	df['activation'] = 'na'
	df['activation_type'] = 'na'

	df = df[['host_protein', 'activation', 'activation_type']].copy()	
	
	df = merge_metadata ('host_protein', df, md)

	df = add_general_fields (df, params)


	return (df)
	



def dh_jm_dti (fname_in, params):
	file_in = fname_in.split (';')[0].strip()
	sn = fname_in.split (';')[1].strip()
	
	df = pd.read_excel (file_in, sheet_name = sn)
	
	
	
	# Specific to DTI

	md = extract_metadata (df, params['metadata_columns'])

	df = df.rename (columns = {
		'Gene': 'host_protein',
		'Known binders': 'drug_name',
	})
	
	df['action_type'] = 'unknown'
	
	df = df[['drug_name', 'host_protein', 'action_type']].copy()
	
	

	
		
	df = merge_metadata ('dti', df, md)
	
		
	df['is_experimental'] = params['is_experimental']
	df['data_source'] = params['data_source']
	df['abbreviated_data_source'] = params['abbreviated_data_source']
	df['acquisition_method'] = params['acquisition_method']
	df['prioritized_for_pathway_analysis'] = params['prioritized_for_pathway_analysis']
	df['do_ppi_expansion'] = params['do_ppi_expansion']
	
	if params['original_score'] != 'na':
		df = df.rename (columns = {
			params['original_score']: 'source_specific_score'
		})
	else:
		df['source_specific_score'] = 'na'
	

	df['source_specific_score_type'] = params['original_score_type']
	df['directed'] = params['directed']

	#df = df.merge (df_drug_ref, on = 'drug_name', how = 'left')
	
	#print (df.head())
	#print (df.columns)

	return (df)
