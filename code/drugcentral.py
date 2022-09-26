# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
# References
#
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html
# Ref: https://pandas.pydata.org/pandas-docs/stable/user_guide/missing_data.html
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.pivot.html
# Ref: https://stackoverflow.com/questions/43831539/how-to-select-rows-with-nan-in-particular-column
# Ref: https://stackoverflow.com/questions/38309729/count-unique-values-with-pandas-per-groups/38309807
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.astype.html


import pandas as pd
import sys

from n4c_commons import *

url_dc_dti = 'http://unmtid-shinyapps.net/download/drug.target.interaction.tsv.gz'
url_dc_structures = 'http://unmtid-shinyapps.net/download/structures.smiles.tsv'

#local_dc_dti = 'drug.target.interaction.tsv.gz'
#local_dc_structures = 'structures.smiles.tsv'

subst_act = 3

def cleanup_dc (df):
	def simplify_phact (phact):
		phact_map = {}

		phact_map['CHELATING AGENT'] = 'undefined'
		phact_map['NEGATIVE ALLOSTERIC MODULATOR'] = 'down'
		phact_map['UNKNOWN'] = 'undefined'
		phact_map['PHARMACOLOGICAL CHAPERONE'] = 'undefined'
		phact_map['ACTIVATOR'] = 'up'
		phact_map['ANTIBODY BINDING'] = 'undefined'
		phact_map['OXIDATIVE ENZYME'] = 'undefined'
		phact_map['CROSS-LINKING AGENT'] = 'undefined'
		phact_map['SEQUESTERING AGENT'] = 'undefined'
		phact_map['POSITIVE ALLOSTERIC MODULATOR'] = 'up'
		phact_map['ANTISENSE INHIBITOR'] = 'down'
		phact_map['NEGATIVE MODULATOR'] = 'down'
		phact_map['SUBSTRATE'] = 'undefined'
		phact_map['PARTIAL AGONIST'] = 'up'
		phact_map['MEMBRANE PERMEABILIZER'] = 'undefined'
		phact_map['ALLOSTERIC ANTAGONIST'] = 'down'
		phact_map['MINIMUM INHIBITORY CONCENTRATION'] = 'undefined'
		phact_map['INVERSE AGONIST'] = 'down'
		phact_map['ANTAGONIST'] = 'down'
		phact_map['RELEASING AGENT'] = 'undefined'
		phact_map['BLOCKER'] = 'down'
		phact_map['AGONIST'] = 'up'
		phact_map['POSITIVE MODULATOR'] = 'up'
		phact_map['INHIBITOR'] = 'down'
		phact_map['ALKYLATING AGENT'] = 'undefined'
		phact_map['HYDROLYTIC ENZYME'] = 'undefined'
		phact_map['OTHER'] = 'undefined'
		phact_map['OPENER'] = 'undefined'
		phact_map['MODULATOR'] = 'undefined'
		phact_map['GATING INHIBITOR'] = 'down'
		phact_map['ALLOSTERIC MODULATOR'] = 'undefined'
		phact_map['BINDING AGENT'] = 'undefined'
		phact_map['PROTEOLYTIC ENZYME'] = 'undefined'
		phact_map['DNA STRAND BREAK'] = 'undefined'


		simplified_phact = phact_map[phact]

		return (simplified_phact)

	def make_small_case (s):
		return (s.lower())


	def filter_single_protein_dti (ids):
		if '|' not in str(ids):
			return (True)
		
		return (False)

	#-> add simplified pharmacological action
	
	df['ACTION_TYPE'] = df['ACTION_TYPE'].fillna ('UNKNOWN')
	df['simplified_pharm_act'] = df.apply (lambda x: simplify_phact(x['ACTION_TYPE']), axis = 1)
	
	#-> make original pharmacological actions small case
	
	df['SC_ACTION_TYPE'] = df.apply (lambda x: make_small_case(x['ACTION_TYPE']), axis = 1)
	df = df.drop(columns  = ['ACTION_TYPE'])
	df = df.rename (columns = {'SC_ACTION_TYPE': 'ACTION_TYPE'})

	
	colnames = [c.lower() for c in df.columns]

	df.columns = colnames
	
	#print (list(set(list(df_dc_dti['action_type']))))


	#-> keep only DTIs where only a single UniProt ID is involved.

	df['single_protein_dti'] = df.apply (lambda x: filter_single_protein_dti (x['gene']), axis = 1)
	df = df[df['single_protein_dti'] == True].copy()
	df = df.drop (columns = ['single_protein_dti'])

	#-> Keep human targets
	
	df = df[df['organism'] == 'Homo sapiens']

	return (df)




def get_dc_dti ():

	### WORKFLOW ###

	# activity value to substitute missing/inactive (below cutoff) activity values in -logM unit
	

	#is_test_mode = False


	if len(sys.argv) > 1:
		if sys.argv[1] == 'test':
			is_test_mode = True
			print ('[*] Test mode: On.')





	# 1. Getting drugs, DTIs and pharmacological actions from DrugCentral


	#if is_test_mode:
	#	df_dc_dti = pd.read_csv (local_dc_dti, compression = 'gzip', sep = '\t', quotechar = '"')
	#	df_dc_structures = pd.read_csv (local_dc_structures, sep = '\t', quotechar = '"')

	#else:
	df_dc_dti = pd.read_csv (url_dc_dti, compression = 'gzip', sep = '\t', quotechar = '"')



	df_dc_dti = cleanup_dc (df_dc_dti)


	#df_pois = df_dc_dti[df_dc_dti['gene'].isin(pois)].copy()








	df = df_dc_dti

	## Sorting by increasing order of activity, at aggregation the mimimum will be considered to err on conservative site. GZK.
	df = df.sort_values (['act_value'])

	df_aggr = df.groupby(['struct_id', 'gene'], as_index = False).aggregate({
					'act_value': 'first',
					'action_type': 'nunique'
					})

	## Define conflicting pharmacological action types as undefined
	
	df_drugs_with_conflicting_phact = df_aggr[df_aggr['action_type'] > 1]
	drugs_with_conflicting_phact = list(df_drugs_with_conflicting_phact['struct_id'])
	
	df = df.groupby(['struct_id', 'gene'], as_index = False).aggregate({
					'act_value': 'first',
					'action_type': 'first',
					'drug_name': 'first'
					})


	
	df.loc[df['struct_id'].isin(drugs_with_conflicting_phact), 'action_type'] = 'undefined'



	df = df[['struct_id', 'gene', 'act_value', 'action_type', 'drug_name']].copy()
	df = df.rename (columns = {'act_value': 'p_chembl'})
	
	
	df['is_activity_known'] = False
	df.loc[~df['p_chembl'].isnull(), 'is_activity_known'] = True
	
	df.loc[df['is_activity_known'] == False, 'p_chembl'] = subst_act
	
	#df_pois = df_pois.pivot (index = 'struct_id', columns = 'gene', values = 'act_value').reset_index()

	#df_pois = df_pois.rename (columns = {'act_value': 'p_chembl'})

	#df_pois = df_pois.merge (df_dc_structures, left_on = 'struct_id', right_on = 'ID', how = 'inner')

	#df_pois = df_pois.fillna(subst_act)

	#df['source'] = 'DrugCentral'
	#df['priority'] = 1

	return(df)

def annotate_drug_structures (df_drugs, join_type):

	df_dc_structures = pd.read_csv (url_dc_structures, sep = '\t', quotechar = '"')
	df_dc_structures = df_dc_structures.drop(columns = ['INN']).copy()
	df_drugs = df_drugs.merge (df_dc_structures, left_on = 'struct_id', right_on = 'ID', how = join_type)
	
	return (df_drugs)
	
	
def annotate_drug_names (df_drugs, join_type):

	df = pd.read_csv (url_dc_dti, sep = '\t', quotechar = '"')
	df = df.groupby (['STRUCT_ID', 'DRUG_NAME'], as_index = False).aggregate ({'TARGET_NAME': 'first'})
	
	df = df[['STRUCT_ID', 'DRUG_NAME']].copy()
	df_drugs = df_drugs.merge (df, left_on = 'struct_id', right_on = 'STRUCT_ID', how = join_type)
	
	return (df_drugs)	

def identify_drugs_by_name (df_drugs, name_col):
	df_dc_dti = pd.read_csv (url_dc_dti, compression = 'gzip', sep = '\t', quotechar = '"')

	df = df_dc_dti.groupby (['DRUG_NAME', 'STRUCT_ID'], as_index = False).aggregate ({'TARGET_NAME': 'nunique'})
	df = df[['DRUG_NAME', 'STRUCT_ID']].copy()
	
	#df = df_dc_dti.groupby (['STRUCT_ID'], as_index = False).aggregate ({'DRUG_NAME': 'nunique'})


	#print (df.sort_values(['DRUG_NAME']))

	df_drugs = df_drugs.merge (df, left_on = name_col, right_on = 'DRUG_NAME', how = 'left')
	df_drugs.STRUCT_ID = df_drugs.STRUCT_ID.fillna(-1)
	df_drugs = df_drugs.astype ({'STRUCT_ID': 'int32'})
	
	return (df_drugs)


def standardize_dc (df):


	# Provide a comma-seoarated list to define the list of columns that should be extracted as metadata
	md_cols = 'is_activity_known'

	md = extract_metadata (df, md_cols)


	df = df.rename (columns = {
		'gene': 'host_protein'
	})

	df = df[['drug_name', 'host_protein', 'action_type', 'p_chembl', 'struct_id']].copy()
	

	drug_names = list(df['drug_name'])
	host_proteins = list(df['host_protein'])
	action_types = list(df['action_type'])
	p_chembl_values = list(df['p_chembl'])	
	#struct_ids = list(df['struct_id'])	
	
	
	#df = pd.DataFrame ({'drug_name':drug_names, 'host_protein':host_proteins,'action_type': action_types, 'p_chembl': p_chembl_values, 'struct_id':struct_ids, 'metadata': md})
	
	df = pd.DataFrame ({'drug_name':drug_names, 'host_protein':host_proteins,'action_type': action_types, 'p_chembl': p_chembl_values, 'metadata': md})
	
	
	# Specific to DTI

	"""
	df = annotate_drug_structures (df, 'left')
	
	df = df.drop(columns = ['struct_id', 'ID']).copy()
	df = df.rename (columns = {
		'InChI': 'inchi',
		'InChIKey': 'inchi_key',
		'SMILES': 'smiles',
		'gene': 'host_protein',
		'p_chembl': 'source_specific_score'
		
	})
	
	
	


	df['smiles'] = df['smiles'].fillna('na')
	df['inchi'] = df['inchi'].fillna('na')
	df['inchi_key'] = df['inchi_key'].fillna('na')
	df['CAS_RN'] = df['CAS_RN'].fillna('na')
	
	df['ns_inchi_key'] = df.apply (lambda x: ik2nsik (x['inchi_key']), axis = 1)
	"""
	
				
	df['is_experimental'] = True
	df['data_source'] = 'https://drugcentral.org/'
	df['abbreviated_data_source'] = 'dti_dc'
	df['acquisition_method'] = 'web' 
	df['prioritized_for_pathway_analysis'] = False
	df['do_ppi_expansion'] = False
	df['source_specific_score_type'] = 'p_chembl'
	df['directed'] = True

	
	agg_cols = ['drug_name', 'host_protein', 'action_type']
	df = aggregate_by_first (df, agg_cols)
	
	
	return (df)



def get_drug_ref ():
	df = get_dc_dti()
	df = annotate_drug_structures (df, 'left')
	
	
	df = df.drop(columns = ['struct_id', 'ID', 'gene', 'p_chembl', 'action_type', 'is_activity_known']).copy()
	df = df.rename (columns = {
		'InChI': 'inchi',
		'InChIKey': 'inchi_key',
		'SMILES': 'smiles'
		
	})

	df['smiles'] = df['smiles'].fillna('na')
	df['inchi'] = df['inchi'].fillna('na')
	df['inchi_key'] = df['inchi_key'].fillna('na')
	df['CAS_RN'] = df['CAS_RN'].fillna('na')
	df['ns_inchi_key'] = df.apply (lambda x: ik2nsik (x['inchi_key']), axis = 1)
	

	df = df.groupby(['drug_name'], as_index = False).agg ('first')

	return (df)