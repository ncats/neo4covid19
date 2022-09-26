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
# Ref; https://stackoverflow.com/questions/60260774/pandas-agg-dropping-columns-lambda-function


import pandas as pd
import sys



def extract_metadata (df, md_cols):
	orig_cols = {}
	length = df.shape[0]
	first = True
	
	md = []
	

	
	if md_cols != 'na':
	
		md_cols = md_cols.split (';')
	
		for col in md_cols:
			col = col.strip()
			orig_cols[col] = list (map(str, df[col]))
			if first:
				first = False

		for i in range(length):
		
			combined_md = ''
			
			for col in orig_cols:

			
				md_val = orig_cols [col][i]
				if md_val == 'nan':
					md_val = ''

	

				combined_md += col + ':' + md_val + ';'
			

			md.append (combined_md[:-1])
	
	else:
		md = ['na' for i in range (length)]
		
	return (md)


def merge_metadata (data_type, df, md):
	if data_type == 'host_protein':
		host_proteins = list(df['host_protein'])
		activations = list(df['activation'])
		activation_types = list(df['activation_type'])

		df = pd.DataFrame ({'host_protein':host_proteins, 'activation':activations, 'activation_type':activation_types, 'metadata': md})
	
	elif data_type == 'pathogen_protein':
		pathogen_proteins = list(df['pathogen_protein'])
		df = pd.DataFrame ({'pathogen_protein':pathogen_proteins, 'metadata': md})
	
	elif data_type == 'hpi':
		pathogen_proteins = list(df['pathogen_protein'])
		host_proteins = list(df['host_protein'])
		interactions = list(df['interaction'])
		mechanisms = list(df['mechanism'])
		df = pd.DataFrame ({ 'pathogen_protein':pathogen_proteins, 'host_protein':host_proteins, 'interaction': interactions, 'mechanism': mechanisms, 'metadata': md})
	
	
	elif data_type == 'ppi':
		host_proteins_a = list(df['host_protein_a'])
		host_proteins_b = list(df['host_protein_b'])
		interactions = list(df['interaction'])
		mechanisms = list(df['mechanism'])
		
		df = pd.DataFrame ({ 'host_proteins_a':host_protein_a, 'host_protein_b':host_proteins_b, 'interaction': interactions, 'mechanism': mechanisms, 'metadata': md})
	
	elif data_type == 'dti':
		drugs = list(df['drug_name'])
		host_proteins = list(df['host_protein'])
		action_types = list(df['action_type'])

		df = pd.DataFrame ({'drug_name':drugs, 'host_protein':host_proteins,'action_type': action_types, 'metadata': md})
		
	return (df)
		



def ik2nsik (ik):
	ik = str(ik)
	if ik != 'na':
		return (ik.split('-')[0].strip())
		
	return ('na')

def aggregate_by_first (df, agg_cols):
	df = df.groupby (agg_cols, as_index = False).agg ('first')
	
	return (df)
	
	

def aggregate_by_concat (df, agg_cols, sep = '||'):
	df = df.groupby (agg_cols, as_index = False).agg (lambda x: sep.join([str(i) for i in x]))
	
	return (df)
	
	

def extract_host_proteins_from_ppi (df):
	df_a = df.drop (columns = ['host_protein_b', 'metadata', 'directed', 'interaction', 'mechanism']).copy()
	df_a = df_a.rename (columns = {'host_protein_a': 'host_protein'})
		
	df_b = df.drop (columns = ['host_protein_a', 'metadata', 'directed', 'interaction', 'mechanism']).copy()
	df_b = df_b.rename (columns = {'host_protein_b': 'host_protein'})
	
	df_a['activation'] = 'na'
	df_a['activation_type'] = 'na'
	df_a['metadata'] = 'source:ppi'

	df_b['activation'] = 'na'
	df_b['activation_type'] = 'na'
	df_b['metadata'] = 'source:ppi'
	
	df = df_a.append(df_b, ignore_index = True)
	
	df = df.groupby (['host_protein', 'abbreviated_data_source'], as_index = False).agg ('first')
	
	
	return (df)

def extract_host_proteins_from_hpi (df):
	df = df.drop (columns = ['pathogen_protein', 'interaction', 'mechanism', 'metadata', 'directed']).copy()
	df['activation'] = 'na'
	df['activation_type'] = 'na'
	df['metadata'] = 'source:hpi'
	
	df = df.groupby (['host_protein', 'abbreviated_data_source'], as_index = False).agg ('first')
	
	return (df)


def extract_pathogen_proteins_from_hpi (df):
	df = df.drop (columns = ['host_protein', 'interaction', 'mechanism', 'metadata', 'directed']).copy()
	df['activation'] = 'na'
	df['activation_type'] = 'na'
	df['metadata'] = 'source:hpi'
	
	df = df.groupby (['pathogen_protein', 'abbreviated_data_source'], as_index = False).agg ('first')
	
	return (df)




def extract_host_proteins_from_dti (df):
	df = df.drop (columns = ['drug_name', 'action_type', 'metadata', 'directed']).copy()
	df['activation'] = 'na'
	df['activation_type'] = 'na'
	df['metadata'] = 'source:dti'
	
	df = df.groupby (['host_protein', 'abbreviated_data_source'], as_index = False).agg ('first')
	
	return (df)

def extract_drugs_from_dti (df):
	df = df[['drug_name', 'abbreviated_data_source', 'is_experimental', 'data_source', 'acquisition_method']].copy()
	df['metadata'] = 'source:dti'
	
	df = df.groupby (['drug_name', 'abbreviated_data_source'], as_index = False).agg ('first')
	
	return (df)


def uniqueness_check (df, agg_cols):
	nrow = df.shape[0]
	nrow_agg = df.groupby(agg_cols, as_index = False).agg ('first').shape[0]
	
	if nrow == nrow_agg:
		#print ('Data set %s is unique' % (dataset_abbreviation))

		return (True)
		
	#print ('Data set %s is *not* unique' % (dataset_abbreviation))
			
	return (False)





	
