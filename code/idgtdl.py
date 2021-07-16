# Author:       Gergely Zahoranszky-Kohalmi, PhD
# 
# Email:        gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advanving Translational Sciences (NCATS/NIH)
#
# References
#
# Uniprot ID Mapping code origin: https://www.uniprot.org/help/api_idmapping
#
# Ref: https://www.uniprot.org/help/api_idmapping
# Ref: https://www.uniprot.org/help/api
#
#

import pandas as pd

import urllib.parse
import urllib.request
import sys

from uniprot_map import *

FILE_tdl_source = '../data/input/TDL_UniProt_TCRD6_rev.xlsx'
FILE_GENE_TO_UNIPROT_TDL = '../data/output/idgtdl_gene2uniprot_mapping.tsv'



def annotate_tdl (df, target_col):

	# df: join key: gene names in 'gene column'

	df_tdl = pd.read_excel (FILE_tdl_source, sheet_name = 'Sheet1')
	df_tdl = df_tdl.rename (columns = {
			'TDL_2019': 'tdl',
			'UniProt': 'uniprot'
	})

	df_tdl = df_tdl[['uniprot', 'tdl']].copy()

	proteins = list(set(list(df_tdl['uniprot'])))

	#print (proteins)

	df_map = uniprot2gene(proteins)


	#print (df_map.head())

	#print (df_map.shape)

	df_tdl = df_tdl.merge(df_map, left_on = 'uniprot', right_on = 'uniprot', how = 'inner')


	df_tdl['tdl'] = df_tdl['tdl'].fillna('unk')

	#print (df)	

	def prioritize_tdl (tdl):
		sc = 99
		if tdl == 'Tclin':
			sc = 1
		elif tdl == 'Tchem':
			sc = 2
		elif tdl == 'Tbio':
			sc = 3
		elif tdl == 'Tdark':
			sc = 4
		elif tdl == 'Tvoid':
			sc = 5
		elif tdl == 'unk':
			sc = 6
		else:
			print ('[ERROR]: Invalid tdl detected. Terminating ...')
			print (tdl)
			sys.exit(-1)	

		return (sc)

	df_tdl['tdl_priority'] = df_tdl.apply (lambda x: prioritize_tdl(x['tdl']), axis = 1)

	df_tdl = df_tdl.sort_values (['tdl_priority'])

	df_tdl = df_tdl.groupby(['gene_symbol'], as_index = False).agg ('first')


	df_tdl.to_csv (FILE_GENE_TO_UNIPROT_TDL, sep = '\t', index = False)

	df = df.merge (df_tdl, left_on = target_col, right_on = 'gene_symbol', how = 'left')

	df['tdl'] = df['tdl'].fillna ('unk')
	
	df = df.drop(columns = ['tdl_priority', 'gene_symbol']).copy()

	#print (df)

	return (df)





