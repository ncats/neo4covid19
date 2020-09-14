# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advanving Translational Sciences (NCATS/NIH)
#
# References
#


## 1. Filter Reactom to keep only edges that are in holistic network.
## 2. Flag edges that are in Reactome
## 3. Get regdir and annotation from Reactome
## 4. Create reverse edges for those that are in Reactome but redir is undefined, preserve Reactome annotation


import pandas as pd
import sys

def is_in_reactome (edge, reactome_edges):
	if edge in reactome_edges:
		return (True)
	
	return (False)


def get_reactome_regdir (in_reactome, edge_label, df_reactome):
	result = ''

	if in_reactome:
		df = df_reactome[df_reactome['edge_label'] == edge_label]
		regdir = list(df['regdir'])
		if len(regdir) > 1:
			print ('[ERROR] Violation of 1:1 mapping between holistic network and Reactome PPI edges. Terminating ...')
			print ('Problematic edges:')
			print (df)
			sys.exit(-3)
		else:
			result = regdir[0]
	else:
		result = 'undefined' 
		
	return (result)


def get_reactome_annotation (in_reactome, edge_label, df_reactome):
	result = ''

	if in_reactome:
		df = df_reactome[df_reactome['edge_label'] == edge_label]
		annotation = list(df['Annotation'])
		if len(annotation) > 1:
			print ('[ERROR] Violation of 1:1 mapping between holistic network and Reactome PPI edges. Terminating ...')
			sys.exit(-3)
		else:
			result = annotation[0]
	else:
		result = '' 
		
	return (result)



def create_reverse_edges_based_on_reactome (df):
	df_rev = df[df['is_in_reactome'] == True]
	df_rev = df_rev[df_rev['reactome_regdir'] == 'undefined']

	df_rev['new_source'] = df_rev['target_node']
	df_rev['new_target'] = df_rev['source_node']

	df_rev = df_rev.drop (columns = ['source_node', 'target_node', 'edge_label']).copy()
	
	df_rev = df_rev.rename (columns = {
		'new_source': 'source_node',
		'new_target': 'target_node'
	})

	def create_edge_label (sn, tn):
	        return (sn + '_' + tn)

	df_rev['edge_label'] = df_rev.apply (lambda x: create_edge_label(x['source_node'], x['target_node']), axis = 1)

	df = df.append (df_rev, ignore_index = True)

	return (df)




def get_reactome_score (in_reactome, edge_label, df_reactome):
	result = 0.0

	if in_reactome:
		df = df_reactome[df_reactome['edge_label'] == edge_label]
		score = list(df['Score'])
		if len(score) > 1:
			print ('[ERROR] Violation of 1:1 mapping between holistic network and Reactome PPI edges. Terminating ...')
			sys.exit(-3)
		else:
			result = score[0]
	else:
		result = 0.0
		
	return (result)





def overlay_reactome_ppi (df_all_unique_ppi, df_reactome):

	df_a = df_reactome.groupby(['source_node', 'target_node'], as_index = False).aggregate({'edge_label': 'count'})
	m = max (list(df_a['edge_label']))
	if m > 1:
		print (df_a.sort_values(['edge_label']))
		print ('[ERROR]')
		sys.exit(-7)


	#df_all_unique_ppi = pd.read_csv ('../holistic.tab', sep = '\t')
	#df_reactome = pd.read_csv ('../reactome.tab', sep = '\t')

	df_reactome = df_reactome[df_reactome['edge_label'].isin(list(df_all_unique_ppi['edge_label']))]
 
 

	reactome_edges = list(df_reactome['edge_label'])


	df_all_unique_ppi['is_in_reactome'] = df_all_unique_ppi.apply (lambda x: is_in_reactome (x['edge_label'], reactome_edges), axis = 1)


	df_all_unique_ppi['reactome_regdir'] = df_all_unique_ppi.apply (lambda x: get_reactome_regdir (x['is_in_reactome'], x['edge_label'], df_reactome), axis = 1)
	df_all_unique_ppi['reactome_mechanism'] = df_all_unique_ppi.apply (lambda x: get_reactome_annotation (x['is_in_reactome'], x['edge_label'], df_reactome), axis = 1)
	df_all_unique_ppi['reactome_score'] = df_all_unique_ppi.apply (lambda x: get_reactome_score (x['is_in_reactome'], x['edge_label'], df_reactome), axis = 1)


	df_all_unique_ppi = create_reverse_edges_based_on_reactome (df_all_unique_ppi)

	df_all_unique_ppi['comment'] = df_all_unique_ppi['comment'].fillna('')

	#df_all_unique_ppi.to_csv ('../out_res.tab', sep = '\t', index = False)

	df_agg = df_all_unique_ppi.groupby(['edge_label']).aggregate({'reactome_regdir': 'count'})



	#print (df_agg.shape[0])
	#print (df_all_unique_ppi.shape[0])


	#if df_agg.shape[0] != df_all_unique_ppi.shape[0]:
	#	df_agg = df_agg.sort_values (['reactome_regdir'])
	#	print (df_agg)
	#	print ('[ERROR] Duplicate edges introduced over Reactome overlay process. Please check code. Terminating ...')
	#	sys.exit (-4)

	## Some edges were present in holistic network with both directions even though they are associated with undefined regulation direction according to Reactome.
	## To remedy this and to tidy-up data frame, a final aggregation and column reordering is done.


	df_all_unique_ppi = df_all_unique_ppi.sort_values (['priority'])

	df_all_unique_ppi = df_all_unique_ppi.groupby (['edge_label'], as_index = False).aggregate ({
		'source_node': 'first',
		'target_node': 'first',
		'inxtype': 'first',
		'interaction': 'first',
		'mechanism': 'first',
		'reactome_mechanism': 'first',
		'reactome_regdir': 'first',
		'reactome_score': 'min',
		'metadata': 'first',
		'comment': 'first',
		'pmids': 'first',
		'priority': 'first',
		'source': 'first',
		'data_origin': 'first',
		'relationship': 'first',
		'source_specific_score': 'first',
		'is_in_preprint': 'first',
		'is_in_phipster': 'first',
		'is_in_string': 'first',
		'is_in_hats': 'first',
		'is_in_reactome': 'first'

	})

	df_all_unique_ppi = df_all_unique_ppi[['source_node', 'target_node', 'inxtype', 'interaction', 'mechanism', 'reactome_mechanism', 'reactome_regdir', 'reactome_score', 'metadata', 'comment', 'pmids', 'priority', 'source', 'data_origin', 'relationship', 'source_specific_score', 'is_in_preprint', 'is_in_phipster', 'is_in_string', 'is_in_hats', 'is_in_reactome']].copy()

	if df_agg.shape[0] != df_all_unique_ppi.shape[0]:
		df_agg = df_agg.sort_values (['reactome_regdir'])
		print (df_agg)
		print ('[ERROR] Duplicate edges introduced over Reactome overlay process. Please check code. Terminating ...')
		sys.exit (-4)

	
	
	#print (df_all_unique_ppi.head())

	return (df_all_unique_ppi)
