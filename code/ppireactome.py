# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advanving Translational Sciences (NCATS/NIH)
#
# References
#

import pandas as pd
import sys


FILE_reactome = '../data/input/FIsInGene_020720_with_annotations.txt'



#arrows = list(set(df['Direction']))
#print(arrows)

arrows_left = ['<-', '|-']
arrows_right = ['->', '-|']



arrows_reversed = {}
arrows_reversed['<-'] = '->'
arrows_reversed['|-'] = '-|'

arrows_simplified_right = {}
arrows_simplified_right['|->'] = '->'
arrows_simplified_right['<-|'] = '-|'
arrows_simplified_right['<->'] = '->'
arrows_simplified_right['|-|'] = '-|'
arrows_simplified_right['-'] = '-'



arrows_simplified_left = {}
arrows_simplified_left['|->'] = '|-'
arrows_simplified_left['<-|'] = '<-'
arrows_simplified_left['<->'] = '<-'
arrows_simplified_left['|-|'] = '|-'
arrows_simplified_left['-'] = '-'







# These are only valid after correcting for direction
up_arrows = ['->']
down_arrows = ['-|']





def create_edge_label (gene1, gene2):
	return (gene1 + '_' + gene2)



def determine_new_arrow (direction):
	new_arrow = ''
	if direction in arrows_to_reverse_edge:
		if direction not in list(arrows_reversed.keys()):
			print ('[ERROR] Problem with arrow/direction. Terminating ...')
			sys.exit (-2)

		else:
			new_arrow = arrows_reversed[direction]
	else:
		new_arrow = direction
	
	return(new_arrow)



def get_regdir (direction):
	if direction in up_arrows:
		return ('up')
	elif direction in down_arrows:
		return ('down')
	else:
		return ('undefined')

	return ('undefined')


def split_df_by_right_arrows (df):
	def mark_edge (direction):
		if direction in arrows_right:
			return (True)
		
		return (False)

	df ['split'] = df.apply (lambda x: mark_edge(x['direction']), axis = 1)
	
	df_right = df[df['split'] == True].copy()
	df = df[df['split'] == False].copy()

	df_right = df_right.drop (columns = ['split']).copy()
	df = df.drop (columns = ['split']).copy()



	return (df, df_right)


def split_df_by_left_arrows (df):
	def mark_edge (direction):
		if direction in arrows_left:
			return (True)
		
		return (False)

	df ['split'] = df.apply (lambda x: mark_edge(x['direction']), axis = 1)
	
	df_left = df[df['split'] == True].copy()
	df = df[df['split'] == False].copy()

	df_left = df_left.drop (columns = ['split']).copy()
	df = df.drop (columns = ['split']).copy()


	return (df, df_left)

def swap_nodes (df_mirror):
	df_mirror['source_node'] = df_mirror['Gene2']
	df_mirror['target_node'] = df_mirror['Gene1']

	df_mirror = df_mirror.drop(columns = ['Gene1', 'Gene2', 'edge_label']).copy()

	df_mirror['edge_label'] = df_mirror.apply (lambda x: create_edge_label (x['source_node'], x['target_node']), axis = 1)

	df_mirror = df_mirror[['source_node', 'target_node', 'direction', 'Annotation', 'edge_label', 'Score']].copy()


#	df_mirror['new_direction'] = df_mirror.apply (lambda x: determine_new_arrow (x['direction']), axis = 1)

#	df_mirror = df_mirror.drop(columns = ['direction']).copy()
#	df_mirror = df_mirror.rename (columns = {'new_direction': 'direction'})

	
	return(df_mirror)

def swap_arrows (df):
	def reverse_arrow(direction):
		new_arrow = ''
		if direction not in list(arrows_reversed.keys()):
			print ('[ERROR] Problem with arrow/direction. Terminating ...')
			sys.exit (-2)

		else:
			new_arrow = arrows_reversed[direction]
		
		return(new_arrow)

	
	df['new_direction'] = df.apply (lambda x: reverse_arrow (x['direction']), axis = 1)

	df = df.drop(columns = ['direction']).copy()
	df = df.rename (columns = {'new_direction': 'direction'})

	df = df[['source_node', 'target_node', 'direction', 'Annotation', 'edge_label', 'Score']]


	return(df)


def simplify_arrows_to_right (df):
	def simplify(direction):
		new_arrow = ''
		if direction not in list(arrows_simplified_right.keys()):
			print ('[ERROR] Problem with arrow/direction. Terminating ...')
			sys.exit (-2)

		else:
			new_arrow = arrows_simplified_right[direction]
		
		return(new_arrow)

	
	df['new_direction'] = df.apply (lambda x: simplify (x['direction']), axis = 1)

	df = df.drop(columns = ['direction']).copy()
	df = df.rename (columns = {'new_direction': 'direction'})

	#df = df[['Gene1', 'Gene2', 'direction', 'Annotation', 'edge_label', 'Score']]


	return(df)


def simplify_arrows_to_left (df):
	def simplify(direction):
		new_arrow = ''
		if direction not in list(arrows_simplified_left.keys()):
			print ('[ERROR] Problem with arrow/direction. Terminating ...')
			sys.exit (-2)

		else:
			new_arrow = arrows_simplified_left[direction]
		
		return(new_arrow)

	
	df['new_direction'] = df.apply (lambda x: simplify (x['direction']), axis = 1)

	df = df.drop(columns = ['direction']).copy()
	df = df.rename (columns = {'new_direction': 'direction'})

	#df = df[['Gene1', 'Gene2', 'direction', 'Annotation', 'edge_label', 'Score']]


	return(df)



### Logic of reversing edges or creating duplicate edges:
## 1. Create edge labels by concatenating original source and target nodes (Gene1, Gene2, respectively)
## 2. Check if Reactome input has any duplicate edges based on the edge labels just created.
## 3. If no duplicate edges are present, proceed.
## 4. Separate edges by left-directional, right-directional or bidirectional.
## 5. Reverse left-directional edges and arrows.
## 6. Simplify bidirectional arrows in the original bidirectional copy to right (e.g. '<-|' to '-|', or '<->' to '->'), and to the left in the reverse copy.
##    Of note, '-' will be 'simplified' to remain the same in both directions, i.e. '-' to '-' in right and left simplification.
## 7. Reverse edges (and arrows) in the bidirectional reverse copy.
## 8. Merge edge sets.
## 9. Annotation is preserved in any case, because some of them are complicated (contain information for both direction, etc.)
## 10. Check for duplicate edges, see if there is any conflict.

# Create edge labels





def get_reactome_ppi ():

	df = pd.read_csv (FILE_reactome, sep = '\t')


	unique_inxs = set(list(df['Annotation']))


	df['edge_label'] = df.apply(lambda x: create_edge_label (x['Gene1'], x['Gene2']), axis = 1)

	## Unique edge test passed.


	#print (df.shape)

	df_unique = df.groupby(['edge_label']).aggregate ({'Annotation': 'first'})

	#print (df_unique.shape)

	nr_edges = df.shape[0]
	nr_unique_edges = df_unique.shape[0]

	if nr_edges != nr_unique_edges:
		print ('[ERROR] Reactome input contains duplicate edges. Please check your input. Terminating ...')
		sys.exit (-1)

	df = df.rename (columns = {'Direction': 'direction'})

	#print (df.columns)

	(df, df_right) = split_df_by_right_arrows(df)
	(df, df_left) = split_df_by_left_arrows(df)

	## At this point we have only bidirectional edges. GZK
	df_bidirect = df

	df_bidirect_reverse = df_bidirect.copy()


	df_left = swap_nodes (df_left)
	df_left = swap_arrows(df_left)
	


	df_bidirect = simplify_arrows_to_right(df_bidirect)
	df_bidirect_reverse = simplify_arrows_to_left(df_bidirect_reverse)


	df_bidirect_reverese = swap_nodes (df_bidirect_reverse)
	
	df_right = df_right.rename (columns = {
			'Gene1': 'source_node',
			'Gene2': 'target_node'
		})

	df_bidirect = df_bidirect.rename (columns = {
			'Gene1': 'source_node',
			'Gene2': 'target_node'
		})



	df = df_right.append(df_left, ignore_index = True)
	df = df.append(df_bidirect, ignore_index = True)
	df = df.append(df_bidirect_reverse, ignore_index = True)

	df['regdir'] = df.apply (lambda x: get_regdir (x['direction']), axis = 1)


	df = df.drop (columns = ['Gene1', 'Gene2', 'edge_label']).copy()

	df['edge_label'] = df.apply(lambda x: create_edge_label(x['source_node'], x['target_node']), axis = 1)

	#print (df.head())

	df_aggr = df.groupby (['source_node', 'target_node'], as_index = False).aggregate ({
		'edge_label': 'count'
	})

	df_aggr = df_aggr.sort_values (['edge_label'])

	#print (df_aggr)

	m = max (list(df_aggr['edge_label']))

	if m > 1:
		print ('[ERROR] Reactome processing led to edges that appear more than once. This should not be the case, check input or code. Terminating ...')


	return (df)


def standardize_reactome (df):
	df = df.rename (columns = {
		'source_node': 'host_protein_a',
		'target_node': 'host_protein_b',
		'Annotation': 'ref_annotation',
		'Score': 'ref_score',
		'regdir': 'ref_interaction',
		'direction': 'ref_direction'
	})

	return (df)



