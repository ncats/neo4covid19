# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NIH/NCATS)
#
# References
#
# Ref: https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers
# Ref: https://www.w3schools.com/python/python_json.asp
# Ref: https://pypi.org/project/requests/2.7.0/
# Ref: https://stackoverflow.com/questions/37453423/python-requests-view-before-sending
# Ref: https://stackoverflow.com/questions/39954475/post-request-works-in-postman-but-not-in-python

# Code origin: https://string-db.org/cgi/help.pl?subpage=api%23getting-all-the-string-interaction-partners-of-the-protein-set and from Lars J. Jensen, PhD
# 
#
# Cite: https://www.ncbi.nlm.nih.gov/pubmed/30476243
# Szklarczyk et al., Nucleic Acids Res. 2019 Jan 8;47(D1):D607-D613. doi: 10.1093/nar/gky1131

# Ref: https://string-db.org/cgi/help.pl?subpage=api%23getting-all-the-string-interaction-partners-of-the-protein-set
# Ref: https://stackoverflow.com/questions/57458804/attributeerror-module-requests-has-no-attribute-post-is-it-deprecated-and
# Ref: https://stackoverflow.com/questions/32490629/getting-todays-date-in-yyyy-mm-dd-in-python
# Ref: https://www.w3schools.com/python/python_try_except.asp

import requests ## python -m pip install requests
import pandas as pd
from datetime import datetime
import json
from idgtdl import uniprot2gene

date_str = datetime.today().strftime('%Y-%m-%d')

#genes = ["p53", "BRCA1", "cdk2", "Q99835"]



genes_ignore_due_problem = ['ELOC', 'EP300', 'SLC25A5', 'TUBA1A', 'STAT1', 'ELOB', 'RBX1', 'CREBBP', 'SKP1']


def filter_problematic_genes (genes):
	genes_ok = []

	for g in genes:
		if g not in genes_ignore_due_problem:
			genes_ok.append(g)

	return (genes_ok)


def stringify (genes = [], species_ncbi = 9606, limit_of_mapped_genes = 1, max_interactor = 100, score_cutoff = 0, alpha = 0.5, request_id = 'Neo4COVID19', comment = '', priority = 0):

	#print ("[*] genes for stringifying:")
	#print ()
	#print (genes)

	genes = filter_problematic_genes (genes)

	string_api_url = "https://string-db.org/api"
	output_format = "tsv-no-header"
	method = "get_string_ids"



	df = pd.DataFrame()
	df_all = pd.DataFrame ()
	
	sn = []
	tn = []
	scores = []

	##
	## Construct the request
	##

	request_url = "/".join([string_api_url, output_format, method])

	##
	## Set parameters
	##
	params = {

		"identifiers" : "\n".join(genes), # your protein list
		"species" : species_ncbi,	# species NCBI identifier 
		"limit" : limit_of_mapped_genes,	# only one (best) identifier per input protein
		"echo_query" : 1, # see your input identifiers in the output
		"caller_identity" : request_id # your app name

	}

	data_source = 'STRING and stringAPP API ' + date_str + ' | '
	data_source += ' max_interactor = ' + str(max_interactor) + ' score_cutoff = ' + str(score_cutoff)
	data_source += ' species_ncbi = ' + str(species_ncbi)
	
	if species_ncbi == 9606:
		data_source += ' (Homo Sapiens)'

	##
	## Call STRING
	##

	response = requests.post(request_url, data=params)
	#print (response.text)


	##
	## Read and parse the results
	##

	string_tids = []
	tids = {}


	for line in response.text.strip().split("\n"):

		l = line.strip().split("\t")

		tid = l[0]

		string_tid = l[2]
		string_tids.append(string_tid)
		
		tids[string_tid] = tid


	### StringApp API ###
	### 
	###

	string_app_api_url = "https://api.jensenlab.org"
	method = "network"

	
	
	request_url = "/".join([string_app_api_url, method])

	params = {
		"entities" : "\n".join(string_tids), # your protein list
		"additional" : max_interactor,        # maximal number of proteins to add
		"alpha" : alpha, # see your input identifiers in the output
		#"species" : species_ncbi,	# species NCBI identifier 
	
	}


	response = requests.post(request_url, data = params)

	#print(response.text)

	data = json.loads(response.text)

	edges_json = data['edges']

	nodes_json = data['nodes']

	uniprot_ids = []
	string_ids = []


	for n in nodes_json:
		string_id = n['@id']
		
		try:
			uniprot_id = n['canonical']
			string_id = string_id.replace('stringdb:', '')
			if (string_id.startswith('9606.')):
				string_ids.append(string_id)
				uniprot_ids.append(uniprot_id)

		except:
			print ('Non-human protein %s ignored.' % (string_id))

	df_ids = pd.DataFrame ({
			'uniprot_id': uniprot_ids,
			'string_id': string_ids,
	})


	
	edge_sources = []
	edge_targets = []
	scores = []

	for e in edges_json:
		source = e['source']
		edge_sources.append(source)	
		target = e['target']
		edge_targets.append(target)	
		score = e['scores']['stringdb::score']
		scores.append(score)

	df = pd.DataFrame ({
			'source_node': edge_sources,
			'target_node': edge_targets,
			'source_specific_score': scores
	})



	### Filter out edges containing non-human protein

	df = df[df['source_node'].isin(string_ids)].copy()
	df = df[df['target_node'].isin(string_ids)].copy()

	df = df.merge (df_ids, left_on = 'source_node', right_on = 'string_id', how = 'inner')
	df = df.rename (columns = {'uniprot_id': 'source_uniprot_id'})
	df = df.merge (df_ids, left_on = 'target_node', right_on = 'string_id', how = 'inner')
	df = df.rename (columns = {'uniprot_id': 'target_uniprot_id'})
	df = df.drop(columns = ['source_node', 'target_node']).copy()
	df = df.rename (columns = {'source_uniprot_id': 'source_node'})
	df = df.rename (columns = {'target_uniprot_id': 'target_node'})



	### STRING is not directed, so reverse edges are created

	df_rev = df.rename (columns = {
		'source_node': 'new_target_node',
		'target_node': 'new_source_node'
	}).rename (columns = {
		'new_source_node': 'source_node',
		'new_target_node': 'target_node'
	})	
	
	df = df.append (df_rev, ignore_index = True)
	df = df.groupby(['source_node', 'target_node'], as_index = False).aggregate ({'source_specific_score': 'first'})
	
	### Convert UniProt IDs to gene symbols


	
	proteins = []
	s = list(df['source_node'])
	t = list(df['target_node'])
	
	for n in s:
		if n not in proteins:
			proteins.append(n)
	for n in t:
		if n not in proteins:
			proteins.append(n)
	

	#proteins = list(set(proteins))

	df_genes = uniprot2gene (proteins)

	print (df_genes.columns)


	df = df.merge (df_genes, left_on = 'source_node', right_on = 'uniprot', how = 'inner')
	df = df.drop (columns = ['source_node']).copy()	
	df = df.rename (columns = {'gene': 'source_node'})
	
	df = df.merge (df_genes, left_on = 'target_node', right_on = 'uniprot', how = 'inner')
	df = df.drop (columns = ['target_node']).copy()	
	df = df.rename (columns = {'gene': 'target_node'})
	
	df = df.groupby (['source_node', 'target_node'], as_index = False).aggregate ({
								'source_specific_score': 'first'
							})



	df_all = df
	df_all['source'] = 'STRING'
	df_all['data_origin'] = data_source
	df_all['mechanism'] = 'unknown'
	df_all['interaction'] = 'undefined'
	df_all['priority'] = priority
	df_all['comment'] = comment
	df_all['metadata'] = ''
	df_all['pmids'] = ''	

	df_all = df_all[['source_node', 'target_node', 'interaction', 'mechanism', 'source_specific_score', 'source', 'data_origin', 'priority', 'metadata', 'pmids', 'comment']].copy()

	

#	print (df_all)


	return (df_all)

#stringify(genes)
