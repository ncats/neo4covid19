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


FILE_tdl_source = '../data/input/TDL_UniProt_TCRD6_rev.xlsx'
FILE_uniprot2gene = '../data/input/uniprot2gene.tab'
FILE_gene2uniprot = '../data/input/gene2uniprot.tab'



def uniprot2gene (proteins):

	df_checkpoint = pd.DataFrame({'uniprot':proteins})
	df_checkpoint.to_csv ('../data/output/uniprot_mapping_checkpoint.tsv', sep = '\t', index = False)

	url = 'https://www.uniprot.org/uploadlists/'

	query = ''


	for p in proteins:
		query += p + ' '

	query = query.strip()

	params = {
		'from': 'ACC+ID',
		'to': 'GENENAME',
		'format': 'tab',
		'query': query
	}

	print ('[*] Starting ID mapping via UniProt API. This make take a while ...')

	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')
	req = urllib.request.Request(url, data)
	#print(req)
	with urllib.request.urlopen(req) as f:
		response = f.read()
	#print  (response.decode('utf-8'))

	results = response.decode('utf-8')

	print ('[*] ID mapping retrieved from UniProt API.')


	lines = results.split ('\n')

	tmp = []

	genes = []
	
	uniprots = []

	lines.pop()

	first = True

	for i in range (0, len(lines)):
		line = lines[i]
		line = line.strip()
		#print (line)
		tmp = line.split ('\t')
		
		if first:
			first = False
		else:
			uniprots.append (tmp[0].strip())
			genes.append (tmp[1].strip())

	
	df = pd.DataFrame ({'uniprot': uniprots, 'gene': genes})

	return (df)



def get_tdl_by_gene_name (df, test = False):

	# df: join key: gene names in 'gene column'

	df_tdl = pd.read_excel (FILE_tdl_source, sheet_name = 'Sheet1')
	df_tdl = df_tdl.rename (columns = {
			'TDL_2019': 'tdl',
			'UniProt': 'uniprot'
	})

	df_tdl = df_tdl[['uniprot', 'tdl']].copy()

	proteins = list(set(list(df_tdl['uniprot'])))

	#print (proteins)

	if not test:
		df_map = uniprot2gene(proteins)
		df_map.to_csv (FILE_uniprot2gene, sep = '\t', index = False)
	else:
		df_map = pd.read_csv (FILE_uniprot2gene, sep = '\t')

	#print (df_map.head())

	#print (df_map.shape)

	df_tdl = df_tdl.merge(df_map, left_on = 'uniprot', right_on = 'uniprot', how = 'inner')


	df = df.merge (df_tdl, left_on = 'gene', right_on = 'gene', how = 'left')

	df['tdl'] = df['tdl'].fillna('unk')

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

	df['tdl_priority'] = df.apply (lambda x: prioritize_tdl(x['tdl']), axis = 1)

	df = df.sort_values (['tdl_priority'])

	df = df.groupby(['gene'], as_index = False).aggregate ({
					'uniprot': 'first',
					'tdl': 'first',
					'tdl_priority': 'first'
				})

	#print (df)

	return (df)


def gene2uniprot (proteins):
	
	df_checkpoint = pd.DataFrame({'gene_name':proteins})
	df_checkpoint.to_csv ('../data/output/gene_name_mapping_checkpoint.tsv', sep = '\t', index = False)


	url = 'https://www.uniprot.org/uploadlists/'
	
	query = ''

	print ("OK")
	print (proteins)

	for p in proteins:
		query += p + ' '

	query = query.strip()

	params = {
		'from': 'GENENAME',
		'to': 'ACC',
		'format': 'tab',
		'query': query
	}

	print ('[*] Starting ID mapping via UniProt API. This make take a while ...')

	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')
	req = urllib.request.Request(url, data)
	#print(req)
	with urllib.request.urlopen(req) as f:
		response = f.read()
	#print  (response.decode('utf-8'))

	results = response.decode('utf-8')

	print ('[*] ID mapping retrieved from UniProt API.')


	lines = results.split ('\n')

	tmp = []

	genes = []
	
	uniprots = []

	lines.pop()

	first = True

	for i in range (0, len(lines)):
		line = lines[i]
		line = line.strip()
		#print (line)
		tmp = line.split ('\t')
		
		if first:
			first = False
		else:
			genes.append (tmp[0].strip())
			uniprots.append (tmp[1].strip())

	
	df = pd.DataFrame ({'gene': genes, 'uniprot': uniprots})

	return (df)



def get_uniprot_from_gene_name (proteins, test = False):
	#print (proteins)
	# df: join key: gene names in 'gene column'

	#print (proteins)

	if not test:
		df_map = gene2uniprot(proteins)
		df_map.to_csv (FILE_gene2uniprot, sep = '\t', index = False)
	else:
		df_map = pd.read_csv (FILE_gene2uniprot, sep = '\t')

	#print (df_map.head())

	#print (df_map.shape)
	#print (df_map)

	return (list(set(list(df_map['uniprot']))))



