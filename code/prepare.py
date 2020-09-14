# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advanving Translational Sciences (NCATS/NIH)
#
# References
#

# Ref: https://stackoverflow.com/questions/4470523/create-a-branch-in-git-from-another-branch

import pandas as pd
from stringify import stringify
from idgtdl import get_uniprot_from_gene_name #, get_tdl_by_gene_name

#from ppireactome import get_reactome_ppi

from drugcentral import get_dc_dti, annotate_drug_structures, annotate_drug_names, identify_drugs_by_name

#from overlayreactome import overlay_reactome_ppi

from harmonize_drugs import harmonize_dtis

#from build_neo4j import build_neo4j

from crispr import load_crispr

import sys

### Input files

FILE_ppi_krogan = '../data/input/SARS-CoV2.tsv'
FILE_ppi_phipster = '../data/input/PHIPSTER_SARS_COV2.xlsx'
FILE_predicted_targets_taiml = '../data/input/Merged.xlsx'
#FILE_HATS_SmartGraph = '../data/input/HATS_PPI_dist_3_conf_0.15.smartgraph.edges.and.attributes.tab'

# Shared by Jordi Mestres, PhD and Tudor I. Oprea, MD, PhD
FILE_DTI_by_collaborator = '../data/input/COVID19_TargetList_ChemSpace_mod_GZK.xlsx'

FILE_nat_drugtargets = '../data/input/nat_drugtargets.tsv'
#FILE_dbcfg = '../cfg/db.cfg'

FILE_uhps_prestring = '../data/output/unique_host_proteins_prestring.txt'

FILE_chembl_uniprot = '../data/input/chembl_version_27_18052020_uniprot_mapping.tsv'



### Priorities:

PRIOR_DCDRUGS = 0
PRIOR_PREPRINT = 1
PRIOR_NAT = 2
PRIOR_CRISPR = 3
PRIOR_TAIML = 4
PRIOR_PHIPSTER = 5
PRIOR_SG = 6
PRIOR_JDRUGS = 7


def annotate_ds_origin (df, source_name, acquisition_method, priority):
	df['source'] = source_name
	df['data_origin'] = acquisition_method
	df['priority'] = priority

	return(df)



def modify_virus_protein_name (name):
	new_name = ''
	name = name.replace ('SARS-CoV2 ', '')
	new_name = 'SARS-CoV2 ' + name

	return (new_name)

"""

def process_hats_nw (df):

	def extract_mechanism (edge_info):
		mechanism = ''
		tmp = []
		tmp2 = []
		
		tmp = edge_info.split (';')
		for i in range(0, len(tmp)):
			tmp2 = tmp[i].strip().split ('|')
			mechanism += tmp2[0].strip() + '|'

		return (mechanism[:-1])

	def extract_refs (edge_info):
		refs = ''
		tmp = []
		tmp2 = []
		
		tmp = edge_info.split (';')
		for i in range(0, len(tmp)):
			tmp2 = tmp[i].strip().split ('|')
			refs += tmp2[2].strip().replace('pubmed:', '') + '|'

		return (refs[:-1])

	def extract_confidience (edge_info):
		confidence_values = []
		tmp = []
		tmp2 = []
		
		tmp = edge_info.split (';')
		for i in range(0, len(tmp)):
			tmp2 = tmp[i].strip().split ('|')
			confidence_values.append (float(tmp2[3].strip()))

		return (min(confidence_values))
	
	df['mechanism'] = df.apply (lambda x: extract_mechanism(x['edgeInfo']), axis = 1)
	df['pmids'] = df.apply (lambda x: extract_refs(x['edgeInfo']), axis = 1)
	df['source_specific_score'] = df.apply (lambda x: extract_confidience(x['edgeInfo']), axis = 1)

	df = df.rename (columns = {
			'start_node': 'source_node',
			'end_node': 'target_node'
		})

	df['interaction'] = 'undefined' 

	df['metadata'] = ''
	df['comment'] = ''

	df = df[['source_node', 'target_node', 'interaction', 'mechanism', 'source_specific_score', 'metadata', 'pmids', 'comment']].copy()

	return (df)

"""

def update_unique_proteins (df, unique_proteins):
	for p in list(df['source_node']):
		unique_proteins.append (p)

	for p in list(df['target_node']):
		unique_proteins.append (p)

	unique_proteins = list(set(unique_proteins))

	return (unique_proteins)



def create_edge_label (sn, tn):
	return (sn + '_' + tn)


"""
def edges_present_in_dataset (df, df_all, dataset_name, new_column_name):
	df_all = df_all[df_all['source'] == dataset_name]
	reference_edges = list(df_all['edge_label'])
	
	def is_edge_present (edge, reference_edges):
		if edge in reference_edges:
			return (True)
		
		return (False)

	df[new_column_name] = df.apply (lambda x: is_edge_present(x['edge_label'], reference_edges), axis = 1)

	return (df) 	
"""

def record_node_origin (genes, data_source, df_node_registry):
	df = pd.DataFrame({'gene': genes})
	df['data_source'] = data_source
	
	if df_node_registry.empty:
		df_node_registry = df
	else:
		df_node_registry = df_node_registry.append (df)


	return (df_node_registry)


def record_nodes (df_n, data_source, df_node_provenance):

	ds = data_source
	df_n = df_all_ppi[df_all_ppi['source'] == ds]
	unodes = []
	unodes = list(df_n['source_node'])
	unodes.extend(list(df_n['target_node']))
	unodes = list(set(unodes))
	df_node_provenance = record_node_origin (unodes, data_source, df_node_provenance) 

	return (df_node_provenance)



def record_nodes_from_nodes_only_data_source (df, source_column, data_source, df_node_registry):
	genes = list(df[source_column])
	df = pd.DataFrame ({'gene': genes})
	df['data_source'] = data_source

	if df_node_registry.empty:
		df_node_registry = df
	else:
		df_node_registry = df_node_registry.append (df)

	return (df_node_registry)

	

"""

def is_node_present_in_dataset (df_unique_nodes, data_set, df_node_provenance, new_column_name):
	df_x = df_node_provenance[df_node_provenance['data_source'] == data_set]
	reference_nodes = list(df_x['gene'])

	def is_node_in_dataset (n, reference_nodes):
		if n in reference_nodes:
			return (True)

		return (False)


	df_unique_nodes[new_column_name] = df_unique_nodes.apply(lambda x: is_node_in_dataset (x['gene'], reference_nodes), axis = 1)

	return (df_unique_nodes)


def is_edge_in_drug_dataset (df, df_all, dataset, new_column):
	df_all = df_all[df_all['source'] == dataset]

	edges = list(set(list(df_all['edge_label'])))
	df.loc[df['edge_label'].isin(edges), new_column] = True
	df.loc[~df['edge_label'].isin(edges), new_column] = False

	return (df)


def is_node_in_drug_dataset (df, df_all, dataset, new_column):
	df_all = df_all[df_all['source'] == dataset]

	nodes = list(set(list(df_all['struct_id'])))
	df.loc[df['struct_id'].isin(nodes), new_column] = True
	df.loc[~df['struct_id'].isin(nodes), new_column] = False

	return (df)



def is_in_chemical_space_dataset (df, df_all, cs, new_column):
	df_all = df_all[df_all['Chemical space'] == cs]

	nodes = list(set(list(df_all['Known binders'])))
	df.loc[df['drug_name'].isin(nodes), new_column] = True
	df.loc[~df['drug_name'].isin(nodes), new_column] = False

	return (df)

"""


####### WORKFLOW ###


is_test = False


if len(sys.argv) > 1:
	if sys.argv[1] == 'test':
		is_test = True
		print ('[*] Test mode: ON.')


df_node_provenance = pd.DataFrame()

#####
##### 1. Parse input
#####

## PHI-PPI by Krogan et al. (https://www.biorxiv.org/content/10.1101/2020.03.22.002386v1)
# Pathogen - Host and Host-Host Protein Interactions (PHI and PPI)

df_ppi = pd.read_csv (FILE_ppi_krogan, sep = '\t')
df_phipster = pd.read_excel (FILE_ppi_phipster, sheet_name = 'Sheet1')



### Predicted Viral Targets by Meta Path predictions (Tudor I. Oprea, MD PhD)
# Tagets only

df_taiml = pd.read_excel(FILE_predicted_targets_taiml, sheet_name = 'TDL_from_SARS-CoV-2_human')


## Histon Acetyl Transferase (HATs) network by SmartGraph based on Tudor and Tim Wilson's hypothesis
# Human PPIs only, not subject to STRING extension

#df_hats = pd.read_csv (FILE_HATS_SmartGraph, sep ='\t')


# Drug targets from Nature paper (Bojkova et al., DOI: 10.1038/s41586-020-2332-7)
# File created from Supplementary Table 1: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2332-7/MediaObjects/41586_2020_2332_MOESM2_ESM.xlsx
#
# Targets only

df_nat = pd.read_csv(FILE_nat_drugtargets, sep = '\t')
df_crispr = load_crispr ()


### DrugCentral
## Filter DTIs on proteins present in df_human_proteins (unique human proteins)
df_drugs = get_dc_dti()

### DTIs from collaborator
df_drugs_add_1 = pd.read_excel (FILE_DTI_by_collaborator, sheet_name = 'NHC')
df_drugs_add_2 = pd.read_excel (FILE_DTI_by_collaborator, sheet_name = 'HCQ')
df_drugs_add_3 = pd.read_excel (FILE_DTI_by_collaborator, sheet_name = 'CAM')

df_drugs_add = df_drugs_add_1.append(df_drugs_add_2, ignore_index = True)
df_drugs_add = df_drugs_add.append(df_drugs_add_3, ignore_index = True)


#####
##### 2. Consolidate IDs
#####


## Viral Protein ID Mapping between PHIPSTER and Krogan data

df_map = pd.read_excel (FILE_predicted_targets_taiml, sheet_name = 'ID_Mapping')


#print (df_phipster.shape)

df_phipster = df_phipster.merge (df_map, left_on = 'Virus Protein', right_on = 'orig_id', how = 'inner')

#print (df_phipster.shape)
# 1:1 mapping check passed. GZK



#####
##### 3. Harmonizing data structure
#####



df_ppi = df_ppi.rename (columns = {
				'Bait': 'virus_protein',
				'PreyGene': 'human_protein'
			})

df_ppi = df_ppi[['virus_protein', 'human_protein']].copy()

df_ppi['interaction'] = 'undefined'
df_ppi['mechanism'] = 'unknown'
df_ppi['source_specific_score'] = 0
df_ppi['metadata'] = ''
df_ppi['pmids'] = ''
df_ppi['comment'] = ''







df_phipster = df_phipster.rename (columns = {
				'Virus Protein': 'virus_protein',
				'Human Protein': 'human_protein'
			})

df_phipster = df_phipster[['virus_protein', 'human_protein']].copy()

df_ppi['interaction'] = 'undefined'
df_ppi['mechanism'] = 'unknown'
df_ppi['source_specific_score'] = 0
df_ppi['metadata'] = ''
df_ppi['pmids'] = ''
df_ppi['comment'] = ''



df_taiml = df_taiml.rename (columns = {
				'Symbol': 'human_protein'
			})

df_taiml = df_taiml[['human_protein']].copy()

df_nat = df_nat.rename(columns = {'Gene Symbol01': 'human_protein'})
df_nat = df_nat[['human_protein']].copy()

df_crispr = df_crispr.rename (columns = {'gene_symbol': 'human_protein'})

#df_hats = process_hats_nw (df_hats)

df_ppi['interaction'] = 'undefined'
df_ppi['metadata'] = ''
df_ppi['comment'] = ''


df_drugs = df_drugs.rename (columns = {'gene': 'human_protein'})  


# This database has targets implicated
df_drugs_add = df_drugs_add.rename (columns = {'Gene': 'human_protein'})  
df_drugs_add_nodes = df_drugs_add[['human_protein']].copy()


#####
##### 4. Data provenance
#####

# Only Drugs
df_drugs = annotate_ds_origin (df_drugs, source_name = 'DrugCentral', acquisition_method = 'Web download', priority = PRIOR_DCDRUGS)

# Only human target nodes, no PPI edges
df_taiml = annotate_ds_origin (df_taiml, source_name = 'Meta Path AIML', acquisition_method = 'Data share on 2020/03/23', priority = PRIOR_TAIML)
df_drugs_add_nodes = annotate_ds_origin (df_drugs_add_nodes, source_name = 'DTI - Predicted', acquisition_method = 'Data share on 2020/04/15 by collaborators', priority = PRIOR_JDRUGS)
df_nat = annotate_ds_origin (df_nat, source_name = 'Targets - Bojkova et al.', acquisition_method = 'Processed supplementary data, 2020/05/21', priority = PRIOR_NAT)
df_crispr = annotate_ds_origin (df_crispr, source_name = 'CRISPR - Wei et al.', acquisition_method = 'Processed data from text, 2020/09/04', priority = PRIOR_CRISPR)


# PPI edges
df_ppi = annotate_ds_origin (df_ppi, source_name = 'SarS-CoV-2 - Human Interactome', acquisition_method = 'Preprint/Nature, experimental data', priority = PRIOR_PREPRINT)
df_phipster = annotate_ds_origin (df_phipster, source_name = 'PHIPSTER', acquisition_method = 'Predictions by PHIPSTER webapp', priority = PRIOR_PHIPSTER)
#df_hats = annotate_ds_origin (df_hats, source_name = 'HATs', acquisition_method = 'HATs network by SmartGraph network pharmacology platform, max. dist 3, min.conf 0.15 ', priority = PRIOR_SG)


df_drugs_add =  annotate_ds_origin (df_drugs_add, source_name = 'DTI - Predicted', acquisition_method = 'Data share on 2020/04/15 by collaborators.', priority = PRIOR_JDRUGS)

#print (df_hats.head())

#####
##### 5. Merge PPIs
#####

df_all_ppi = df_ppi.append(df_phipster, ignore_index = True)

#####
##### 6. Modify virus protein names to start with 'SARS-CoV2 ' prefix
#####

df_all_ppi['new_name'] = df_all_ppi.apply (lambda x: modify_virus_protein_name (x['virus_protein']), axis = 1)
df_all_ppi = df_all_ppi.drop(columns = ['virus_protein']).copy()
df_all_ppi = df_all_ppi.rename (columns = {'new_name': 'virus_protein'})



#####
##### 7. Get unique human proteins of PPIs so far
#####


unique_human_proteins = []
unique_virus_proteins = []



for p in list(df_all_ppi['human_protein']):
	unique_human_proteins.append (p)

for p in list(df_taiml['human_protein']):
	unique_human_proteins.append (p)

#for p in list(df_drugs['human_protein']):
#	unique_human_proteins.append (p)

for p in list(df_drugs_add_nodes['human_protein']):
	unique_human_proteins.append (p)

for p in list(df_nat['human_protein']):
	unique_human_proteins.append (p)

for p in list(df_crispr['human_protein']):
	unique_human_proteins.append (p)




unique_human_proteins = list(set(unique_human_proteins))


for p in list(df_all_ppi['virus_protein']):
	unique_virus_proteins.append (p)


unique_virus_proteins = list(set(unique_virus_proteins))

unique_human_proteins = get_uniprot_from_gene_name (unique_human_proteins, is_test)



df = pd.DataFrame ({'uniprot_id': unique_human_proteins})


df_chembl = pd.read_csv (FILE_chembl_uniprot, sep = '\t')

df = df.merge (df_chembl, left_on = 'uniprot_id', right_on = 'uniprot_id', how = 'inner')

df.to_csv (FILE_uhps_prestring, sep = '\t', index = False)

print ()
print ('NEXT STEP: SmartGraph analysis (https://smartgraph.ncats.io).')
print ()

print ("[Done.]")

