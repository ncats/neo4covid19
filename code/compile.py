# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advanving Translational Sciences (NCATS/NIH)
#
# References
#

# Ref: https://stackoverflow.com/questions/4470523/create-a-branch-in-git-from-another-branch
# Ref: https://git-scm.com/docs/git-branch
# Ref: https://git-scm.com/docs/git-checkout
# Ref: https://github.github.com/training-kit/downloads/github-git-cheat-sheet.pdf

import pandas as pd
from stringify import stringify
from idgtdl import get_tdl_by_gene_name

from ppireactome import get_reactome_ppi

from drugcentral import get_dc_dti, annotate_drug_structures, annotate_drug_names, identify_drugs_by_name

from overlayreactome import overlay_reactome_ppi

from harmonize_drugs import harmonize_dtis

from build_neo4j import build_neo4j

from sg import process_sg

from crispr import load_crispr

import sys

### Input files

FILE_ppi_krogan = '../data/input/SARS-CoV2.tsv'
FILE_ppi_phipster = '../data/input/PHIPSTER_SARS_COV2.xlsx'
FILE_predicted_targets_taiml = '../data/input/Merged.xlsx'

FILE_HATS_SmartGraph = '../data/input/SG_HATs_dist_3_conf_0.00.json'
FILE_HATS_reverse_SmartGraph = '../data/input/SG_HATs_reverse_dist_3_conf_0.00.json'

#FILE_HATS_SmartGraph = '../data/input/HATS_PPI_dist_3_conf_0.15.smartgraph.edges.and.attributes.tab'


# Shared by Jordi Mestres, PhD and Tudor I. Oprea, MD, PhD
FILE_DTI_by_collaborator = '../data/input/COVID19_TargetList_ChemSpace_mod_GZK.xlsx'

FILE_nat_drugtargets = '../data/input/nat_drugtargets.tsv'
FILE_dbcfg = '../cfg/db.cfg'


## Output files
FILE_edges_ppi ='../data/output/edges_ppi.tab'
FILE_nodes_proteins = '../data/output/nodes_proteins.tab'
FILE_edges_DTI = '../data/output/edges_dtis.tab'
FILE_nodes_drugs = '../data/output/nodes_drugs.tab'


### Data Priorities

PRIOR_DCDRUGS = 0
PRIOR_PREPRINT = 1
PRIOR_NAT = 2
PRIOR_CRISPR = 3
PRIOR_TAIML = 4
PRIOR_PHIPSTER = 5
PRIOR_STRING = 6
PRIOR_SG = 7
PRIOR_JDRUGS = 8

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



def process_hats_nw (fname):

	df = process_sg (fname)

	#df['interaction'] = 'undefined' 

	df['metadata'] = ''
	df['comment'] = ''
	df['pmids'] = ''

	df = df[['source_node', 'target_node', 'interaction', 'mechanism', 'source_specific_score', 'metadata', 'pmids', 'comment']].copy()

	return (df)



def update_unique_proteins (df, unique_proteins):
	for p in list(df['source_node']):
		unique_proteins.append (p)

	for p in list(df['target_node']):
		unique_proteins.append (p)

	unique_proteins = list(set(unique_proteins))

	return (unique_proteins)



def create_edge_label (sn, tn):
	return (sn + '_' + tn)



def edges_present_in_dataset (df, df_all, dataset_name, new_column_name):
	df_all = df_all[df_all['source'] == dataset_name]
	reference_edges = list(df_all['edge_label'])
	
	def is_edge_present (edge, reference_edges):
		if edge in reference_edges:
			return (True)
		
		return (False)

	df[new_column_name] = df.apply (lambda x: is_edge_present(x['edge_label'], reference_edges), axis = 1)

	return (df) 	

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



def add_conditional_metadata (df, cond_col, cond_val, md_col, md_key, md_val):
	def set_val (cond_col, cond_val, md_orig, md_key, md_val):
                md_rec = ''
                if cond_col == cond_val:
                        if len(md_orig) == 0:
                                        md_rec = "%s:%s" % (md_key, md_val)
                        else:
                                        md_rec = md_orig
                                        md_rec += ";%s:%s" % (md_key, md_val)
                else:   
                        md_rec = md_orig

                return (md_rec)

	df[md_col] = df.apply (lambda x: set_val(x[cond_col], cond_val, x[md_col], md_key, md_val), axis = 1)

	return (df)






def add_node_metadata (df, node_id, df_met, md_node_id_col, md_col, md_orig, md_key):
	def set_node_val (node_id, df, md_node_id_col, md_col, md_orig, md_key):
                md_rec = ''

                if node_id in list(df[md_node_id_col]):
                        df = df[df[md_node_id_col] == node_id].copy()
                        md_val = list(df[md_col])[0]



                        if len(md_orig) == 0:
                                        md_rec = "%s:%s" % (md_key, md_val)
                        else:
                                        md_rec = md_orig
                                        md_rec += ";%s:%s" % (md_key, md_val)
                else:   
                        md_rec = md_orig

                return (md_rec)

	df[md_orig] = df.apply(lambda x: set_node_val (x[node_id], df_met, md_node_id_col, md_col, x[md_orig], md_key), axis = 1)

	return (df)


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

# The forward and reverse networks will be imported and harmonized the same time in the "Harmonization" section the same.

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
df_phipster['inxtype'] = 'HPI'

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
df_nat_orig = df_nat.copy()
df_nat = df_nat[['human_protein']].copy()

df_crispr = df_crispr.rename (columns = {'gene_symbol': 'human_protein'})

df_hats = process_hats_nw (FILE_HATS_SmartGraph)
df_hats['inxtype'] = 'HHI'

df_hats_rev = process_hats_nw (FILE_HATS_reverse_SmartGraph)
df_hats_rev['inxtype'] = 'HHI'

df_hats = df_hats.append(df_hats_rev, ignore_index = True)

df_hats = df_hats.groupby(['source_node', 'target_node'], as_index = False).aggregate ({
										'interaction': 'first',
										'mechanism': 'first',
										'source_specific_score': 'first',
										'metadata': 'first',
										'pmids': 'first',
										'comment': 'first',
										'inxtype': 'first'
									})





df_ppi['interaction'] = 'undefined'
df_ppi['metadata'] = ''
df_ppi['comment'] = ''

df_ppi['inxtype'] = 'HPI'

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
df_hats = annotate_ds_origin (df_hats, source_name = 'HATs', acquisition_method = 'HATs network by SmartGraph network pharmacology platform, max. dist 3, min.conf 0.00 ', priority = PRIOR_SG)


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






#####
##### 8. STRINGify network
#####


print ('[*] STRINGifying network via STRING API, this may take a while, please hang on ...')

#print ("Nr of unique human protiens for stringifying: %d" % (len(unique_human_proteins)))
#print ()

df_string = stringify (unique_human_proteins, species_ncbi = 9606, limit_of_mapped_genes = 1, max_interactor = 100, score_cutoff = 0, alpha = 0.5, request_id = 'Neo4COVID19', comment = '', priority = PRIOR_STRING)
df_string['inxtype'] = 'HHI'


#exit()

print ('[*] STRINGifying done.')

sources = list(df_string.source_node)
targets = list(df_string.target_node)
all_targets = {}

for n in sources:
	all_targets[n] = 1
for n in targets:
	all_targets[n] = 1

targets = all_targets.keys()


#print ("Nr of unique human protiens after stringifying: %d" % (len(targets)))

#exit()


#print (df_string.shape)


## Updating unique human protein list

unique_human_proteins = update_unique_proteins(df_string, unique_human_proteins)



## Direction of virus-human protein interactions are not known at this point, for sake of consistency virus proteins were appointed as source nodes. GZK

df_all_ppi = df_all_ppi.rename (columns = {
		'virus_protein': 'source_node',
		'human_protein': 'target_node'
	})



df_all_ppi = df_all_ppi [['source_node', 'target_node', 'interaction', 'mechanism', 'source_specific_score', 'source', 'data_origin', 'priority', 'metadata', 'pmids', 'comment', 'inxtype']].copy()

df_all_ppi = df_all_ppi.append(df_string, ignore_index = True)



#####
##### 9. Add aditional networks without STRINGifying them
#####

# HATs network by Tudor, Tim Wilson and SmartGraph

df_all_ppi = df_all_ppi.append(df_hats, ignore_index = True)


## Updating unique human protein list

unique_human_proteins = update_unique_proteins(df_hats, unique_human_proteins)



#####
##### 10. Deduplicating edges by highest priority
#####


df_all_ppi = df_all_ppi.sort_values(['priority'])

df_all_ppi['relationship'] = 'INTERACTS'

#df_all_ppi.to_csv ('../out.tab', sep = '\t', index = False)

df_all_ppi['edge_label'] = df_all_ppi.apply (lambda x: create_edge_label(x['source_node'], x['target_node']), axis = 1)

df_all_unique_ppi = df_all_ppi.groupby(['source_node', 'target_node'], as_index = False).aggregate ({
											'interaction': 'first',
											'mechanism': 'first',
											'source_specific_score': 'first',
											'source': 'first',
											'data_origin': 'first',
											'priority': 'first',
											'metadata': 'first',
											'pmids': 'first',
											'comment': 'first',
											'relationship': 'first',
											'edge_label': 'first',
											'inxtype': 'first'
										})

df_all_unique_ppi = edges_present_in_dataset (df_all_unique_ppi, df_all_ppi, 'SarS-CoV-2 - Human Interactome', 'is_in_preprint')
df_all_unique_ppi = edges_present_in_dataset (df_all_unique_ppi, df_all_ppi, 'PHIPSTER', 'is_in_phipster')
df_all_unique_ppi = edges_present_in_dataset (df_all_unique_ppi, df_all_ppi, 'STRING', 'is_in_string')
df_all_unique_ppi = edges_present_in_dataset (df_all_unique_ppi, df_all_ppi, 'HATs', 'is_in_hats')







#print (df_all_unique_ppi.head())

#df_all_unique_ppi.to_csv ('../holistic.tab', sep = '\t', index = False)



#####
##### 11. Track node provenance
#####


df_node_provenance = record_nodes (df_all_ppi, 'SarS-CoV-2 - Human Interactome', df_node_provenance)


df_node_provenance = record_nodes_from_nodes_only_data_source (df_taiml, 'human_protein', 'Meta Path AIML', df_node_provenance)

df_node_provenance = record_nodes (df_all_ppi, 'PHIPSTER', df_node_provenance)

df_node_provenance = record_nodes (df_all_ppi, 'STRING', df_node_provenance)

df_node_provenance = record_nodes (df_all_ppi, 'HATs', df_node_provenance)

df_node_provenance = record_nodes_from_nodes_only_data_source (df_drugs, 'human_protein', 'DrugCentral', df_node_provenance)

df_node_provenance = record_nodes_from_nodes_only_data_source (df_drugs_add_nodes, 'human_protein', 'DTI - Predicted', df_node_provenance)

df_node_provenance = record_nodes_from_nodes_only_data_source (df_nat, 'human_protein', 'Targets - Bojkova et al.', df_node_provenance)

df_node_provenance = record_nodes_from_nodes_only_data_source (df_crispr, 'human_protein', 'CRISPR - Wei et al.', df_node_provenance)


#print (df_node_provenance)

#print (df_drugs.columns)



#####
##### 12. Add TDL categories
#####

df_human_proteins = pd.DataFrame ({'gene': unique_human_proteins})

df_human_proteins = get_tdl_by_gene_name (df_human_proteins, is_test)
df_human_proteins['target_type'] = 'human'
df_human_proteins['node_type'] = 'target'


print (df_human_proteins.head())



df_virus_proteins = pd.DataFrame ({'gene': unique_virus_proteins})
df_virus_proteins['tdl'] = 'virus'
df_virus_proteins['target_type'] = 'virus'
df_virus_proteins['node_type'] = 'target'

#####
##### 13. Finalize node information
#####


df_human_proteins = is_node_present_in_dataset (df_human_proteins, 'SarS-CoV-2 - Human Interactome', df_node_provenance, 'is_in_preprint') 
df_human_proteins = is_node_present_in_dataset (df_human_proteins, 'Meta Path AIML', df_node_provenance, 'is_in_taiml') 
df_human_proteins = is_node_present_in_dataset (df_human_proteins, 'PHIPSTER', df_node_provenance, 'is_in_phipster') 
df_human_proteins = is_node_present_in_dataset (df_human_proteins, 'STRING', df_node_provenance, 'is_in_string') 
df_human_proteins = is_node_present_in_dataset (df_human_proteins, 'HATs', df_node_provenance, 'is_in_hats')
df_human_proteins = is_node_present_in_dataset (df_human_proteins, 'DTI - Predicted', df_node_provenance, 'is_in_jdti')
df_human_proteins = is_node_present_in_dataset (df_human_proteins, 'Targets - Bojkova et al.', df_node_provenance, 'is_in_natdt')
df_human_proteins = is_node_present_in_dataset (df_human_proteins, 'DrugCentral', df_node_provenance, 'is_in_drugcentral')
df_human_proteins = is_node_present_in_dataset (df_human_proteins, 'CRISPR - Wei et al.', df_node_provenance, 'is_in_crispr')





df_virus_proteins = is_node_present_in_dataset (df_virus_proteins, 'SarS-CoV-2 - Human Interactome', df_node_provenance, 'is_in_preprint') 
df_virus_proteins = is_node_present_in_dataset (df_virus_proteins, 'Meta Path AIML', df_node_provenance, 'is_in_taiml') 
df_virus_proteins = is_node_present_in_dataset (df_virus_proteins, 'PHIPSTER', df_node_provenance, 'is_in_phipster') 
df_virus_proteins = is_node_present_in_dataset (df_virus_proteins, 'STRING', df_node_provenance, 'is_in_string') 
df_virus_proteins = is_node_present_in_dataset (df_virus_proteins, 'HATs', df_node_provenance, 'is_in_hats')
df_virus_proteins = is_node_present_in_dataset (df_virus_proteins, 'DTI - Predicted', df_node_provenance, 'is_in_jdti')
 

print (df_human_proteins.head())
print (df_virus_proteins.head())

#df_human_proteins.to_csv ('../out_human_proteins.tab', sep = '\t', index = False)


#####
##### 14. Overlay Reactome data
#####

## 1. Filter Reactom to keep only edges that are in holistic network.
## 2. Flag edges that are in Reactome
## 3. Get regdir and annotation from Reactome
## 4. Create reverse edges for those that are in Reactome but redir is undefined, preserve Reactome annotation

df_reactome = get_reactome_ppi()

#df_reactome.to_csv ('../reactome.tab', sep = '\t', index = False)

df_all_unique_ppi = overlay_reactome_ppi (df_all_unique_ppi, df_reactome)

df_all_unique_ppi['comment'] = df_all_unique_ppi['comment'].fillna('')
df_all_unique_ppi['metadata'] = df_all_unique_ppi['metadata'].fillna('')
df_all_unique_ppi['edge_label'] = df_all_unique_ppi.apply (lambda x: create_edge_label (x['source_node'], x['target_node']), axis = 1)



#####
##### 15. Process DTIs and Drug nodes
#####

#print (df_drugs.columns)

#print (df_human_proteins.columns)

df_drugs = df_drugs[df_drugs['human_protein'].isin(list(df_human_proteins['gene']))]

df_drugs_orig = df_drugs



#df_drugs_orig.to_csv ('../out_orig_dc_dti.tab', sep ='\t', index = False)



#df_drugs_add.to_csv ('../out_jdti_all.tab', sep = '\t', index = False)




df_drugs_jdti = identify_drugs_by_name (df_drugs_add, 'Known binders')

df_drugs_jdti = df_drugs_jdti.rename (columns = {
		'STRUCT_ID': 'struct_id'
	})



df_drugs_jdti = df_drugs_jdti[['struct_id', 'human_protein', 'Chemical space', 'source', 'data_origin', 'priority', 'Known binders']].copy()




#df_drugs_jdti.to_csv ('../out_jdti.tab', sep = '\t', index = False)
#df_drugs.to_csv ('../out_dcdti.tab', sep = '\t', index = False)

# Harmonize drug sets

# For JDTI add missing columns ('p_chembl', etc as unkown), then add back those jDTIs that are not in df_drugs_jdti



# Concatenate DC DTI and complete JDTI
# Deduplicate existing DTI by priority

(df_dti_all, df_dti_all_aggr) = harmonize_dtis (df_drugs_jdti, df_drugs)





# Record edge and node provenance
#df_dti_all = df_dti_all
df_dti = df_dti_all_aggr

df_dti = is_edge_in_drug_dataset(df_dti, df_dti_all, 'DrugCentral', 'is_in_drugcentral')
df_dti = is_edge_in_drug_dataset(df_dti, df_dti_all, 'DTI - Predicted', 'is_in_jdti')



df_dti['comment'] = ''
df_dti['pmids'] = ''
df_dti['metadata'] = ''

df_drugs = df_dti.groupby('struct_id', as_index = False).aggregate ({'drug_name': 'first'})

df_drugs = annotate_drug_names(df_drugs, 'left')



#exit()


df_drugs['drug_name'] = df_drugs['drug_name'].fillna('')

def fix_drug_name (dn, inn):
	if dn == '':
		return (inn)

	return (dn)

df_drugs['new_drug_name'] = df_drugs.apply (lambda x: fix_drug_name (x['drug_name'], x['DRUG_NAME']), axis = 1)
df_drugs = df_drugs.drop (columns = ['drug_name', 'DRUG_NAME', 'STRUCT_ID']).copy()
df_drugs = df_drugs.rename (columns = {'new_drug_name': 'drug_name'})


df_drugs = annotate_drug_structures (df_drugs, 'left')




def is_in_drugs (df, df_all, ds, new_column):
	df_all = df_all[df_all['source'] == ds]

	nodes = list(set(list(df_all['struct_id'])))
	df.loc[df['struct_id'].isin(nodes), new_column] = True
	df.loc[~df['struct_id'].isin(nodes), new_column] = False

	return (df)



df_drugs = is_node_in_drug_dataset (df_drugs, df_dti_all, 'DrugCentral', 'is_in_drugcentral')
df_drugs = is_node_in_drug_dataset (df_drugs, df_dti_all, 'DTI - Predicted', 'is_in_jdti')





df_drugs = is_in_chemical_space_dataset (df_drugs, df_drugs_add, 'HCQ', 'is_in_hcq')
df_drugs = is_in_chemical_space_dataset (df_drugs, df_drugs_add, 'NHC', 'is_in_nhc')
df_drugs = is_in_chemical_space_dataset (df_drugs, df_drugs_add, 'CAM', 'is_in_cam')
df_drugs = is_in_drugs (df_drugs, df_dti_all, 'DrugCentral', 'is_in_drugs')

df_drugs['node_type'] = 'Compound'



df_dti = df_dti.rename (columns = {
		'struct_id': 'source_node',
		'human_protein': 'target_node'
	})


nr_prior = df_dti.shape[0]
#df_dti = df_dti.drop (columns = ['chemical_space'])
df_j = df_drugs [['struct_id', 'drug_name']].copy()
df_dti = df_dti.drop(columns = ['drug_name'])




df_dti = df_dti.merge (df_j, left_on = 'source_node', right_on = 'struct_id', how = 'inner')


#exit()


nr_post = df_dti.shape[0]

if nr_prior != nr_post:
	print ('[ERROR] 1:1 mapping violated when joning the final drugs and dti data. Teminating ...')
	sys.exit (-7)


print (df_dti)

print (df_drugs)



# DTI: df_dti
# PPI: df_all_unique_ppi

df_proteins = df_human_proteins
df_proteins = df_proteins.append(df_virus_proteins, ignore_index = True)


#####
##### 16. Add metadata
#####

#df_all_unique_ppi['metadata'] = ''
#df_dti['metadata'] = ''
df_proteins['metadata'] = ''
df_drugs['metadata'] = ''


#### Edge direction
df_all_unique_ppi = add_conditional_metadata (df=df_all_unique_ppi, cond_col='is_in_hats', cond_val=True, md_col='metadata', md_key='orig_direction_hats', md_val='directed')
df_all_unique_ppi = add_conditional_metadata (df=df_all_unique_ppi, cond_col='is_in_string', cond_val=True, md_col='metadata', md_key='orig_direction_string', md_val='undirected')
df_dti = add_conditional_metadata (df=df_dti, cond_col='is_in_drugcentral', cond_val=True, md_col='metadata', md_key='orig_direction_drugcentral', md_val='manually_assigned_drug_to_target')
df_dti = add_conditional_metadata (df=df_dti, cond_col='is_in_jdti', cond_val=True, md_col='metadata', md_key='orig_direction_jdti', md_val='manually_assigned_drug_to_target')
df_all_unique_ppi = add_conditional_metadata (df=df_all_unique_ppi, cond_col='is_in_preprint', cond_val=True, md_col='metadata', md_key='orig_direction_preprint', md_val='manually_assigned_virus_to_host')
df_all_unique_ppi = add_conditional_metadata (df=df_all_unique_ppi, cond_col='is_in_phipster', cond_val=True, md_col='metadata', md_key='orig_direction_phipster', md_val='manually_assigned_virus_to_host')

#### Notes for CRIPSR data

df_proteins = add_node_metadata (df=df_proteins, node_id='gene', df_met=df_crispr, md_node_id_col='human_protein', md_col='notes', md_orig='metadata', md_key='notes_crispr')
df_proteins = add_node_metadata (df=df_proteins, node_id='gene', df_met=df_crispr, md_node_id_col='human_protein', md_col='role', md_orig='metadata', md_key='role_crispr')
df_proteins = add_node_metadata (df=df_proteins, node_id='gene', df_met=df_nat_orig, md_node_id_col='human_protein', md_col='Quantity_Change', md_orig='metadata', md_key='quantity_change')





#####
##### 17. Build Neo4j
#####

df_all_unique_ppi.to_csv (FILE_edges_ppi, sep = '\t', index = False)
df_proteins.to_csv (FILE_nodes_proteins, sep = '\t', index = False)
df_dti.to_csv (FILE_edges_DTI, sep = '\t', index = False)
df_drugs.to_csv (FILE_nodes_drugs, sep = '\t', index = False)


build_neo4j (FILE_dbcfg, df_proteins, df_drugs, df_all_unique_ppi, df_dti)

## Deduplicate DTIs

#print (df_ppi.head())
#print (df_phipster.head())
#print (df_taiml.head())
#print (df_map.head())
#print (df_all_ppi.head())
#print (df_string.head())
