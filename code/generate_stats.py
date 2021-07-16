# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advancing Translational Sciences (NCATS/NIH)
#
# References
#

import pandas as pd
from dbutils import *
import py2neo

# DB Config

FILE_DB_CONFIG = '../cfg/db.cfg'

# Output

FILE_STATS = '../data/output/db_stats.tsv'
FILE_STATS2 = '../data/output/db_stats2.tsv'


graph =  open_neo4j_connection (FILE_DB_CONFIG)


# Host targets

host_targets = {}


host_targets['preprint'] = graph.run("MATCH (h:HostProtein {is_in_hpi_krogan:True}) RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']
host_targets['taiml'] = graph.run("MATCH (h:HostProtein {is_in_hostprot_taiml:True}) RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']
host_targets['hats'] = graph.run("MATCH (h:HostProtein) where h.is_in_hostprot_hats=True or h.is_in_ppi_sg=True RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']


ht_jm_nhc = graph.run("MATCH (h:HostProtein {is_in_dti_jm_nhc:True}) RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']
ht_jm_hcq = graph.run("MATCH (h:HostProtein {is_in_dti_jm_hcq:True}) RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']
ht_jm_cam = graph.run("MATCH (h:HostProtein {is_in_dti_jm_cam:True}) RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']

host_targets['jdti'] = ht_jm_nhc + ht_jm_hcq + ht_jm_cam

host_targets['phipster'] = graph.run("MATCH (h:HostProtein {is_in_hpi_phipster:True}) RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']
host_targets['string'] = graph.run("MATCH (h:HostProtein {is_in_ppi_string:True}) RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']
host_targets['natdt'] = graph.run("MATCH (h:HostProtein {is_in_hostprot_nat:True}) RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']
host_targets['crispr'] = graph.run("MATCH (h:HostProtein {is_in_hostprot_crispr:True}) RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']

host_targets['dc'] = graph.run("MATCH (d:Drug)-[r:DTI]->(h:HostProtein) where r.is_in_dti_dc=True RETURN count(distinct(h.gene_symbol)) as res").data()[0]['res']



# Viral targets

viral_targets = {}


viral_targets['preprint'] = graph.run("MATCH (p:PathogenProtein {is_in_hpi_krogan:True}) RETURN count(distinct(p.name)) as res").data()[0]['res']
viral_targets['taiml'] = '-'
viral_targets['hats'] = '-'	
viral_targets['jdti'] = '-'
viral_targets['phipster'] = graph.run("MATCH (p:PathogenProtein {is_in_hpi_phipster:True}) RETURN count(distinct(p.name)) as res").data()[0]['res']
viral_targets['string'] = '-'
viral_targets['natdt'] = '-'
viral_targets['dc'] = '-'
viral_targets['crispr'] = '-'





# Drugs

drugs = {}


drugs['dc'] = graph.run("MATCH (d:Drug {is_in_dti_dc:True}) RETURN count(distinct(d.drug_name)) as res").data()[0]['res']

nr_drug_nhc = graph.run("MATCH (d:Drug {is_in_dti_jm_nhc:True}) RETURN count(distinct(d.drug_name)) as res").data()[0]['res']
nr_drug_hcq = graph.run("MATCH (d:Drug {is_in_dti_jm_hcq:True}) RETURN count(distinct(d.drug_name)) as res").data()[0]['res']
nr_drug_cam = graph.run("MATCH (d:Drug {is_in_dti_jm_cam:True}) RETURN count(distinct(d.drug_name)) as res").data()[0]['res']

drugs['jdti'] = nr_drug_nhc + nr_drug_hcq + nr_drug_cam

drugs['hats'] = '-'
drugs['taiml'] = '-'
drugs['phipster'] = '-'
drugs['string'] = '-'
drugs['natdt'] = '-'
drugs['preprint'] = '-'
drugs['crispr'] = '-'



# HPIs

hpis = {}

hpis['preprint'] = graph.run("MATCH (p:PathogenProtein)-[r:HPI {is_in_hpi_krogan:True}]->(h:HostProtein) RETURN count(distinct(r.edge_label)) as res").data()[0]['res']
hpis['taiml'] = '-'	
hpis['hats'] = '-'
hpis['jdti'] = '-'
hpis['phipster'] = graph.run("MATCH (p:PathogenProtein)-[r:HPI {is_in_hpi_phipster:True}]->(h:HostProtein) RETURN count(r) as res").data()[0]['res']
hpis['string'] = '-'
hpis['natdt'] = '-'
hpis['dc'] = '-'
hpis['crispr'] = '-'





# PPIs 

ppis = {}

ppis['preprint'] = '-'
ppis['taiml'] = '-'
ppis['hats'] = graph.run("MATCH (t1:HostProtein)-[r:PPI {is_in_ppi_sg:True}]->(t2:HostProtein) RETURN count(distinct(r.edge_label)) as res").data()[0]['res']
ppis['jdti'] = '-'
ppis['phipster'] = '-'
ppis['string'] = graph.run("MATCH (t1:HostProtein)-[r:PPI {is_in_ppi_string:True}]->(t2:HostProtein) RETURN count(distinct(r.edge_label)) as res").data()[0]['res']
ppis['natdt'] = '-'
ppis['dc'] = '-'
ppis['crispr'] = '-'



# DTIs

dtis = {}

dtis['preprint'] = '-'
dtis['taiml'] = '-'
dtis['hats'] = '-'

jdti_nhc = graph.run("MATCH (d:Drug)-[r:DTI {is_in_dti_jm_nhc:True}]->(h:HostProtein) RETURN count(distinct(r.edge_label)) as res").data()[0]['res']

jdti_hcq = graph.run("MATCH (d:Drug)-[r:DTI {is_in_dti_jm_hcq:True}]->(h:HostProtein) RETURN count(distinct(r.edge_label)) as res").data()[0]['res']

jdti_cam = graph.run("MATCH (d:Drug)-[r:DTI {is_in_dti_jm_cam:True}]->(h:HostProtein) RETURN count(distinct(r.edge_label)) as res").data()[0]['res']

dtis['jdti'] = jdti_nhc + jdti_hcq + jdti_cam


dtis['phipster'] = '-'
dtis['string'] = '-'
dtis['natdt'] = '-'
dtis['dc'] = graph.run("MATCH (d:Drug)-[r:DTI {is_in_dti_dc:True}]->(h:HostProtein) RETURN count(distinct(r.edge_label)) as res").data()[0]['res']
dtis['crispr'] = '-'

ds = ['preprint', 'natdt', 'taiml', 'phipster', 'hats', 'jdti', 'string', 'dc', 'crispr']
idx = [6, 1, 3, 7, 5, 8, 4, 9, 2]


labels = {
	'preprint': 'Interactome Study',
	'taiml': 'Meta Path AI/ML',
	'hats': 'SmartGraph / HATs',
	'jdti': 'Predicted DTIs',
	'phipster': 'P-HIPSter',
	'string': 'STRING',
	'natdt': 'Proteomics Study',
	'dc': 'DrugCentral',
	'crispr': 'CRISPR' 
}

stats_host_targets = []
stats_viral_targets = []
stats_drugs = []
stats_hpis = []
stats_ppis = []
stats_dtis = []
ds_labels = []

for i in range(len(ds)):

	stats_host_targets.append(host_targets[ds[i]])
	stats_viral_targets.append(viral_targets[ds[i]])
	stats_drugs.append(drugs[ds[i]])
	stats_hpis.append(hpis[ds[i]])
	stats_ppis.append(ppis[ds[i]])
	stats_dtis.append(dtis[ds[i]])
	ds_labels.append(labels[ds[i]])


df = pd.DataFrame ({
	'Dataset': ds_labels,
	'Host Targets': stats_host_targets,
	'Viral Targets': stats_viral_targets,
	'Drugs': stats_drugs,
	'HPIs': stats_hpis,
	'PPIs': stats_ppis,
	'DTIs': stats_dtis,
	'order': idx
})

df = df.sort_values (by = ['order'])
df = df.drop (columns = ['order']).copy()

#print (host_targets)
#print (viral_targets)
#print (compounds)
#print (phis)
#print (ppis)
#print (dtis)

print (df)

df.to_csv (FILE_STATS, sep = '\t', index = False)



stats = {}

stats['nr_human_targets'] = graph.run("MATCH (h:HostProtein) return count(distinct(h.gene_symbol)) as res").data()[0]['res']
stats['nr_viral_targets'] = graph.run("MATCH (p:PathogenProtein) return count(distinct(p.name)) as res").data()[0]['res']
stats['nr_proteins'] = stats['nr_human_targets'] + stats['nr_viral_targets']



stats['nr_drugs'] = graph.run("MATCH (d:Drug) return count(distinct(d.drug_name)) as res").data()[0]['res']

stats['nr_hpi'] = graph.run("match ()-[r:HPI]->() return count(distinct(r.edge_label)) as res").data()[0]['res']
stats['nr_ppi'] = graph.run("match ()-[r:PPI]->() return count(distinct(r.edge_label)) as res").data()[0]['res']
stats['nr_dti'] = graph.run("match ()-[r:DTI]->() return count(distinct(r.edge_label)) as res").data()[0]['res']


cats = []
vals = []

for k in stats.keys():
	cats.append(k)

for k in cats:
	vals.append(stats[k])

df_stat = pd.DataFrame ({
		'category': cats,
		'number': vals
	})

df_stat.to_csv (FILE_STATS2, sep = '\t', index = False)

print ()

print (df_stat)

print ('[Done.]')
