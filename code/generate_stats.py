# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advanving Translational Sciences (NCATS/NIH)
#
# References
#

import pandas as pd
from dbutils import *
import py2neo


FILE_DB_CONFIG = '../cfg/db.cfg'
FILE_STATS = '../data/output/db_stats.tsv'
FILE_STATS2 = '../data/output/db_stats2.tsv'


graph =  open_neo4j_connection (FILE_DB_CONFIG)


# Host targets

host_targets = {}


host_targets['preprint'] = graph.run("MATCH (t:Target {is_in_preprint:'True', target_type:'human'}) RETURN count(t) as res").data()[0]['res']
host_targets['taiml'] = graph.run("MATCH (t:Target {is_in_taiml:'True', target_type:'human'}) RETURN count(t) as res").data()[0]['res']
host_targets['hats'] = graph.run("MATCH (t:Target {is_in_hats:'True', target_type:'human'}) RETURN count(t) as res").data()[0]['res']
host_targets['jdti'] = graph.run("MATCH (t:Target {is_in_jdti:'True', target_type:'human'}) RETURN count(t) as res").data()[0]['res']
host_targets['phipster'] = graph.run("MATCH (t:Target {is_in_phipster:'True', target_type:'human'}) RETURN count(t) as res").data()[0]['res']
host_targets['string'] = graph.run("MATCH (t:Target {is_in_string:'True', target_type:'human'}) RETURN count(t) as res").data()[0]['res']
host_targets['natdt'] = graph.run("MATCH (t:Target {is_in_natdt:'True', target_type:'human'}) RETURN count(t) as res").data()[0]['res']
host_targets['dc'] = graph.run("MATCH (t:Target {is_in_drugcentral:'True', target_type:'human'}) RETURN count(t) as res").data()[0]['res']
host_targets['crispr'] = graph.run("MATCH (t:Target {is_in_crispr:'True', target_type:'human'}) RETURN count(t) as res").data()[0]['res']



# Viral targets

viral_targets = {}


viral_targets['preprint'] = graph.run("MATCH (t:Target {is_in_preprint:'True', target_type:'virus'}) RETURN count(t) as res").data()[0]['res']
viral_targets['taiml'] = '-'	# graph.run("MATCH (t:Target {is_in_taiml:'True', target_type:'virus'}) RETURN count(t) as res").data()[0]['res']
viral_targets['hats'] = '-'	# graph.run("MATCH (t:Target {is_in_hats:'True', target_type:'virus'}) RETURN count(t) as res").data()[0]['res']
viral_targets['jdti'] = '-'	# graph.run("MATCH (t:Target {is_in_jdti:'True', target_type:'virus'}) RETURN count(t) as res").data()[0]['res']
viral_targets['phipster'] = graph.run("MATCH (t:Target {is_in_phipster:'True', target_type:'virus'}) RETURN count(t) as res").data()[0]['res']
viral_targets['string'] = '-'	# graph.run("MATCH (t:Target {is_in_string:'True', target_type:'virus'}) RETURN count(t) as res").data()[0]['res']
viral_targets['natdt'] = '-'	# graph.run("MATCH (t:Target {is_in_natdt:'True', target_type:'virus'}) RETURN count(t) as res").data()[0]['res']
viral_targets['dc'] = '-'
viral_targets['crispr'] = '-'





# Compounds

compounds = {}


compounds['dc'] = graph.run("MATCH (c:Compound {is_in_drugcentral:'True'}) RETURN count(c) as res").data()[0]['res']
compounds['jdti'] = graph.run("MATCH (c:Compound {is_in_jdti:'True'}) RETURN count(c) as res").data()[0]['res']
compounds['hats'] = '-'
compounds['taiml'] = '-'
compounds['phipster'] = '-'
compounds['string'] = '-'
compounds['natdt'] = '-'
compounds['preprint'] = '-'
compounds['crispr'] = '-'



# PHIs (HPIs)

phis = {}

phis['preprint'] = graph.run("MATCH (v:Target {target_type:'virus'})-[r:INTERACTS {is_in_preprint:'True'}]->(t:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
phis['taiml'] = '-'	# graph.run("MATCH (v:Target {target_type:'virus'})-[r:INTERACTS {is_in_taiml:'True'}]->(t:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
phis['hats'] = '-'
phis['jdti'] = '-'
phis['phipster'] = graph.run("MATCH (v:Target {target_type:'virus'})-[r:INTERACTS {is_in_phipster:'True'}]->(t:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
phis['string'] = '-'
phis['natdt'] = '-'
phis['dc'] = '-'
phis['crispr'] = '-'





# PPIs (HHIs)

ppis = {}

ppis['preprint'] = '-'	# graph.run("MATCH (t1:Target {target_type:'human'})-[r:INTERACTS {is_in_preprint:'True'}]->(t2:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
ppis['taiml'] = '-'	# graph.run("MATCH (t1:Target {target_type:'human'})-[r:INTERACTS {is_in_taiml:'True'}]->(t2:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
ppis['hats'] = graph.run("MATCH (t1:Target {target_type:'human'})-[r:INTERACTS {is_in_hats:'True'}]->(t2:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
ppis['jdti'] = '-'	# graph.run("MATCH (t1:Target {target_type:'human'})-[r:INTERACTS {is_in_jdti:'True'}]->(t2:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
ppis['phipster'] = '-'	# graph.run("MATCH (t1:Target {target_type:'human'})-[r:INTERACTS {is_in_phipster:'True'}]->(t2:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
ppis['string'] = graph.run("MATCH (t1:Target {target_type:'human'})-[r:INTERACTS {is_in_string:'True'}]->(t2:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
ppis['natdt'] = '-'	# graph.run("MATCH (t1:Target {target_type:'human'})-[r:INTERACTS {is_in_natdt:'True'}]->(t2:Target {target_type:'human'}) RETURN count(r) as res").data()[0]['res']
ppis['dc'] = '-'
ppis['crispr'] = '-'


# DTIs

dtis = {}

dtis['preprint'] = '-' # graph.run("MATCH (c:Compound)-[r:DTI {is_in_preprint:'True'}]->(t:Target) RETURN count(r) as res").data()[0]['res']
dtis['taiml'] = '-'    # graph.run("MATCH (c:Compound)-[r:DTI {is_in_taiml:'True'}]->(t:Target) RETURN count(r) as res").data()[0]['res']
dtis['hats'] = '-'	# graph.run("MATCH (c:Compound)-[r:DTI {is_in_hats:'True'}]->(t:Target) RETURN count(r) as res").data()[0]['res']
dtis['jdti'] = graph.run("MATCH (c:Compound)-[r:DTI {is_in_jdti:'True'}]->(t:Target) RETURN count(r) as res").data()[0]['res']
dtis['phipster'] = '-' # graph.run("MATCH (c:Compound)-[r:DTI {is_in_phipster:'True'}]->(t:Target) RETURN count(r) as res").data()[0]['res']
dtis['string'] = '-'	# graph.run("MATCH (c:Compound)-[r:DTI {is_in_string:'True'}]->(t:Target) RETURN count(r) as res").data()[0]['res']
dtis['natdt'] = '-'    # graph.run("MATCH (c:Compound)-[r:DTI {is_in_natdt:'True'}]->(t:Target) RETURN count(r) as res").data()[0]['res']
dtis['dc'] = graph.run("MATCH (c:Compound)-[r:DTI {is_in_drugcentral:'True'}]->(t:Target) RETURN count(r) as res").data()[0]['res']
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
stats_compounds = []
stats_phis = []
stats_ppis = []
stats_dtis = []
ds_labels = []

for i in range(len(ds)):

	stats_host_targets.append(host_targets[ds[i]])
	stats_viral_targets.append(viral_targets[ds[i]])
	stats_compounds.append(compounds[ds[i]])
	stats_phis.append(phis[ds[i]])
	stats_ppis.append(ppis[ds[i]])
	stats_dtis.append(dtis[ds[i]])
	ds_labels.append(labels[ds[i]])


df = pd.DataFrame ({
	'Dataset': ds_labels,
	'Host Targets': stats_host_targets,
	'Viral Targets': stats_viral_targets,
	'Compounds': stats_compounds,
	'HPIs': stats_phis,
	'HHIs': stats_ppis,
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

stats['nr_targets'] = graph.run("MATCH (t:Target) return count(t) as res").data()[0]['res']
stats['nr_human_targets'] = graph.run("MATCH (t:Target) where t.target_type='human' return count(t) as res").data()[0]['res']
stats['nr_viral_targets'] = graph.run("MATCH (t:Target) where t.target_type='virus' return count(t) as res").data()[0]['res']

stats['nr_compounds'] = graph.run("MATCH (d:Compound) return count(d) as res").data()[0]['res']

stats['nr_hpi'] = graph.run("match ()-[r:INTERACTS]->() where r.interaction_type='HPI' return count(distinct(r.edge_label)) as res").data()[0]['res']
stats['nr_hhi'] = graph.run("match ()-[r:INTERACTS]->() where r.interaction_type='HHI' return count(distinct(r.edge_label)) as res").data()[0]['res']

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
