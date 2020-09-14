# Ref: https://stackoverflow.com/questions/27884268/return-pandas-dataframe-from-postgresql-query-with-sqlalchemy
# Ref: http://initd.org/psycopg/docs/module.html
# Ref: https://gist.github.com/00krishna/9026574
# Ref: https://medium.com/neo4j/py2neo-v4-2bedc8afef2
# Ref: https://py2neo.org/v4/
# Ref: https://stackoverflow.com/questions/43711680/i-want-to-return-the-edge-list-of-a-a-path-using-cypher-and-neo4j-how-to-do-it
# Ref: https://pandas.pydata.org/pandas-docs/version/0.18/generated/pandas.DataFrame.empty.html
# Ref: https://neo4j.com/docs/cypher-manual/current/clauses/delete/
# Ref: https://neo4j.com/docs/cypher-manual/current/clauses/create/
# Ref: https://stackoverflow.com/questions/21979782/how-to-push-values-to-property-array-cypher-neo4j
# Ref: https://github.com/ncats/smartgraph_backend


#from sqlalchemy import create_engine
#import psycopg2 as pg
from py2neo import Graph
from dfply import *


### Neo4j ####




def read_neo4j_config(db_config_file):
        fp = open (db_config_file, 'r')

        db_par = {}

        db_par['host'] = fp.readline().strip()
        db_par['port'] = fp.readline().strip()
        db_par['user'] = fp.readline().strip()
        db_par['pwd'] = fp.readline().strip()

        fp.close()

        #print (db_par)

        return (db_par)



def open_neo4j_connection (db_config_file):

	db_par = read_neo4j_config(db_config_file)

	
	host_url = "bolt://" + db_par['host'] + ":" + db_par['port']
	
	graph = Graph(host_url, auth=(db_par['user'], db_par['pwd']))
        
	return (graph)
	



def erase_db (conn):
	query = 'MATCH (n) DETACH DELETE n'
	
	conn.run(query)


def create_target_node(conn, gene, is_in_hats, is_in_jdti, is_in_phipster, is_in_preprint, is_in_crispr, is_in_string, is_in_taiml, is_in_natdt, is_in_drugcentral, node_type, target_type, tdl, uniprot, uuid, metadata):
	query = "CREATE (t:Target {node_type: '%s', " % (node_type)
	query += "gene_symbol: '%s', " % (gene)
	query += "target_type: '%s', " % (target_type)
	query += "tdl: '%s', " % (tdl)
	query += "uniprot: '%s', " % (uniprot)
	query += "is_in_preprint: '%s', " % (str(is_in_preprint))	
	query += "is_in_crispr: '%s', " % (str(is_in_crispr))	
	query += "is_in_phipster: '%s', " % (str(is_in_phipster))	
	query += "is_in_taiml: '%s', " % (str(is_in_taiml))
	query += "is_in_drugcentral: '%s', " % (str(is_in_drugcentral))
	query += "is_in_jdti: '%s', " % (str(is_in_jdti))
	query += "is_in_string: '%s', " % (str(is_in_string))
	query += "is_in_hats: '%s', " % (str(is_in_hats))
	query += "is_in_natdt: '%s', " % (str(is_in_natdt))
	query += "metadata: '%s', " % (metadata)
	query += "uuid: '%s'})" % (uuid)

	
	
	#print (query)
	
	conn.run(query)


def create_target_indices (conn):
	query = 'CREATE CONSTRAINT ON (t:Target) ASSERT t.gene_symbol IS UNIQUE'

	#print (query)
	
	conn.run(query)

	query = 'CREATE CONSTRAINT ON (t:Target) ASSERT t.uuid IS UNIQUE'

	#print (query)
	
	conn.run(query)


def create_compound_node(conn, struct_id, drug_name, SMILES, InChI, InChIKey, CAS_RN, is_in_drugcentral, is_in_jdti, is_in_hcq, is_in_nhc, is_in_cam, is_in_drugs, node_type, uuid, metadata):
	query = "CREATE (c:Compound {node_type: '%s', " % (node_type)
	query += "struct_id: '%s', " % (struct_id)
	query += "drug_name: '%s', " % (drug_name)
	query += "smiles: '%s', " % (SMILES.replace ('\\', '\\\\'))
	query += "inchi: '%s', " % (InChI.replace ('\\', '\\\\'))
	query += "inchi_key: '%s', " % (InChIKey)
	query += "ns_inchi_key: '%s', " % (InChIKey.split('-')[0].strip())
	query += "CAS_RN: '%s', " % (CAS_RN)	
	query += "is_in_drugcentral: '%s', " % (str(is_in_drugcentral))	
	query += "is_in_jdti: '%s', " % (str(is_in_jdti))
	query += "is_in_hcq: '%s', " % (str(is_in_hcq))
	query += "is_in_nhc: '%s', " % (str(is_in_nhc))
	query += "is_in_cam: '%s', " % (str(is_in_cam))
	query += "is_in_drugs: '%s', " % (str(is_in_drugs))
	query += "metadata: '%s', " % (metadata)
	query += "uuid: '%s'})" % (uuid)

	
	
	#print (query)
	
	conn.run(query)




def create_compound_indices (conn):
	query = 'CREATE CONSTRAINT ON (c:Compound) ASSERT c.drug_name IS UNIQUE'

	#print (query)
	
	conn.run(query)

	query = 'CREATE CONSTRAINT ON (c:Compound) ASSERT c.uuid IS UNIQUE'

	#print (query)
	
	conn.run(query)

	query = 'CREATE CONSTRAINT ON (c:Compound) ASSERT c.struct_id IS UNIQUE'

	#print (query)
	
	conn.run(query)	


def create_ppi_edge (conn, source_node, target_node, interaction, mechanism, reactome_mechanism, reactome_regdir, reactome_score, metadata, comment, pmids, priority, source, data_origin, relationship, source_specific_score, is_in_preprint, is_in_phipster, is_in_string, is_in_hats, is_in_reactome, edge_label, source_node_uuid, target_node_uuid, uuid, inxtype):
	query = "MATCH (t1:Target { gene_symbol: '" + source_node + "' }) MATCH (t2:Target { gene_symbol: '" + target_node + "' }) MERGE (t1)-[:INTERACTS "
	query += "{ edge_label: '%s', " % (edge_label)
	
	query += "source_node: '%s', " % (source_node) 
	query += "target_node: '%s', " % (target_node)
	query += "interaction_type: '%s', " % (inxtype)
	query += "interaction: '%s', " % (interaction)
	query += "mechanism: '%s', " % (mechanism)
	query += "reactome_mechanism: '%s', " % (reactome_mechanism)
	query += "reactome_regdir: '%s', " % (reactome_regdir)
	query += "reactome_score: '%.2f', " % (reactome_score)
	query += "metadata: '%s', " % (metadata)
	query += "comment: '%s', " % (comment)
	query += "pmids: '%s', " % (pmids)
	query += "priority: '%d', " % (priority)
	query += "source: '%s', " % (source)
	query += "data_origin: '%s', " % (data_origin)
	query += "relationship: '%s', " % (relationship)
	query += "source_specific_score: '%s', " % (source_specific_score)
	query += "is_in_preprint: '%s', " % (str(is_in_preprint))
	query += "is_in_phipster: '%s', " % (str(is_in_phipster))
	query += "is_in_string: '%s', " % (str(is_in_string))
	query += "is_in_hats: '%s', " % (str(is_in_hats))
	query += "is_in_reactome: '%s', " % (str(is_in_reactome))
	query += "source_node_uuid: '%s', " % (source_node_uuid) 
	query += "target_node_uuid: '%s', " % (target_node_uuid)
	query += "uuid: '%s' }]->(t2)" % (uuid)
															

	#print (query)
	
	conn.run(query)
	
def create_dti_edge (conn, action_type, data_origin, drug_name, target_node, is_activity_known, p_chembl, priority, source, source_node, edge_label, relationship, is_in_drugcentral, is_in_jdti, comment, pmids, metadata, source_node_uuid, target_node_uuid, uuid):
	query = "MATCH (c:Compound { struct_id: '" + str(source_node) + "' }) MATCH (t:Target { gene_symbol: '" + target_node + "' }) MERGE (c)-[:DTI "
	query += "{ edge_label: '%s', " % (drug_name + '_' + target_node)
	query += "drug_name: '%s', " % (drug_name) 	
	query += "source_node: '%d', " % (source_node) 
	query += "target_node: '%s', " % (target_node)
	query += "action_type: '%s', " % (action_type)
	query += "p_chembl: '%.3f', " % (p_chembl)
	query += "is_activity_known: '%s', " % (str(is_activity_known))
	query += "priority: '%d', " % (priority)	
	query += "source: '%s', " % (source)	
	query += "relationship: '%s', " % (relationship)
	query += "comment: '%s', " % (comment)
	query += "pmids: '%s', " % (pmids)	
	query += "metadata: '%s', " % (metadata)	
	query += "is_in_drugcentral: '%s', " % (str(is_in_drugcentral))
	query += "is_in_jdti: '%s', " % (str(is_in_jdti))
	query += "source_node_uuid: '%s', " % (source_node_uuid) 
	query += "target_node_uuid: '%s', " % (target_node_uuid)
	query += "uuid: '%s' }]->(t)" % (uuid)

	#print (query)
	
	conn.run(query)	



