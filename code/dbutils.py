# Author:       Gergely Zahoranszky-Kohalmi, PhD
#
# Email:        gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
# References
#
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
# Ref: https://stackoverflow.com/questions/29095795/how-to-import-csv-containing-boolean-values-into-neo4j/42425004


from py2neo import Graph


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


def generate_data_source_fields (abbreviated_data_source, unique_abbr_data_sources):
	dss = []
	subquery = ''
	
	dss = abbreviated_data_source.split ('||')
	
	for ds in unique_abbr_data_sources:
		subquery += '_'.join(['is_in', ds])
		if ds in dss:
			subquery += ': True, '
		else:
			subquery += ': False, '
	
	return (subquery)
	

def create_host_protein_node(conn, host_protein, abbreviated_data_source, unique_abbr_data_sources, is_experimental, data_source, acquisition_method, source_specific_score, source_specific_score_type, activation, activation_type, metadata, uniprot, tdl, uuid):
	query = "CREATE (h:HostProtein {"
	query += "name: '%s', " % (host_protein)
	query += "node_type: 'host_protein', "
	query += "gene_symbol: '%s', " % (host_protein)
	query += "abbreviated_data_source: '%s', " % (abbreviated_data_source)
	query += "data_source: '%s', " % (data_source)
	query += "is_experimental: '%s', " % (is_experimental)
	query += "acquisition_method: '%s', " % (acquisition_method)	
	query += "source_specific_score: '%s', " % (source_specific_score)	
	query += "source_specific_score_type: '%s', " % (source_specific_score_type)	
	query += "activation: '%s', " % (activation)
	query += "activation_type: '%s', " % (activation_type)
	query += "metadata: '%s', " % (metadata)
	query += "uniprot: '%s', " % (uniprot)
	query += "tdl: '%s', " % (tdl)
	query += generate_data_source_fields (abbreviated_data_source, unique_abbr_data_sources)
	query += "uuid: '%s'})" % (uuid)

	
	
	#print (query)
	
	conn.run(query)




def create_pathogen_protein_node(conn, pathogen_protein, abbreviated_data_source, unique_abbr_data_sources, is_experimental, data_source, acquisition_method, source_specific_score, source_specific_score_type, activation, activation_type, metadata, uuid):
	query = "CREATE (p:PathogenProtein {"
	query += "name: '%s', " % (pathogen_protein)
	query += "node_type: 'pathogen_protein', "
	query += "abbreviated_data_source: '%s', " % (abbreviated_data_source)
	query += "data_source: '%s', " % (data_source)
	query += "is_experimental: '%s', " % (is_experimental)
	query += "acquisition_method: '%s', " % (acquisition_method)	
	query += "source_specific_score: '%s', " % (source_specific_score)	
	query += "source_specific_score_type: '%s', " % (source_specific_score_type)	
	query += "activation: '%s', " % (activation)
	query += "activation_type: '%s', " % (activation_type)
	query += "metadata: '%s', " % (metadata)
	query += generate_data_source_fields (abbreviated_data_source, unique_abbr_data_sources)
	query += "uuid: '%s'})" % (uuid)

	
	
	#print (query)
	
	conn.run(query)

def create_host_protein_indices (conn):
	query = 'CREATE CONSTRAINT ON (h:HostProtein) ASSERT h.gene_symbol IS UNIQUE'

	#print (query)
	
	conn.run(query)
	
	
	query = 'CREATE CONSTRAINT ON (h:HostProtein) ASSERT h.name IS UNIQUE'

	#print (query)
	
	conn.run(query)

	query = 'CREATE CONSTRAINT ON (h:HostProtein) ASSERT h.uuid IS UNIQUE'

	#print (query)
	
	conn.run(query)
	
	
def create_pathogen_protein_indices (conn):
	query = 'CREATE CONSTRAINT ON (p:PathogenProtein) ASSERT p.name IS UNIQUE'

	#print (query)
	
	conn.run(query)

	query = 'CREATE CONSTRAINT ON (p:PathogenProtein) ASSERT p.uuid IS UNIQUE'

	#print (query)
	
	conn.run(query)	





def create_drug_node(conn, drug_name, abbreviated_data_source, unique_abbr_data_sources, data_source, acquisition_method, is_experimental, smiles, inchi, inchi_key, ns_inchi_key, CAS_RN, metadata, uuid):
	query = "CREATE (c:Drug {"
	query += "drug_name: '%s', " % (drug_name)
	query += "name: '%s', " % (drug_name)
	query += "node_type: 'drug', "
	query += "abbreviated_data_source: '%s', " % (abbreviated_data_source)
	query += "data_source: '%s', " % (data_source)
	query += "is_experimental: '%s', " % (is_experimental)
	query += "acquisition_method: '%s', " % (acquisition_method)	
	query += "smiles: '%s', " % (smiles.replace ('\\', '\\\\'))
	query += "inchi: '%s', " % (inchi.replace ('\\', '\\\\'))
	query += "inchi_key: '%s', " % (inchi_key)
	query += "ns_inchi_key: '%s', " % (ns_inchi_key)
	query += "CAS_RN: '%s', " % (CAS_RN)
	query += "metadata: '%s', " % (metadata)
	query += generate_data_source_fields (abbreviated_data_source, unique_abbr_data_sources)	
	query += "uuid: '%s'})" % (uuid)

	
	
	#print (query)
	
	conn.run(query)




def create_compound_indices (conn):
	query = 'CREATE CONSTRAINT ON (d:Drug) ASSERT d.drug_name IS UNIQUE'

	#print (query)
	
	conn.run(query)

	query = 'CREATE CONSTRAINT ON (d:Drug) ASSERT d.uuid IS UNIQUE'

	#print (query)
	
	conn.run(query)





def create_ppi_edge (conn, source_node, target_node, interaction, mechanism, source_specific_score, is_experimental, data_source, abbreviated_data_source, unique_abbr_data_sources, acquisition_method, source_specific_score_type, directed, edge_label, ref_annotation, ref_direction, ref_score, ref_interaction, metadata, source_node_uuid, target_node_uuid, uuid):
	query = "MATCH (t1:HostProtein { gene_symbol: '" + source_node + "' }) MATCH (t2:HostProtein { gene_symbol: '" + target_node + "' }) MERGE (t1)-[:PPI "
	query += "{ edge_label: '%s', " % (edge_label)
	query += "source_node: '%s', " % (source_node) 
	query += "target_node: '%s', " % (target_node)
	query += "edge_type: 'ppi', "
	query += "interaction: '%s', " % (interaction)
	query += "mechanism: '%s', " % (mechanism)
	query += "abbreviated_data_source: '%s', " % (abbreviated_data_source)
	query += "data_source: '%s', " % (data_source)
	query += "is_experimental: '%s', " % (is_experimental)
	query += "acquisition_method: '%s', " % (acquisition_method)	
	query += "source_specific_score: '%s', " % (source_specific_score)	
	query += "source_specific_score_type: '%s', " % (source_specific_score_type)	
	query += "directed: '%s', " % (directed)	
	query += "ref_annotation: '%s', " % (ref_annotation)	
	query += "ref_direction: '%s', " % (ref_direction)		
	query += "ref_score: '%s', " % (ref_score)	
	query += "ref_interaction: '%s', " % (ref_interaction)
	query += "metadata: '%s', " % (metadata)	
	query += generate_data_source_fields (abbreviated_data_source, unique_abbr_data_sources)
	query += "source_node_uuid: '%s', " % (source_node_uuid)
	query += "target_node_uuid: '%s', " % (target_node_uuid)	
	query += "uuid: '%s' }]->(t2)" % (uuid)

			

	#print (query)
	
	conn.run(query)
	
	


def create_hpi_edge (conn, source_node, target_node, interaction, mechanism, source_specific_score, is_experimental, data_source, abbreviated_data_source, unique_abbr_data_sources, acquisition_method, source_specific_score_type, directed, edge_label, metadata, source_node_uuid, target_node_uuid, uuid):
	query = "MATCH (p:PathogenProtein { name: '" + source_node + "' }) MATCH (t:HostProtein { gene_symbol: '" + target_node + "' }) MERGE (p)-[:HPI "
	query += "{ edge_label: '%s', " % (edge_label)
	query += "source_node: '%s', " % (source_node) 
	query += "target_node: '%s', " % (target_node)
	query += "edge_type: 'hpi', "
	query += "interaction: '%s', " % (interaction)
	query += "mechanism: '%s', " % (mechanism)
	query += "abbreviated_data_source: '%s', " % (abbreviated_data_source)
	query += "data_source: '%s', " % (data_source)
	query += "is_experimental: '%s', " % (is_experimental)
	query += "acquisition_method: '%s', " % (acquisition_method)	
	query += "source_specific_score: '%s', " % (source_specific_score)	
	query += "source_specific_score_type: '%s', " % (source_specific_score_type)	
	query += "directed: '%s', " % (directed)	
	query += "metadata: '%s', " % (metadata)	
	query += generate_data_source_fields (abbreviated_data_source, unique_abbr_data_sources)
	query += "source_node_uuid: '%s', " % (source_node_uuid)
	query += "target_node_uuid: '%s', " % (target_node_uuid)	
	query += "uuid: '%s' }]->(t)" % (uuid)

			

	#print (query)
	
	conn.run(query)	

	

def create_dti_edge (conn, drug_name, target_node, action_type, p_chembl, is_experimental, data_source, abbreviated_data_source, unique_abbr_data_sources, acquisition_method, source_specific_score_type, directed, source_specific_score, edge_label,metadata, source_node_uuid, target_node_uuid, uuid):
	query = "MATCH (d:Drug { drug_name: '" + drug_name + "' }) MATCH (t:HostProtein { gene_symbol: '" + target_node + "' }) MERGE (d)-[:DTI "
	query += "{ edge_label: '%s', " % (edge_label)
	query += "source_node: '%s', " % (drug_name) 
	query += "edge_type: 'dti', "
	query += "target_node: '%s', " % (target_node)
	query += "action_type: '%s', " % (action_type)
	query += "p_chembl: '%s', " % (p_chembl)
	query += "is_experimental: '%s', " % (is_experimental)
	query += "data_source: '%s', " % (data_source)
	query += "abbreviated_data_source: '%s', " % (abbreviated_data_source)
	query += "acquisition_method: '%s', " % (acquisition_method)	
	query += "source_specific_score: '%s', " % (source_specific_score)	
	query += "source_specific_score_type: '%s', " % (source_specific_score_type)	
	query += "directed: '%s', " % (directed)	
	query += "metadata: '%s', " % (metadata)	
	query += generate_data_source_fields (abbreviated_data_source, unique_abbr_data_sources)
	query += "source_node_uuid: '%s', " % (source_node_uuid)
	query += "target_node_uuid: '%s', " % (target_node_uuid)
	query += "uuid: '%s' }]->(t)" % (uuid)

	#print (query)
	
	conn.run(query)	



