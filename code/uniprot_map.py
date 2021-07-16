# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advancing Translational Sciences (NCATS/NIH)
#
# References
#

import pandas as pd

url_uniprot_mapping = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'


df_map =  pd.read_csv (url_uniprot_mapping, compression = 'gzip', sep = '\t', quotechar = '"', header = None)
df_map.columns = ['uniprot', 'id_type', 'id_val']

df_map = df_map[df_map['id_type'] == 'Gene_Name'].copy()
df_map = df_map.drop (columns = ['id_type']).copy()
df_map = df_map.rename (columns = {
	'id_val': 'gene_symbol'
})





def gene2uniprot(genes):
	df = pd.DataFrame ({'gene_symbol': genes})

	df = df.merge (df_map, left_on = 'gene_symbol', right_on = 'gene_symbol', how = 'inner')
	

	


	df = df.groupby (['gene_symbol', 'uniprot'], as_index = False).agg ('first')
	


	return (df)
	
	

def uniprot2gene (proteins):
        df = pd.DataFrame ({'uniprot': proteins})

        df = df.merge (df_map, left_on = 'uniprot', right_on = 'uniprot', how = 'inner')
        
        df = df.rename (columns = {
                'genes': 'gene_symbol'
        })

        df = df.groupby (['gene_symbol', 'uniprot'], as_index = False).agg ('first')



        return (df)


