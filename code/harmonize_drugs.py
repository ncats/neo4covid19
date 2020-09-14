# Author:       Gergely Zahoranszky-Kohalmi, PhD
#
# Email:        gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advanving Translational Sciences (NCATS/NIH)
#
# References
#

import pandas as pd


#FILE_jdti = '../out_jdti.tab'
#FILE_full_jdti = '../out_jdti_all.tab'
#FILE_dcdti = '../out_dcdti.tab'


#df_jdti = pd.read_csv (FILE_jdti, sep = '\t')
#df_full_jdti = pd.read_csv (FILE_full_jdti, sep = '\t')

#df_dcdti = pd.read_csv (FILE_dcdti, sep = '\t')


#print (df_jdti.columns)
#print (df_full_jdti.columns)
#print (df_dcdti.columns)



def assign_missing_struct_id (df):

	id_map = {}

	names = list(df['drug_name'])
	ids = list(df['struct_id'])



	next_idx = -1

	for i in range(len(ids)):
		d_id = ids[i]
		name = names[i]

		if d_id == -1:
			if name not in id_map.keys():
				id_map [name] = next_idx
				next_idx -= 1
		else:
			id_map[name] = d_id

	ser_idx = []
	ser_name = []

	for k in id_map.keys():
		ser_name.append(k)
		ser_idx.append(id_map[k])

	df = df.drop(columns = ['struct_id']).copy()

	df_map = pd.DataFrame ({'drug_name': ser_name, 'struct_id': ser_idx})

	df = df.merge (df_map, left_on = 'drug_name', right_on = 'drug_name', how = 'inner')




	return(df)

def create_edge_label (sn, tn):
	return(sn + '_' + tn)


def harmonize_dtis (df_jdti, df_dcdti):
	### Concatenate JDTIs


	#df_full_jdti = df_full_jdti.rename (columns = {
	#		'Known binders': 'drug_name'
	#	})

	df_jdti = df_jdti.rename (columns = {
			'Known binders': 'drug_name'
		})


	#df_full_jdti = df_full_jdti[~df_full_jdti['drug_name'].isin(list(df_jdti['drug_name']))].copy()




	#df_full_jdti = df_full_jdti[['struct_id', 'human_protein', 'Chemical space', 'source', 'data_origin', 'priority', 'drug_name']].copy()
	#print (df_full_jdti)

	#df_jdti = df_jdti.append(df_full_jdti, ignore_index = True)
	#print (df_jdti)


	### Harmonizing data structure

	df_jdti['p_chembl'] = 3
	df_jdti['action_type'] = 'unknown'
	df_jdti['is_activity_known'] = False
	df_jdti = df_jdti[['struct_id', 'human_protein', 'p_chembl', 'action_type', 'is_activity_known', 'source', 'data_origin', 'priority', 'drug_name']].copy()

	df_dcdti['chemical_space'] = 'approved drugs'



	df_all_dti = df_dcdti.append(df_jdti, ignore_index = True)



	df_all_dti = df_all_dti.sort_values(['priority'])

	df_all_dti = assign_missing_struct_id(df_all_dti)




	df_all_dti['edge_label'] = df_all_dti.apply(lambda x: create_edge_label (str(x['struct_id']), x['human_protein']), axis = 1)


	df_all_dti['relationship'] = 'DTI' 
	
		
	df_all_dti_aggr = df_all_dti.groupby(['struct_id', 'human_protein'], as_index = False).aggregate({
			'drug_name': 'first',
			'p_chembl': 'first',
			'action_type': 'first',
			'is_activity_known': 'first',
			'source': 'first',
			'data_origin': 'first',
			'priority': 'first',
			'relationship': 'first',
			'edge_label': 'first'
		})

	#print (df_all_dti_aggr)

	return (df_all_dti, df_all_dti_aggr)


