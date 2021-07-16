# Author:	Gergely Zahoranszky-Kohalmi, PhD
#
# Email:	gergely.zahoranszky-kohalmi@nih.gov
#
# Organization:	National Center for Advancing Translational Sciences (NCATS/NIH)
#
# References
#

import pandas as pd

from sg import process_sg
from n4c_commons import *


# Input

FILE_SG_FORWARD = '../data/input/SG_HATs_dist_3_conf_0.00.json'

FILE_SG_REVERSE = '../data/input/SG_HATs_reverse_dist_3_conf_0.00.json'


# Input/Output
FILE_STD_PPI = '../data/input/standardized/STD_ppi.tsv'

# Output

FILE_STD_SG = '../input/standardized/'


print ('\n\n[*] Processing SmartGraph input started ...')

df_fw = process_sg (FILE_SG_FORWARD)
df_rev = process_sg (FILE_SG_REVERSE)

df = df_fw.append(df_rev, ignore_index = True)

df = df[['source_node', 'target_node', 'interaction', 'mechanism', 'source_specific_score']].copy()




									
# Generic to all datasets
df['is_experimental'] = True
df['data_source'] = 'SmartGraph'
df['abbreviated_data_source'] = 'ppi_sg'
df['acquisition_method'] = 'SmartGraph analysis max. dist=3, min. conf=0'
df['prioritized_for_pathway_analysis'] = False
df['do_ppi_expansion'] = False
df['source_specific_score_type'] = 'confidence'
df['directed'] = True
df['metadata'] = ''



agg_cols = ['source_node', 'target_node', 'interaction', 'mechanism', 'metadata']
df = aggregate_by_first (df, agg_cols)



#agg_cols = ['source_node', 'target_node']
#df = aggregate_by_concat (df, agg_cols)

df = df.rename (columns = {
		'source_node': 'host_protein_a',
		'target_node': 'host_protein_b'
	})




df_ppi = pd.read_csv (FILE_STD_PPI, sep = '\t')

df_ppi = df_ppi.append (df, ignore_index = True)

df_ppi.to_csv (FILE_STD_PPI, sep = '\t', index = False)

print ('\n[Done.]\n')
