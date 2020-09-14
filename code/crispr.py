# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
# Script: crispr.py
#
# Aim: Parse SARS-CoV-2 related genes implicated by CRIPSR study.
#
# References: Wei et al., bioRxiv, 2020, https://doi.org/10.1101/2020.06.16.155101
#
#

import pandas as pd

FILE_crispr = '../data/input/SARS-CoV-2-CRISPR-paper_mod.xlsx'

def load_crispr ():
	df = pd.read_excel (FILE_crispr, sheet_name = 'Data')
	
	return (df)

