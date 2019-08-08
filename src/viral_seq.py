
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
# must import explicitely since still experimental
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
import pandas as pd
import numpy as np

################### VIRAL SEQ MODULE

# Aim: perform computations on viral sequencing data


################### FUNCTIONS

def collapse_variants(lst, pos=True):
	"""Transform lst into a list of position in the viral genome.
	If pos=False, returns the list of genes.
	Ex: lst = ['gene_X_pos_001_M', 'gene_X_pos_001_L', 'gene_S_pos_0030_L']
	Will return ['gene_X_pos_001', 'gene_S_pos_030'], or
	['gene_X', 'gene_S'] if pos=False.

	Input: lst, list of strings"""

	positions = []
	substr = '_'
	# If we're only interested in the gene names,
	# we want to get rid of the substr '_pos_XXX'
	if pos==False: 
		substr = '_pos_'

	for i in lst :
		# Warning: don't use splicing[-1], because of the STOP codon
		idx = i.rfind(substr)
		positions.append(i[0:idx])

	# Remove duplicates
	positions = list(dict.fromkeys(positions))
	return positions


def index_by_pos(df, gene, pos):
	pass


def pca_impute(df, n_components, impute_strategy):
	"""Returns a sklearn PCA object. 
	impute_strategy: mean or most frequent implies a SimpleImputer.
					 iterative implies an iterative imputer
	Note: for the iterative imputation"""
	# Impute missing values
	df_pca = df.transpose()

	if impute_strategy == 'iterative':
		imputer = IterativeImputer(min_value=0, max_value=1)
	else:
		imputer = SimpleImputer(strategy=impute_strategy)

	imputed = imputer.fit_transform(df_pca)

	# Rebuild the DataFrame with the imputed values
	df_pca = pd.DataFrame(data=imputed, columns=df_pca.columns, 
						  index=df_pca.index)

	# PCA
	pca = PCA(n_components = n_components)
	pcs = pca.fit(df_pca.transpose())

	return pcs


###################################

if __name__ == '__main__' :

	# Test collapse_variants
	lst = ['gene_PC_C_pos_0001_M', 'gene_PC_C_pos_0020_S',
			'gene_PC_C_pos_0020_L']
	print(collapse_variants(lst))
	print()
	print(collapse_variants(lst, pos=False))

	# Test pca_impute
	df = pd.DataFrame([[0,1],[np.nan,0],[2,3],[2,np.nan]])
	print()
	print(df)

	pca = pca_impute(df, 2, impute_strategy='mean')
	