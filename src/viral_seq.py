
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

###################################

if __name__ == '__main__' :

	lst = ['gene_PC_C_pos_0001_M', 'gene_PC_C_pos_0020_S',
			'gene_PC_C_pos_0020_L']
	print(collapse_variants(lst))
	print()
	print(collapse_variants(lst, pos=False))