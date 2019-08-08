import os.path
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

####### MODULE COMMON

# Aim: put together functions that can be useful in several contexts


####### FUNCTIONS

def manage_pickle(path, fun, args=None, verbose=True) :
	"""Allows to store results of hard-to-compute processes. Note: this can only manage computations
	that are directly returned by a function (input fun).
	Uses pickle to write and read binary files.
	If the file specified by 'path' exists, returns the content of the file.
	Otherwise, the result of the function fun is returned AND stored in a binary file.

	WARNING: keep in mind that the only way to redo computations is to delete the binary file."""

	# Test if the file already exists
	if os.path.exists(path) :
		# The file exists: load with pickle, and return
		if verbose==True:
			print("manage_pickle: the file '{}' has been loaded"
				.format(path))
		with open(path, 'rb') as file:
			return pickle.load(file)
	else:
		if args==None: result = fun()
		else: result = fun(args)
		with open(path, 'wb') as file:
			pickler = pickle.Pickler(file)
			pickler.dump(result)
		return result



def plot_pca(pcs, pc_plot_height):
	sns.scatterplot(x=pcs.components_[0],
					y=pcs.components_[1])

	# Plot singular values
	sns.barplot(x=np.arange(1, len(pcs.singular_values_)+1), 
				y=pcs.singular_values_);
	plt.xlabel('principal components')
	plt.ylabel('singular values');
	# Plot data w.r.t. the PCs
	ratio = pcs.explained_variance_ratio_[0]/pcs.explained_variance_ratio_[1]
	plt.subplots(1,1,figsize=(ratio*pc_plot_height, pc_plot_height))

	sns.scatterplot(x=pcs.components_[0], y=pcs.components_[1]);
	plt.xlabel("PC1 ({:.3}%)".
		format(pcs.explained_variance_ratio_[0]*100))
	plt.ylabel("PC2 ({:.3}%)".
		format(pcs.explained_variance_ratio_[1]*100));

##############################


if __name__ == '__main__' :

	n = 4
	path = "test_manage_pickle"

	def foo(x):
		return x**2;

	result_of_computation = manage_pickle(path, foo, n)

	def foo2():
		return 0
	result_of_computation = manage_pickle('test_manage_pickle2',
										foo2)