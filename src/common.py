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



def plot_pca(pcs, pc_plot_height, n_plots=1, 
			 plt_ratio=True, figsize=(10,3)):

	# Define a function that writes axes labels
	def get_labels(k, x_ratio, y_ratio):
		# Set axis labels
		# Note: we use 2*k+1 and 2*k+2 since humans prefer
		# 		to count from 1 rather than 0.
		return ( "PC{} ({:.3}%)".format(2*k+1, x_ratio*100),
		 		 "PC{} ({:.3}%)".format(2*k+2, y_ratio*100))

	################# Core of plot_pca function

	# Plot singular values
	sns.barplot(x=np.arange(1, len(pcs.singular_values_)+1), 
				y=pcs.singular_values_);
	plt.xlabel('principal components')
	plt.ylabel('singular values');

	# Manage the number of plots that must be drawn (parameter)
	fig, axs = plt.subplots(1,n_plots, figsize=figsize)
	if n_plots == 1: axs = [axs]

	# Scatter plot of the data's principal components
	k = 0
	while k < n_plots :	
		sns.scatterplot(x=pcs.components_[2*k],
						y=pcs.components_[2*k+1], ax = axs[k])
		xlabel, ylabel = get_labels(k, 
				   x_ratio=pcs.explained_variance_ratio_[2*k], 
				   y_ratio=pcs.explained_variance_ratio_[2*k+1])
		axs[k].set_xlabel(xlabel)
		axs[k].set_ylabel(ylabel)
		k += 1

	# Make some visual adjustments
	if n_plots > 1: 
		plt.subplots_adjust(wspace=0.4)

	#### Plot with a consistent aspect ratio
	if plt_ratio == False: # but only if user want it
		return

	k = 0
	while k < n_plots :
		# Compute this ratio
		ratio = pcs.explained_variance_ratio_[2*k]/pcs.explained_variance_ratio_[2*k+1]
		# Set figure size
		fig, ax = plt.subplots(1,1,figsize=(ratio*pc_plot_height, pc_plot_height))
		# Plot
		sns.scatterplot(x=pcs.components_[2*k], y=pcs.components_[2*k+1]);
		xlabel, ylabel = get_labels(k, 
				   x_ratio=pcs.explained_variance_ratio_[2*k], 
				   y_ratio=pcs.explained_variance_ratio_[2*k+1])
		ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)

		k += 1


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