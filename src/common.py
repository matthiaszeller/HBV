
####### MODULE COMMON

# Aim: put together functions that can be useful in several contexts


######## IMPORTS
import os.path
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ != '__main__' :
	from src import setup

####### FUNCTIONS

def write_binary(path, obj):
	with open(path, 'wb') as file:
		pickler = pickle.Pickler(file)
		pickler.dump(obj)

	print("write_binary(): successfully written in '{}'".format(path))

def manage_pickle(path, fun, args=None, verbose=True) :
	"""Store results of hard-to-compute processes. 

	Note: this can only manage computations that are directly returned by a function (input fun).
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
			 plt_ratio=True, figsize=(10,3), 
			 hue=None, scaled_only=False, bbox_to_anchor=None,
			 singular_values=True):
	"""TO DO: implement bbox_to_anchor for n_plots>1 
	if scaled_only=False"""

	# Define a function that writes axes labels
	def get_labels(k, x_ratio, y_ratio):
		# Set axis labels
		# Note: we use 2*k+1 and 2*k+2 since humans prefer
		# 		to count from 1 rather than 0.
		return ( "PC{} ({:.3}%)".format(2*k+1, x_ratio*100),
		 		 "PC{} ({:.3}%)".format(2*k+2, y_ratio*100))

	################# Core of plot_pca function
	axes_labels_fontsize = 14

	# Plot singular values
	if singular_values == True:
		sns.barplot(x=np.arange(1, len(pcs.singular_values_)+1), 
					y=pcs.singular_values_);
		plt.xlabel('principal components')
		plt.ylabel('singular values');

	if scaled_only == False:
		# NOT ON SCALEs
		# Manage the number of plots that must be drawn (parameter)
		fig, axs = plt.subplots(1,n_plots, figsize=figsize)
		if n_plots == 1: axs = [axs]

		# Scatter plot of the data's principal components
		k = 0
		while k < n_plots :	
			#axs[k].scatter(x=pcs.components_[2*k],
			#				y=pcs.components_[2*k+1])
			sns.scatterplot(x=pcs.components_[2*k],
							y=pcs.components_[2*k+1], 
							ax=axs[k], hue=hue)
			xlabel, ylabel = get_labels(k, 
					   x_ratio=pcs.explained_variance_ratio_[2*k], 
					   y_ratio=pcs.explained_variance_ratio_[2*k+1])
			axs[k].set_xlabel(xlabel, fontsize=axes_labels_fontsize);
			axs[k].set_ylabel(ylabel, fontsize=axes_labels_fontsize);

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
		#ax.scatter(x=pcs.components_[2*k], y=pcs.components_[2*k+1]);
		# Put a legend only if it's the last plot
		lg = False if k<n_plots-1 else 'brief'
		sns.scatterplot(x=pcs.components_[2*k], y=pcs.components_[2*k+1],
						ax=ax, hue=hue, legend=lg)
		xlabel, ylabel = get_labels(k, 
				   x_ratio=pcs.explained_variance_ratio_[2*k], 
				   y_ratio=pcs.explained_variance_ratio_[2*k+1])
		ax.set_xlabel(xlabel, fontsize=axes_labels_fontsize); 
		ax.set_ylabel(ylabel, fontsize=axes_labels_fontsize);
		k += 1

	# Reset position of the legend
	if bbox_to_anchor is not None and hue is not None:
		plt.gca().legend(bbox_to_anchor=bbox_to_anchor, 
			loc=2, borderaxespad=0.);


def convert_to_plink_phenotype(path, out, id_col_name, cols) :
	"""Convert a pickled DataFrame into a plink-readable format. 
	Inputs: - path to the pickled dataframe (binary file written with pickle package)
			- out: output file path
			- id_col_name: the column in the DataFrame that contains the relevant ID
			- cols: columns of the dataframe to include in the phenotype file
	Output: None, but write a text file in the specified output path.
	Format: the file begins with two columns (identical in this study) FID and IID, with extra columns depending on the 'cols' provided."""
	# Check required values
	if not os.path.exists(path):
		print("ERROR: convert_to_plink_phenotype, file path does not exist")

	with open(path, 'rb') as file:
		df_pickle = pickle.load(file)

	df = pd.DataFrame()

	for col in cols:
		df.insert(column=col, value=df_pickle[col], loc=0)
	# Create the two id columns (FID and IID)
	df.insert(column='IID', value=df_pickle[id_col_name], loc=0)
	df.insert(column='FID', value=df_pickle[id_col_name], loc=0)
	# Convert to text
	txt = df.to_csv(sep='\t', index=False, 
							 na_rep=setup.PLINK_NA_VALUE_REPRESENTATION)

	with open(out, 'w') as file:
		file.write(txt)
	print("Successfully written '{}'".format(out))



##############################


if __name__ == '__main__' :

	###################################### MANAGE PICLLE
	n = 4
	path = "test_manage_pickle"

	def foo(x):
		return x**2;

	result_of_computation = manage_pickle(path, foo, n)

	def foo2():
		return 0
	result_of_computation = manage_pickle('test_manage_pickle2',
										foo2)

	###################################### CONVERT TO PLINK FORMAT
	import pandas as pd
	import setup

	# Create the dataframe
	def test_convert_to_plink_phenotype():
		df = pd.DataFrame()
		rd = np.random.rand(4)
		rd[2] = np.nan
		df.insert(loc=0, column='var2', value=rd)
		df.insert(loc=0, column='var1', value=np.random.rand(4))
		df.insert(loc=0, column='id', value=['word_id2', 'word_id10', 'word_id000', 'word_id3'])
		return df
	# Pickle the file
	os.system('mkdir ../data/test')
	path_test = '../'+setup.PATH_TEST
	df = manage_pickle(path_test+'convert_to_plink_phenotype', fun=test_convert_to_plink_phenotype)

	# Use the function
	print(df)
	print()
	convert_to_plink_phenotype(path=path_test+'convert_to_plink_phenotype', out=path_test+'converted',
							   id_col_name='id', cols=['var2', 'var1'])