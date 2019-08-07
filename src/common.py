import os.path
import pickle

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