import subprocess
import os.path
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def run_shell_command(cmd_lst) :
    """Returns a tuple: (stdout, stderr).
    Note: you must separate the options from their values.
    Note: the output is returned when all calculations are done
    Example: you want 'plink --bfile <file> ...'. Must provide cmd_list=['plink', '--bfile', '<file>', ...]"""
    process = subprocess.Popen(cmd_lst, stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, universal_newlines=True)
    return process.communicate()

def run_plink(command, file, out, extension, plink2=True, verbose_err=True, force=False):
    """Returns a tuple: (stdout, stderr).
    Note: command input must have the options/values separated by a single space"""
    # Check is result already exist
    stdout, stderr = '', ''
    if extension[0] == '.': # avoid unexpected behaviour
        stderr = 'ERROR: run_plink function, extension must not have a dot at beginning.'
    if force==False and os.path.exists(out+'.'+extension): 
        print("run_plink: command '{}', the file '{}' already exists (force is set to False)."
             .format(command, out+'.'+extension))
        # ABORT !
        return ('run_plink was aborted (not an error)', '')
    if not stderr :
        # Construct the shell commands
        lst = []
        # Determine plink version to use
        lst.append('plink2' if plink2 else 'plink')
        # Input and output
        lst.append('--bfile'); lst.append(file)
        lst.append('--out');   lst.append(out)
        # Must separate the options from their values (see run_shell_command)
        commands = command.split(' ')
        for c in commands :
            lst.append(c)
        # Provide values to the run_shell_command function
        stdout, stderr = run_shell_command(lst)
    # Display errors
    if stderr and verbose_err : print(stderr)
    # Return the outputs
    return stdout, stderr


def join_position(df):
    pass



def plot_plink_pca(path, n_pcs=0, scaled=True, h=3, hue_col=None,
                   path_hue=None):
    """Plot PCA whose computation was done with plink.
    Input: - path to the plink .eigenval and .eigenvec files.
           - n_pcs determines the number of principal components to plot (must be multiple of 2).
                   by default, this value is 0 (i.e. only eigenvalues are plotted)
           - scaled: True makes the plot to have consistent ratio with explained variances by each PC
           - h: height of the plot(s). If scaled==True: width consistent with ratio, otherwise width=height
           - hue_col: column name of a categorical pandas.Series to color the data points.
           - path_hue: path to the pickled DataFrame binary file containing hue_col
    Output: None"""
    ext1, ext2 = '.eigenval', '.eigenvec'

    ######################## ERRORS

    if not os.path.exists(path+ext1) or not os.path.exists(path+ext2):
        raise ValueError("plot_plink_pca: the path {} does not have any .eigenval or .eigenval files.".format(path))
    if hue != None and path_hue == None:
        raise ValueError("hue_col provided, but no path_hue specified to load a DataFrame")
    if hue == None and path_hue != None:
        raise ValueError("path_hue provided, but no hue_col specified to select a column of the DataFrame")

    ######################## EIGENVALUES

    # Load the eigenvalues
    with open(path+ext1, 'r') as file :
        content = file.read().split('\n') # 1 eigenvalue per line
    # Process the eigenvalues
    vals = []
    for v in content[0:-1]: # the last item is empty (the file ends with '\n')
        vals.append(float(v))
    # Make the plot of eigenvalues
    sns.barplot(x=np.arange(1, len(vals)+1), y = vals)
    plt.xlabel('Principal components'); plt.ylabel('Eigenvalue')
    # Compute ratios of explained variance
    vals_sum = np.sum(vals)
    vars_expl = []
    for v in vals :
        vars_expl.append(v / vals_sum)

    ######################## EIGENVECTORS
    if n_pcs == 0: return
    if n_pcs % 2 != 0: 
        raise ValueError("ERROR plot_plink_pca: the number of pcs must be a multiple of 2")
    N_plots = int(n_pcs / 2)

    # Load a dataframe from the .eigenvec file
    # Load only the required PCs (usecols)
    txt_pcs = [ 'PC'+str(i) for i in np.arange(1, n_pcs+1) ]
    df_pca = pd.read_csv(path+ext2, sep='\s+', usecols= ['IID']+txt_pcs )

    ###### LOADING HUE
    if path_hue != None : 
        # in this block, hue_col != None
        with open(path_hue, 'rb') as file:
            df = pickle.load(file)
        df = df[hue_col]

    ###### PLOTTING
    def get_labels(i): # returns a 2-tuple
        return ("{} ({:.3}%)".format(txt_pcs[i],   vars_expl[i]*100), 
                "{} ({:.3}%)".format(txt_pcs[i+1], vars_expl[i+1]*100))

    for i_p in np.arange(N_plots):
        i = 2*i_p # more intuitive and convinient this way
        # Compute width of subplot
        width = vars_expl[i]/vars_expl[i+1] * h if scaled else h
        fig, ax = plt.subplots(1,1, figsize=(width, h))
        # Scatter plot
        sns.scatterplot(x=df_pca[txt_pcs[i]], y=df_pca[txt_pcs[i+1]])
        # Axe labels
        xlabel, ylabel = get_labels(i)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel);


########################################
########################################
########################################

if __name__ == '__main__' :

    # plot_plink_pca
    plot_plink_pca("../data/plink/host_geno_clean_pca", n_pcs=2)
    plt.show()
    print(help(plot_plink_pca))
