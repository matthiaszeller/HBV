
# Functions dedicated to host data analysis.

# ############################################################################### #
# =================================== IMPORTS =================================== #
# ############################################################################### #

import subprocess
import os.path
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import re # RegExps
from matplotlib.colors import hex2color
from assocplots import manhattan
from assocplots.qqplot import *

if __name__ != '__main__':
    from src import setup

# ############################################################################### #
# ================================== FUNCTIONS ================================== #
# ############################################################################### #

# =============================================================================== #
# -------------------------- SHELL COMMANDS MANAGEMENT -------------------------- #
# =============================================================================== #

def run_shell_command(cmd_lst, verbose=False) :
    """Returns a tuple: (stdout, stderr).
    Note: you must separate the options from their values.
    Note: the output is returned when all calculations are done
    Example: you want 'plink --bfile <file> ...'. Must provide cmd_list=['plink', '--bfile', '<file>', ...]"""
    if verbose:
        t = ''
        for c in cmd_lst: t += c+' '
        print("Running '{}'".format(t))
    process = subprocess.Popen(cmd_lst, stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, universal_newlines=True)
    return process.communicate()

# ############################################################################### #

def run_plink(command, file, out, extension, plink2=True, verbose=True, force=setup.DEFAULT_FORCE_PLINK_COMPUTATIONS,
              log_name=None):
    """Returns a tuple: (stdout, stderr). If the output file already exists, does not run plink
    Note: command input must have the options/values separated by a single space"""
    stdout, stderr = '', ''
    if extension[0] == '.': # avoid unexpected behaviour
        stderr = 'ERROR: run_plink function, extension must NOT have a dot at beginning.'
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
        stdout, stderr = run_shell_command(lst, verbose)
    # Display errors
    if stderr and verbose : print(stderr)
    # Return the outputs
    if log_name != None:
        path_log = setup.PATH_PLINK_LOG+log_name+'.log'
        with open(path_log, 'w') as file:
            file.write(stdout)
        if verbose: print("Log written to '{}'".format(path_log))
    return stdout, stderr

# ############################################################################### #

def run_assoc(phenotypes=setup.PATH_WORKING_PHENOTYPES, exclude_chrs=setup.DEFAULT_CHROMOSOME_EXCLUSION_GWAS,
              file=setup.PATH_HOST_CLEAN_DATA, out=setup.PATH_HOST_CLEAN_DATA, 
              binary=False, plink2=True, covariates=None, extra=None, hide_covars=False, verbose=True):
    """Call run_plink function together with --keep, --not-chr, or --linear. This function is used by run_gwas.
    Note that computations are always 'forced' since extension is set to ' '
    Inputs: - all the same as run_plink, except phenotypes.
            - phenotypes: path to the file that contains plink phenotypes (also provided to --keep).
            - binary: boolean, if the phenotype is a binary var, we use --logistic and --1 to make a case/control study.
            - covariates: None, or path to the covariates file
            - extra: string, specify extra commands (user should not use it, this param is for run_gwas function)
    Output: (stdout, stderr)"""
    t = ""
    if exclude_chrs != None:
        t += " " + exclude_chrs + " "
    if binary:
        #t += " --make-pheno " + phenotypes + " '*'"
        #t += " --allow-no-sex "
        
        t += " --1 --logistic "
    else:
        t += " --linear "
    # THIS MUST BE RIGHT AFTER --logistic and --linear
    if hide_covars == True:
        t += " hide-covar"
    t += " --pheno " + phenotypes
    t += " --keep " + phenotypes

    if covariates != None:
        t += " --covar " + covariates
    if extra != None:
        t += " " + extra

    t = t.replace('  ', ' ')
    if t[0] == ' ': t = t[1:]

    return run_plink(file=file, out=out, extension=' ', command=t, plink2=plink2,
                     log_name = "assoc", verbose=verbose)

# ############################################################################### #

# TODO, IMPROVEMENT : 
#   1. manually loop over each amino acid to have a better management of what's happenining
#   2. reduce the data (output files) after each association
#   3. parallelize to reduce data whenever a file is generated

def run_gwas(path_phenotypes=setup.PATH_ASIANS_GWAS_PHENOTYPES, path_covariates=setup.PATH_ASIANS_GWAS_COVARIATES,
             use_pheno=None, path_out=setup.PATH_ASIANS_GWAS_RESULTS, hide_covars=False, verbose=True):
    """Calls run_assoc with default parameters for a GWAS with viral amino acids. 
    Get the paths of the output files (for each pheno) as well as stdout, stderr. 
    The list of file paths are STORED in setup.PATH_ASIANS_GWAS_RESULTS_LIST

    Inputs: - use_pheno: None, or string being the column name provided to --pheno-name
                         Warning: if None, you probably want to run write_phenotypes(phenotype='all') BEFORE RUNNING THIS FUNCTION.
            - hide_covars: True by default -> covariates-specific lines in output are dropped (False makes much bigger files !)

    Parameters passed to run_assoc() function :
        - exclude_chrs = None
        - binary = True
        - path_covariates = setup.PATH_ASIANS_GWAS_COVARIATES by default
        - path_phenotypes = setup.PATH_ASIANS_GWAS_PHENOTYPES by default
        - plink2 = True
        - file = setup.PATH_ASIANS_GWAS
        - path_out = setup.PATH_ASIANS_GWAS_RESULTS

    Output: 3-tuple (stdout, stderr)""" 

    ################################
    # ------- SPECIAL CASE ------- #  The user only want a specific amino acid
    ################################  Or we're in the loop
    extra = None
    # If user want to use only 1 phenotype when many are present in the file
    if use_pheno != None:
        # Create extra command to select that specific phenotype
        extra = "--pheno-name " + use_pheno

        o, e = run_assoc(phenotypes=path_phenotypes, exclude_chrs=None, binary=True, 
                         covariates=path_covariates, plink2=True, extra=extra,
                         file=setup.PATH_ASIANS_GWAS, out=path_out, hide_covars=hide_covars,
                         verbose=verbose)
        return (o, e)

    ################################    Test all amino acids
    # ---- SERIOUS THINGS NOW ---- #    MORE INFORMATION: see notebook, features are exposed
    ################################    (G2G asian individuals)
    # We basically want to loop over each possible amino acid manually.
    # 1. Detect all phenotypes / amino acids that will be tested (look in the phenotype file)
    with open(path_phenotypes, 'r') as file:
        first_line = file.readline()
    phenos = first_line.split('\t')[2:] # exclude FID, IID

    # 2. We loop over each amino acid
    k = 0; N = len(phenos) # keep track of the remaining computations
    print("Preparing loop over {} amino acids".format(N))

    #####################################3
    ###################################
    warnings = []
    messages = []
    ####################################
    #################################### FIND STH TO DO WITH IT (messages and warnings)

    file_suffix = '.glm.logistic'
    # Warning: must set float even if real type is int, in order to allow NaN values to be read
    dtype={'#CHROM':int, 'POS':int, 'REF':str, 'ALT':str, 'A1':str, 'TEST':str, 
            'OBS_CT':int, 'OR':float, 'LOG(OR)_SE':float, 'Z_STAT':float, 'P':float}
    # This DataFrame will store everything
    main_df = pd.DataFrame(columns=dtype)

    step = 5 
    get_file_path = lambda aa: setup.PATH_ASIANS_GWAS_TMP_RESULTS+'.'+aa+file_suffix

    for aa in phenos : # aa can be for example PC_149_A, X_5_K
        k += 1
        if k % step == 1 :
            print('.', end='')

        # 1. Call that function for a single amino acid, store in temp dir
        o, e = run_gwas(use_pheno=aa, path_out=setup.PATH_ASIANS_GWAS_TMP_RESULTS, verbose = False)
        warnings.append("("+aa+") " + e)
        if os.path.exists(get_file_path(aa)):
            # Append the log to the concatenated log file
            with open(setup.PATH_ASIANS_GWAS_TMP_RESULTS+'.log', 'r') as file:
                content = file.read()
            with open(setup.PATH_ASIANS_GWAS_TMP_LOGS_CONCAT, 'a+') as file:
                file.write("\n" + content)
            # 2. Filter the output
            #df = pd.read_csv(get_file_path(aa), sep='\t', dtype=dtype)
            #df = df[ df.P < setup.GWAS_TMP_FILTERING_P_THRESHOLD ]
            # 3. Add the specific AA that was tested
            #df['AA'] = np.repeat(aa, df.shape[0])
            #print(df)
        else:
            print(aa, "does not exist")
    print()
    # write warnings
    concat_warnings = ""
    for w in warnings:
        concat_warnings += "\n" + w
    with open(setup.PATH_ASIANS_GWAS_TMP_WARNINGS, 'w') as file:
        file.write(concat_warnings)


# ############################################################################### #

def extract_ADD(path, output_path):
    df = pd.read_csv(path, sep='\t')
    df = df[ df.TEST == 'ADD' ]
    split = path.split('.')
    file_path = output_path + '.' + split[-3] + '.glm.logistic.ADD'
    with open(file_path, 'w') as file:
        file.write(df.to_csv(index=False, sep='\t'))
    print('.', end='')
    return file_path
    
# =============================================================================== #
# ----------------------------- PLOTTING MANAGEMENT ----------------------------- #
# =============================================================================== #

def plot_plink_pca(path, n_pcs=0, scaled=True, h=3, hue_col=None,
                   path_hue=None, bbox_to_anchor=None):
    """Plot PCA whose computation was done with plink.
    Input: - path to the plink .eigenval and .eigenvec files.
           - n_pcs determines the number of principal components to plot (must be multiple of 2).
                   by default, this value is 0 (i.e. only eigenvalues are plotted)
           - scaled: True makes the plot to have consistent ratio with explained variances by each PC
           - h: height of the plot(s). If scaled==True: width consistent with ratio, otherwise width=height
           - hue_col: column name of a categorical pandas.Series from the clinical dataframe, or a Series (MUST BE IN THE RIGHT ORDER)
           - path_hue: path to the pickled DataFrame binary file containing hue_col
    Output: None"""
    ext1, ext2 = '.eigenval', '.eigenvec'

    ######################## ERRORS

    if not os.path.exists(path+ext1) or not os.path.exists(path+ext2):
        raise ValueError("plot_plink_pca: the path {} does not have any .eigenval or .eigenval files.".format(path))
    #if hue != None and path_hue == None:
    #    raise ValueError("hue_col provided, but no path_hue specified to load a DataFrame")
    #if hue == None and path_hue != None:
    #    raise ValueError("path_hue provided, but no hue_col specified to select a column of the DataFrame")

    ######################## EIGENVALUES

    # Load the eigenvalues
    with open(path+ext1, 'r') as file :
        content = file.read().split('\n') # 1 eigenvalue per line
    # Process the eigenvalues
    vals = []
    for v in content[0:-1]: # the last item is empty (the file ends with '\n')
        vals.append(float(v))
    # Make the plot of eigenvalues
    plt.figure(figsize=(4,2))
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
    hue = None
    if type(hue_col) == list or type(hue_col) == np.ndarray :
        hue = hue_col
    elif hue_col != None:
        # Load pickled dataframe
        with open(setup.PATH_CLINICAL_DATA, 'rb') as file:
            df_clinical = pickle.load(file)
        df_clinical = df_clinical[[setup.ID_IGM_CLINICAL_DF, hue_col]]
        df_clinical.set_index(setup.ID_IGM_CLINICAL_DF, inplace=True)
        df_pca.set_index('IID', inplace=True)
        # Join the tables with the IDs
        df_pca = df_pca.join(other=df_clinical)
        hue = df_pca[hue_col]

    ###### PLOTTING
    def get_labels(i): # returns a 2-tuple
        return ("{} ({:.3}%)".format(txt_pcs[i],   vars_expl[i]*100), 
                "{} ({:.3}%)".format(txt_pcs[i+1], vars_expl[i+1]*100))

    for i_p in np.arange(N_plots):
        i = 2*i_p # more intuitive and convinient this way
        # Compute width of subplot
        width = vars_expl[i]/vars_expl[i+1] * h if scaled else h
        fig, ax = plt.subplots(1,1, figsize=(width, h))
        # Figure out if we must plot the legend for this subplot
        lg = False if i_p < N_plots-1 else 'brief'
        # Scatter plot
        sns.scatterplot(x=df_pca[txt_pcs[i]], y=df_pca[txt_pcs[i+1]], hue=hue, legend=lg)
        # Axe labels
        xlabel, ylabel = get_labels(i)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel);

    # Reset position of the legend
    if bbox_to_anchor is not None and hue is not None:
        plt.gca().legend(bbox_to_anchor=bbox_to_anchor, 
            loc=2, borderaxespad=0.);


def gen_manhattan(path_list, threshold=None):
    ################# COLORS
    chrs = [str(i) for i in range(1,23)]
    chrs_names = np.array([str(i) for i in range(1,23)])
    chrs_names[1::2] = ''
    colors = ['#1b9e77', "#d95f02", '#7570b3', '#e7298a']
    # Converting from HEX into RGB
    colors = [hex2color(colors[i]) for i in range(len(colors))]
    ################# LOAD DAT
    # TODO: manage better paths and amino acid name
    aa_names = [ file.split('.')[-4] for file in path_list ]
    list_p_series = []
    for i, file in enumerate(path_list):
        df = pd.read_csv(file, sep='\t')
        list_p_series.append(df.P)
        ################# PLOT
        manhattan.manhattan(df['P'], df['POS'], df['#CHROM'].astype(str), '',
               plot_type='single',
               chrs_plot=[str(i) for i in range(1,23)],
               chrs_names=chrs_names,
               cut = 0,
               title=aa_names[i],
               xlabel='chromosome',
               ylabel='-log10(p-value)',
               lines= [],
               #lines_colors=['r'],
               colors = colors)
        if threshold != None :
            minlogthreshold = - np.log10(threshold)
            plt.plot([0, 5000000000], [minlogthreshold,minlogthreshold])
        plt.show()

    ## QQ PLOTS
    N = len(aa_names)
    fill_dens = [0.2] * N 
    # ['#1b9e77', "#d95f02", '#7570b3']
    colors = ['#c0392b', '#884ea0', '#2471a3', '#17a589', '#229954', '#d4ac0d', '#ca6f1e',
              '#839192', '#17202a']
    qqplot(list_p_series, 
       aa_names, 
       color=colors[0:N], 
       fill_dens=fill_dens, 
       error_type='theoretical', 
       distribution='beta',
       title='QQ plots')
    plt.gca().legend(bbox_to_anchor=(1.1,0.1*N), 
            loc=2, borderaxespad=0.);
    plt.grid(linestyle='--')


# ############################################################################### #
# ==================================== TESTS ==================================== #
# ############################################################################### #


if __name__ == '__main__' :

    import setup
    # plot_plink_pca
    plot_plink_pca("../data/plink/host_geno_clean", n_pcs=2,
                   hue_col='GT')
    plt.show()
    print(help(plot_plink_pca))
