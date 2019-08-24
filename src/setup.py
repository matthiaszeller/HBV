#!/usr/bin/env python

# Sets up all constants and some imports.
# Aims to be imported this way: from setup import *

# ############################################################################### #
# =================================== IMPORTS =================================== #
# ############################################################################### #

import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 
import pickle
import os

# ############################################################################### #
# =============================== PATH MANAGEMENT =============================== #
# ############################################################################### #

PATH_PROJECT="/home/mazeller/HBV/"

############ DATA
PATH_DATA = "data/"
### RAW DATA
PATH_RAW_DATA = PATH_DATA+"raw/"
PATH_CLINICAL_RAW_DATA = PATH_RAW_DATA+"to_fellay_manifest.csv"
PATH_VIRAL_RAW_DATA = PATH_RAW_DATA+"viral_seq/OUT_aa_binary_table.txt"
PATH_HOST_RAW_DATA = PATH_RAW_DATA+"wes_plink/"
PATH_HOST_RAW_PLINK_DATA = PATH_HOST_RAW_DATA+"hbv_gilead_bcftools"
### REDUCED DATA
PATH_HOST_RAW_PLINK_REDUCED = PATH_HOST_RAW_DATA+"host_reduced"
### PROCESSED DATA
PATH_VIRAL_DATA = PATH_DATA+"viral_data"
PATH_VIRAL_DATA_PCS = PATH_DATA+"viral_data_pcs"
PATH_CLINICAL_DATA = PATH_DATA+"clinical_data"
PATH_CLINICAL_DATA_TRANSFORMED = PATH_DATA+"clinical_data_transformed"
PATH_VIRAL_GROUPED_DATA = PATH_DATA+"viral_seq_grouped"
### HOST DATA
PATH_PLINK_DATA = PATH_DATA+"plink/"
PATH_SNP_MAP = PATH_DATA+'SNP_dict'
PATH_CLINICAL_PLINK_PHENOTYPE = PATH_PLINK_DATA+'pheno.txt'
PATH_HOST_CLEAN_DATA = PATH_PLINK_DATA+"host_geno_clean"

PATH_PLINK_LOG = PATH_PLINK_DATA+'log/'

### GWAS/G2G DATA
PATH_ASIANS_GWAS = PATH_PLINK_DATA+"asians_gwas"
PATH_ASIANS_GWAS_COVARIATES = PATH_DATA + "covariates_asians_gwas.txt"
PATH_ASIANS_GWAS_PHENOTYPES = PATH_DATA + "phenotypes_asians_gwas.txt"

PATH_GWAS_FINAL_RESULTS = PATH_DATA+"sig_associations.txt"
PATH_GWAS_RESULTS = PATH_PLINK_DATA + "gwas/"
PATH_ASIANS_GWAS_RESULTS = PATH_GWAS_RESULTS + "asians_gwas"
PATH_ASIANS_GWAS_TMP = "/scratch/mazeller/gwas_tmp_results/"
PATH_ASIANS_GWAS_TMP_RESULTS = PATH_ASIANS_GWAS_TMP+"asians_gwas"
PATH_ASIANS_GWAS_TMP_LOGS_CONCAT = PATH_ASIANS_GWAS_TMP+"concatenated_logs"
PATH_ASIANS_GWAS_TMP_WARNINGS = PATH_ASIANS_GWAS_TMP+"concatenated_warnings"
PATH_ASIANS_GWAS_TMP_SUMMARY = PATH_GWAS_RESULTS+"pval_summary"
PATH_ASIANS_GWAS_TMP_SUMMARY_SNPS = PATH_GWAS_RESULTS+"pval_summary_snps"
PATH_ASIANS_GWAS_FILTERED = PATH_GWAS_RESULTS+"filtered_g2g"
PATH_ASIANS_GWAS_RAND = PATH_GWAS_RESULTS+"asians_gwas_random"
PATH_SIGNIFICANT_AAS = PATH_GWAS_RESULTS + "asians_gwas"

#PATH_ASIANS_GWAS_RESULTS_LIST = PATH_GWAS_RESULTS + "output_list"

### INDIVIDUALS
PATH_INTERSECTING_INDIVIDUALS = PATH_DATA+'all_intersecting_ids'
PATH_CLUSTERING_ASIANS = PATH_DATA+"clustering_asians_gwas.txt"

############ SOURCES
PATH_SRC = "src/"
### PLINK SCRIPTS
PATH_SRC_PLINK = PATH_SRC+"plink/"

############ TESTS
PATH_TEST = PATH_DATA+'test/'

############ SPECIAL
PATH_NULL = PATH_DATA+"null"

# ################################################################################ #
# ============================= CURRENT COMPUTATIONS ============================= #
# ################################################################################ #

PATH_WORKING_PHENOTYPES = PATH_DATA+'working_pheno'
PATH_WORKING_COVARIATES = PATH_DATA+'working_covariates'

# ############################################################################### #
# =============================== DATA MANAGEMENT =============================== #
# ############################################################################### #

ID_IGM_CLINICAL_DF = 'IGM_ID'
ID_GS_CLINICAL_DF = 'gilead_id'
ID_GS_VIRAL_DF = 'id'

CLINICAL_COVARIATES = ['AGE', 'SEX']

DEFAULT_FORCE_PLINK_COMPUTATIONS = True

############ VIRAL DATA CLEANING

# MISSING VALUES
INDIVIDUALS_THRESHOLD_MAX_MISSING = 0.1
# VARIANTS FREQUENCY
VARIANTS_THRESHOLD_FREQUENCY = 0.05 # bilateral

############ HOST DATA CLEANING
THRESHOLD_MISSING_INDIVIDUALS = 0.1
THRESHOLD_MISSING_VARIANTS = 0.1	 # 20 august: 0.01->0.1
THRESHOLD_HWE = 1e-6
THRESHOLD_MAF = 0.05

###### MERGING DATA
VIRAL_CLINICAL_VARIABLES = [ID_GS_CLINICAL_DF, 'GT', 'RACE', 'COUNTRY']
HOST_CLINICAL_VARIABLES = [ID_IGM_CLINICAL_DF, 'GT', 'RACE', 'COUNTRY']

###### CHROMOSOMES
DEFAULT_CHROMOSOME_EXCLUSION = '--not-chr 0 6 X Y XY MT'
DEFAULT_CHROMOSOME_EXCLUSION_GWAS = '--not-chr 0 X Y XY MT'

############ DATA ENCODING
# For plink-readable format
PLINK_NA_VALUE_REPRESENTATION = -9
PLINK_PHENOTYPE_NA_REP = "NA" # the default NaN recognition by plink

############ GWAS
GWAS_HOST_NUMBER_PCS = 4
