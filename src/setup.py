#!/usr/bin/env python

########################################
############## IMPORTS #################
########################################

import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 
import pickle
import os

########################################
########### PATH MANAGEMENT ############
########################################

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
PATH_ASIANS_GWAS = PATH_PLINK_DATA+"asians_gwas"

PATH_PLINK_LOG = PATH_PLINK_DATA+'log/'

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

########################################
######## CURRENT COMPUTATIONS ##########
########################################

PATH_WORKING_PHENOTYPES = PATH_DATA+'working_pheno'
PATH_WORKING_COVARIATES = PATH_DATA+'working_covariates'

########################################
########### DATA MANAGEMENT ############
########################################

ID_IGM_CLINICAL_DF = 'IGM_ID'
ID_GS_CLINICAL_DF = 'gilead_id'
ID_GS_VIRAL_DF = 'id'

CLINICAL_COVARIATES = ['AGE', 'SEX']

############ VIRAL DATA CLEANING

# MISSING VALUES
INDIVIDUALS_THRESHOLD_MAX_MISSING = 0.1
# VARIANTS FREQUENCY
VARIANTS_THRESHOLD_FREQUENCY = 0.05 # bilateral

############ HOST DATA CLEANING
THRESHOLD_MISSING_INDIVIDUALS = 0.1
THRESHOLD_MISSING_VARIANTS = 0.01

###### MERGING DATA
VIRAL_CLINICAL_VARIABLES = [ID_GS_CLINICAL_DF, 'GT', 'RACE', 'COUNTRY']
HOST_CLINICAL_VARIABLES = [ID_IGM_CLINICAL_DF, 'GT', 'RACE', 'COUNTRY']
###### CHROMOSOMES
DEFAULT_CHROMOSOME_EXCLUSION = '--not-chr 0 6 X Y XY MT'
DEFAULT_CHROMOSOME_EXCLUSION_GWAS = '--not-chr 0 X Y XY MT'

############ DATA ENCODING
# For plink-readable format
PLINK_NA_VALUE_REPRESENTATION = -9

############ 
