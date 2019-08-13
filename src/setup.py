#!/usr/bin/env python

###### ###### ###### ###### IMPORTS

import pandas as pd

###### ###### ###### ###### PATH MANAGEMENT

PATH_PROJECT="/home/mazeller/HBV/"

### ### DATA
PATH_DATA = "data/"

### RAW DATA
PATH_RAW_DATA = PATH_DATA+"raw/"
PATH_CLINICAL_RAW_DATA = PATH_RAW_DATA+"to_fellay_manifest.csv"
PATH_VIRAL_RAW_DATA = PATH_RAW_DATA+"viral_seq/OUT_aa_binary_table.txt"

PATH_HOST_RAW_DATA = PATH_RAW_DATA+"wes_plink/"
PATH_HOST_RAW_PLINK_DATA = PATH_HOST_RAW_DATA+"hbv_gilead"

### REDUCED DATA
PATH_HOST_RAW_PLINK_REDUCED = PATH_HOST_RAW_DATA+"host_reduced"

### PROCESSED DATA
PATH_VIRAL_DATA = PATH_DATA+"viral_data"
PATH_CLINICAL_DATA = PATH_DATA+"clinical_data"
PATH_VIRAL_GROUPED_DATA = PATH_DATA+"viral_seq_grouped"

# Host data
PATH_PLINK_DATA = PATH_DATA+"plink/"
PATH_SNP_MAP = PATH_DATA+'SNP_dict'
PATH_CLINICAL_PLINK_PHENOTYPE = PATH_PLINK_DATA+'pheno.txt'
PATH_HOST_CLEAN_DATA = PATH_PLINK_DATA+"host_geno_clean"

### ### SOURCES
PATH_SRC = "src/"

### PLINK SCRIPTS
PATH_SRC_PLINK = PATH_SRC+"plink/"

### ### SPECIAL
PATH_NULL = PATH_DATA+"null"

### ### TESTS
PATH_TEST = PATH_DATA+'test/'

###### ###### ###### ###### DATA MANAGEMENT

###### DATA CLEANING

# MISSING VALUES
INDIVIDUALS_THRESHOLD_MAX_MISSING = 0.1

# VARIANTS FREQUENCY
VARIANTS_THRESHOLD_FREQUENCY = 0.05 # bilateral

###### MERGING DATA
VIRAL_CLINICAL_VARIABLES = ['gilead_id', 'GT', 'RACE', 'COUNTRY']


###### DATA ENCODING

# For plink-readable format
PLINK_NA_VALUE_REPRESENTATION = -9