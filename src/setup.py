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
PATH_CLINICAL_DATA = PATH_DATA+"clinical_data"
PATH_VIRAL_GROUPED_DATA = PATH_DATA+"viral_seq_grouped"
### HOST DATA
PATH_PLINK_DATA = PATH_DATA+"plink/"
PATH_SNP_MAP = PATH_DATA+'SNP_dict'
PATH_CLINICAL_PLINK_PHENOTYPE = PATH_PLINK_DATA+'pheno.txt'
PATH_HOST_CLEAN_DATA = PATH_PLINK_DATA+"host_geno_clean"
PATH_HOST_CLEAN_ASIAN = PATH_PLINK_DATA+"host_asian_no_6_x_y"

PATH_PLINK_LOG = PATH_PLINK_DATA+'log/'

### INDIVIDUALS
PATH_INDIVIDUALS = PATH_DATA+'individuals/'
PATH_INTERSECTING_INDIVIDUALS = PATH_INDIVIDUALS+'intersecting_ids'
PATH_INDIVIDUALS_ASIAN = PATH_INDIVIDUALS+'asians'

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

PATH_WORKING_DATASET = PATH_HOST_CLEAN_ASIAN

########################################
########### DATA MANAGEMENT ############
########################################

ID_IGM_CLINICAL_DF = 'IGM_ID'
ID_GS_CLINICAL_DF = 'gilead_id'
ID_GS_VIRAL_DF = 'id'

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

############ DATA ENCODING
# For plink-readable format
PLINK_NA_VALUE_REPRESENTATION = -9

############ 
