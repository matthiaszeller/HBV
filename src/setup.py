

################### PATH MANAGEMENT

# Don't forget path are relative to the root directory

PATH_VIRAL_DATA = "data/viral_data"
PATH_CLINICAL_DATA = "data/clinical_data"



################### DATA MANAGEMENT

###### DATA CLEANING

# MISSING VALUES
INDIVIDUALS_THRESHOLD_MAX_MISSING = 0.1

# VARIANTS FREQUENCY
VARIANTS_THRESHOLD_FREQUENCY = 0.05 # bilateral

###### MERGING DATA
VIRAL_CLINICAL_VARIABLES = ['gilead_id', 'GT', 'RACE', 'COUNTRY']