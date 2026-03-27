"""Configuration for TargetSNAP application."""

import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Paths to TargetScan scripts and data
TARGETSCAN_SCRIPT = os.path.join(BASE_DIR, "scripts", "targetscan_70.pl")
TARGETSCAN_BL_PCT_SCRIPT = os.path.join(
    BASE_DIR, "scripts", "targetscan_70_BL_PCT.pl"
)
MIRNA_FAMILY_FILE = os.path.join(BASE_DIR, "data", "miR_Family_Info.txt")

# Temporary directory for TargetScan input/output
TEMP_DIR = os.path.join(BASE_DIR, "tmp")
os.makedirs(TEMP_DIR, exist_ok=True)

# UCSC DAS server for hg19 UTR sequences
UCSC_DAS_URL = "https://genome.ucsc.edu/cgi-bin/das/hg19/dna"

# Ensembl GRCh37 REST API for SNP data
ENSEMBL_REST_URL = "https://grch37.rest.ensembl.org"
