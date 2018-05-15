#!/usr/bin/env python

import subprocess
import os

## Inputs ##
INPUT_GENELIST = os.environ['INPUT_GENELIST'] 
INPUT_SUMSTAT = os.environ['INPUT_SUMSTAT']
COND_ANNOT_FILE = os.environ['COND_ANNOT_FILE'] # List of files for conditioning
PREFIX = os.environ['PREFIX']
OUTLDSCORE = os.environ['OUTLDSCORE'] # Where to save the ldscore generated from INPUT_GENELIST
OUT = os.environ['OUT']

subprocess.call(['/home/sc_enrichement/sc_enrichement-master/main_ldscore.py',
                    '--main-annot',INPUT_GENELIST,
                    '--summary-stats-files',INPUT_SUMSTAT,
                    '--ldscores-prefix',PREFIX,
                    '--out',OUT,
                    '--condition-annot',COND_ANNOT_FILE,
                    '--export_ldscore_path',OUTLDSCORE,
                    '--windowsize','10000',
                    '--verbose'])
