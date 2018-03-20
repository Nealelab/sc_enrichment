#!/usr/bin/env python

import subprocess
import os

## Inputs ##
INPUT_GENELIST = os.environ['INPUT_GENELIST']
INPUT_SUMSTAT = os.environ['INPUT_SUMSTAT']
COND_ANNOT_FILE = os.environ['COND_ANNOT_FILE']
COND_ANNOT_LDSCORE = os.environ['COND_ANNOT_LDSCORE']
PREFIX = os.environ['PREFIX']
OUTLDSCORE = os.environ['OUTLDSCORE']
OUT = os.environ['OUT']

subprocess.call(['/home/sc_enrichement/sc_enrichement-master/main.py',
                    '--main-annot-file',INPUT_GENELIST,
                    '--summary-stats-files',INPUT_SUMSTAT,
                    '--ldscores-prefix',PREFIX,
                    '--out',OUT,
                    '--condition-annot-file',COND_ANNOT_FILE,
                    '--condition-annot-ldscores',COND_ANNOT_LDSCORE,
                    '--export_ldscore_path',OUTLDSCORE,
                    '--windowsize','10000',
                    '--no_baseline',
                    '--verbose'])
