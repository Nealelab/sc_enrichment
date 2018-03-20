#!/usr/bin/env python

import subprocess
import os

## Inputs ##
INPUT_GENELIST = os.environ['INPUT_GENELIST']
INPUT_SUMSTAT = os.environ['INPUT_SUMSTAT']
COND_ANNOT = os.environ['COND_ANNOT']
PREFIX = os.environ['PREFIX']
OUTLDSCORE = os.environ['OUTLDSCORE']
OUT = os.environ['OUT']

subprocess.call(['/home/sc_enrichement/sc_enrichement-master/main.py',
                    '--main-annot-file',INPUT_GENELIST,
                    '--summary-stats-files',INPUT_SUMSTAT,
                    '--ldscores-prefix',PREFIX,
                    '--out',OUT,
                    '--condition-annot-file',COND_ANNOT,
                    '--export_ldscore_path',OUTLDSCORE,
                    '--windowsize','10000',
                    '--no_baseline',
                    '--verbose'])
