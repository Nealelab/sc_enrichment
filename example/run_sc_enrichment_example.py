#!/usr/bin/env python

import subprocess
import os

## Inputs ##
INPUT_GENELIST = os.environ['INPUT_GENELIST']
INPUT_SUMSTAT = os.environ['INPUT_SUMSTAT']
PREFIX = os.environ['PREFIX']
OUT = os.environ['OUT']

subprocess.call(['/home/sc_enrichement/sc_enrichement-master/main.py',
                    '--main-annot',INPUT_GENELIST,
                    '--summary-stats-files',INPUT_SUMSTAT,
                    '--ldscores-prefix',PREFIX,
                    '--out',OUT,
                    '--verbose'])
