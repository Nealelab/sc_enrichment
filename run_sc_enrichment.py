#!/usr/bin/env python

import subprocess

## Inputs ##
INPUT_GENELIST = os.environ['INPUT_GENELIST']
INPUT_SUMSTAT = os.environ['INPUT_SUMSTAT']
PREFIX = os.environ['PREFIX']
OUT = os.environ['OUT']


subprocess.call(['/home/main.py',
                    '--main-annot-file',INPUT_GENELIST,
                    '--summary-stats-files',INPUT_SUMSTAT,
                    '--ldscores-prefix',PREFIX,
                    '--out',OUT])




