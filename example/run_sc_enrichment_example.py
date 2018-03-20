#!/usr/bin/env python

import subprocess
import os

## Inputs ##
INPUT_GENELIST = os.environ['INPUT_GENELIST']
INPUT_SUMSTAT = os.environ['INPUT_SUMSTAT']
PREFIX = os.environ['PREFIX']
OUT = os.environ['OUT']

subprocess.call(['/home/sc_enrichement/main.py',
                    '--main-annot-file',INPUT_GENELIST,
                    '--summary-stats-files',INPUT_SUMSTAT,
                    '--ldscores-prefix',PREFIX,
                    '--out',OUT,
                    '--verbose'])


# INPUT_GENELIST = "gs://ldscores/example/example.geneset"
# INPUT_SUMSTAT = "gs://ldscores/example/asd_summary_stats.sumstats.gz,gs://ldscores/example/scz_summary_stats.sumstats.gz"
# PREFIX = "example"
# OUT = "gs://ldscores/example/"