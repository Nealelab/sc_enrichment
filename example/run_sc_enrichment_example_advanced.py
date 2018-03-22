#!/usr/bin/env python

import subprocess
import os

## Inputs ##
INPUT_GENELIST = os.environ['INPUT_GENELIST'] 
INPUT_SUMSTAT = os.environ['INPUT_SUMSTAT']
COND_ANNOT_FILE = os.environ['COND_ANNOT_FILE'] # List of files for conditioning
COND_ANNOT_LDSCORE = os.environ['COND_ANNOT_LDSCORE'] # List of alredy-computed ldscores for conditioning
PREFIX = os.environ['PREFIX']
OUTLDSCORE = os.environ['OUTLDSCORE'] # Where to save the ldscore generated from INPUT_GENELIST
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


# INPUT_GENELIST="gs://ldscores/example/example.geneset,gs://ldscores/example/example2.geneset"
# INPUT_SUMSTAT="gs://ldscores/example/asd_summary_stats.sumstats.gz,gs://ldscores/example/scz_summary_stats.sumstats.gz"
# PREFIX="example,example2"
# OUT="gs://ldscores/example/"
# COND_ANNOT_FILE="gs://ldscores/example/example_cond.geneset,gs://ldscores/example/example_cond2.geneset"
# COND_ANNOT_LDSCORE="gs://ldscores/example/anno_ldscore/"
# OUTLDSCORE="gs://ldscores/example/outld/"