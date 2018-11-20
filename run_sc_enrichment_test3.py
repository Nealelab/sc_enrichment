#!/usr/bin/env python

import subprocess

## Inputs ##
MAIN_ANNOT_LDSCORES = 'gs://andrea-ukbb/konrad-pli/ldscore_files100'
INPUT_SUMSTAT_LDSCORE = 'gs://singlecellldscore/PASS/PASS_Schizophrenia.sumstats'
PREFIX = 'ag_test'
COND_LDSCORES='gs://andrea-ukbb/konrad-pli/control_ldscores_geneset/'  
OUT = 'gs://singlecellldscorore/test/' 

subprocess.call(['/home/sc_enrichment/sc_enrichment-master/main_ldscore.py',
				'--summary-stats-files',INPUT_SUMSTAT_LDSCORE,
                    '--prefix',PREFIX,
                    '--main-annot-ldscores-ldcts',MAIN_ANNOT_LDSCORES,
                    '--condition-annot-ldscores',COND_LDSCORES,
                    '--full-report',
                    '--out',OUT,
                    '--verbose'])
