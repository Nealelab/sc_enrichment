#!/usr/bin/env python

import subprocess

## Inputs ##
INPUT_GENELIST = 'gs://singlecellldscore/example/snp_example_anno_noannot'
INPUT_SUMSTAT = 'gs://singlecellldscore/pysch_sumstats/scz_summary_stats.sumstats.gz'
PREFIX = 'kt_test'
OUT = 'gs://singlecellldscore/example/kt_test'
OUT_LDSCORES='gs://singlecellldscore/example/kt_test'


subprocess.call(['/home/sc_enrichement/sc_enrichment-master/main_ldscore.py',
                    '--main-annot-rsids',INPUT_GENELIST,
                    '--summary-stats-files',INPUT_SUMSTAT,
                    '--ldscores-prefix',PREFIX,
		    '--export-ldscore-path',OUT_LDSCORES,
                    '--out',OUT])
