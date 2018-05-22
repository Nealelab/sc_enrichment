#!/usr/bin/env python

import subprocess

## Inputs ##
INPUT_GENELIST = 'gs://singlecellldscore/example/outld/example2/'
INPUT_SUMSTAT = 'gs://singlecellldscore/pysch_sumstats/scz_summary_stats.sumstats.gz'
PREFIX = 'kt_test'
OUT = 'gs://singlecellldscore/example/kt_test'
OUT_LDSCORES='gs://singlecellldscore/example/kt_test'
COND_GENESET='gs://singlecellldscore/control_ldscores_geneset/'
COORD_FILE='gs://singlecellldscore/ENSG_gene_annot.txt'
GENE_NAME='ENSGID'
subprocess.call(['/home/sc_enrichement/sc_enrichment-master/main_ldscore.py',
                    '--main-annot-ldscores',INPUT_GENELIST,
		    '--condition-annot-ldscores',COND_GENESET,
		    '--gene-coord-file',COORD_FILE,
                    '--gene-col-name',GENE_NAME,
		    '--summary-stats-files',INPUT_SUMSTAT,
                    '--ldscores-prefix',PREFIX,
		    '--export-ldscore-path',OUT_LDSCORES,
                    '--out',OUT])




