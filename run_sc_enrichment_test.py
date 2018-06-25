#!/usr/bin/env python

import subprocess

## Inputs ##
INPUT_GENELIST = 'gs://inrich_analyses/inrich_genesets/inrich.00010.GeneSet'
INPUT_SUMSTAT = 'gs://singlecellldscore/PASS/PASS_Schizophrenia.sumstats'
PREFIX = 'kt_test'
OUT = 'gs://inrich_analyses/test/'
OUT_LDSCORES='gs://inrich_analyses/test/'
COND_LDSCORES='gs://singlecellldscorore/entrez_control/'
COORD_FILE='gs://inrich_analyses/ENTREZ_gene_annot.txt'
GENE_NAME='ENTREZ'
subprocess.call(['/home/sc_enrichement/sc_enrichment-master/main_ldscore.py',
                    '--main-annot-genes',INPUT_GENELIST,
		    '--gene-coord-file',COORD_FILE,
                    '--gene-col-name',GENE_NAME,
		    '--summary-stats-files',INPUT_SUMSTAT,
		    '--condition-annot-ldscores',COND_LDSCORES,
                    '--prefix',PREFIX,
		    '--export-ldscore-path',OUT_LDSCORES,
                    '--out',OUT])




