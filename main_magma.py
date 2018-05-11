#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import numpy as np
import scipy.stats as st
import argparse
import subprocess
import glob
import sys
import logging
import os
import random
import string
from pybedtools import BedTool
from argparse import Namespace


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--main-annot', required=True, help = 'File(s) containing the gene list to calculate partition h2 or LDscores(s).  If file(s) are detected, LDscores are generated otherwise LDscore(s) are directly used.')
    parser.add_argument('--summary-stats-files', required=True,  help = 'File(s) (already processed with munge_sumstats.py) where to apply partition LDscore, files should end with .sumstats.gz. If multiple files are used, need a comma-separated list.')
    parser.add_argument('--ldscores-prefix', required=True, help = 'Prefix for main-annot file.')
    parser.add_argument('--out', required=True, help = 'Path to save the results')
    parser.add_argument("--verbose", help="increase output verbosity",action="store_true")
    parser.add_argument('--quantiles', type=int, default=5,required=False, help='If using a continuous annotation,the number of quantiles to split it into for regression.')
    parser.add_argument('--cont-breaks',type=str,required=False,help='Specific boundary points to split your continuous annotation on, comma separated list e.g. 0.1,0.4,0.5,0.6. ATTENTION: if you use negative values add a space in the beginning e.g. <space>-0.1,-0.4,0.5,0.6')

    args = parser.parse_args()
    if not (args.main_annot or args.summary_stats_files or args.ldscores_prefix or args.out):
        parser.error("You have to specify --main-annot-file and --summary-stats-files and --ldscores-prefix and --out")

    if (len(args.main_annot.split(',')) != len(args.ldscores_prefix.split(','))):
        parser.error("--main-annot and --ldscores-prefix should be of the same length")

    if (args.cont_breaks):
        args.quantiles = None

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    return args



def type_of_file(file_input):
    '''Want to return a noun that describes file type: rsid/genelist, binary/continuous combination'''
    x = pd.read_csv(file_input,delim_whitespace=True,header=None)
    if x.shape[1] > 1:
       noun = 'continuous'
    else:   
       noun = 'binary'
    if 'rs' in x.loc[0,0]:
        noun = noun + ' rsids'
    else:
        noun = noun + ' genelist'
    return noun


def run_magma(magma_gwas_resuts,out):
    subprocess.call(['/home/magma',
                            '--gene-results',magma_gwas_resuts,
                            '--set-annot','/mnt/data/gene_list_for_magma',
                            '--out',out])

     
def commonprefix(m):

    """Given a list of pathnames, returns the longest common leading component"""

    if not m: return ''
    s1 = min(m)
    s2 = max(m)
    for i, c in enumerate(s1):
        if c != s2[i]:
            return s1[:i]
    return s1


def download_magma():

    """ Download MAGMA files and do initial gene assignment """

    logging.info('Download 1000 genomes reference panel')
    subprocess.call(['gsutil','cp','gs://singlecellldscore/g1000_eur.zip','/mnt/data/'])
    subprocess.call(['unzip','-o','/mnt/data/g1000_eur.zip','-d','/mnt/data/'])
    subprocess.call(['gsutil','cp','gs://singlecellldscore/NCBI37.3.gene.name.loc','/mnt/data/'])
    subprocess.call(['/home/magma',
                                '--annotate',
                                '--snp-loc','/mnt/data/g1000_eur.bim',
                                '--gene-loc','/mnt/data/NCBI37.3.gene.name.loc',
                                '--out','/mnt/data/magma_annotation_1000g_h37'])

def prepare_magma_binary(args,noun):

    """Download and prepare geneset file for MAGMA analysis for binary genelist"""

    with open("/mnt/data/"+ os.path.basename(args.main_annot)) as input:
        content = input.read().splitlines() 
    content.insert(0,args.ldscores_prefix)
    outlist = (" ".join(map(str, content)))
    with open('/mnt/data/gene_list_for_magma', 'w') as output:
        output.write(outlist)
    logging.info('Wrote geneset for MAGMA: /mnt/data/gene_list_for_magma')

    download_magma()


def prepare_magma_continuous(args,noun):

    """Download and prepare geneset file for MAGMA analysis for continuous genelist"""

    df = pd.read_csv("/mnt/data/"+ os.path.basename(args.main_annot), sep="\t", header=None)

    df = pd.read_csv("/mnt/data/allclusters_all_annot_4m.txt", sep="\t", header=None)


    if args.quantiles:
        df["anno_break"] = pd.qcut(df[1], args.quantiles)
        temp_breaks = pd.unique(df["anno_break"])
        n_breaks = len(temp_breaks)
        labs = [str(x).replace(", ","_").replace("(","").replace("]","").replace("[","") for x in temp_breaks]
        logging.info('MAGMA: using the following breaks: '+ "; ".join([str(i) for i in labs]))

    elif args.cont_breaks:
        max_vec = np.max(df[1])
        min_vec = np.min(df[1])

        quantiles_str=args.cont_breaks
        cut_breaks = [float(x) for x in quantiles_str.split(',')]
        name_breaks = list(cut_breaks)

        if np.all(cut_breaks <= max_vec):
            name_breaks.append(max_vec)
            cut_breaks.append(max_vec+1)

        if np.all(cut_breaks >= min_vec):
            name_breaks.append(min_vec)
            cut_breaks.append(min_vec-1)

        name_breaks.sort()
        cut_breaks.sort()
        n_breaks = len(cut_breaks)

        name_breaks[0] = str(min_vec)
        name_breaks[-1] = str(max_vec)
        name_breaks = [str(x) for x in name_breaks]
        labs = [name_breaks[i]+'_'+name_breaks[i+1] for i in xrange(n_breaks-1)]
        labs = labs[::-1]
        logging.info('MAGMA: using the following breaks: '+ "; ".join([str(i) for i in labs]))

        df["anno_break"] = pd.cut(df[1], bins=cut_breaks, labels=labs)


    for ind,anno in enumerate(pd.unique(df["anno_break"])):
        dfout=df.loc[df["anno_break"] == anno]
        outlist = (" ".join(map(str, dfout[0])))
        outlist = args.ldscores_prefix + "_" + labs[ind] + " " + outlist
        with open('/mnt/data/gene_list_for_magma_'+str(ind), 'w') as output:
            output.write(outlist)
        logging.info('Wrote geneset for MAGMA: /mnt/data/gene_list_for_magma_'+str(ind))

    download_magma()

def run_magma(sumstat,phname):

    """ Run MAGMA analysis for each sumistat and given the geneset """

    df = pd.read_csv(sumstat, compression='gzip', header=0, sep='\t')
    df["P"] = 2*st.norm.cdf(-abs(df.Z))
    df = df.dropna(axis=0, how='any')
    df["N"] = df['N'].astype(int)
    dfout = df[['SNP', 'P', 'N']]
    dfout.to_csv('/mnt/data/tmp/extracted_for_magma_'+phname,index=False,sep='\t')
    subprocess.call(['/home/magma',
                            '--bfile','/mnt/data/g1000_eur',
                            '--pval','/mnt/data/tmp/extracted_for_magma_' + phname,
                            'ncol=N',
                            '--gene-annot','/mnt/data/magma_annotation_1000g_h37.genes.annot',
                            '--out','/mnt/data/tmp/genes_for_magma_'+ phname])

    n_magma_genefiles=len(glob.glob('/mnt/data/gene_list_for_magma*'))
    if n_magma_genefiles==1:
        subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/gene_list_for_magma',
                                '--out','/mnt/data/magma_results_0' + phname])
    elif n_magma_genefiles > 1:
        for quantvalue in range(n_magma_genefiles):
            subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/gene_list_for_magma_' + str(quantvalue),
                                '--out','/mnt/data/magma_results_' + str(quantvalue) + "_" + phname])

    logging.info('MAGMA file generated: '+ '/mnt/data/magma_results_' + phname)


if __name__ == "__main__":

    args = parse_args()
    main_file = args.main_annot

    # Download main annotations
    logging.info('Downloading main annotation file(s):' + main_file)
    subprocess.call(['gsutil','cp',main_file,'/mnt/data/'])

    noun = type_of_file('/mnt/data/' + os.path.basename(main_file))
    logging.info('The type of file that will be used in the analysis: '+noun)

    if noun == 'binary genelist' or noun == 'continuous genelist':

        # Download summary stats
        prefix = args.ldscores_prefix
        ss_list = args.summary_stats_files.split(',')
        logging.info('The summary statistic(s) to download: ' + ':'.join(ss_list))

        logging.info('Downloading summary statistic(s):' + ':'.join(ss_list))
        subprocess.call(['mkdir','/mnt/data/tmp'])
        subprocess.call(['mkdir','/mnt/data/ss'])
        for ss in ss_list:
            subprocess.call(['gsutil','cp',ss,'/mnt/data/ss/'])

        # Summary statistics
        list_sumstats_file=glob.glob("/mnt/data/ss/*")

        #If it is just a gene-list run MAGMA
        if noun=='binary genelist':
            prepare_magma_binary(args,noun)
        elif noun=='continuous genelist':
            prepare_magma_continuous(args,noun)
     
        # Run MAGMA
        for sumstats in list_sumstats_file:
            phname = os.path.basename(sumstats).replace('.sumstats.gz','')
            run_magma(sumstats,phname)

        # Writing the results
        logging.info('Results copied to ' + str(args.export_ldscore_path))

        subprocess.call(['gsutil','cp','/mnt/data/magma_results_*',os.path.join(args.out,"")])

        logging.info('FINITO!')
    else:
        logging.info('MAGMA analysis not done because --main-annot was not compatible with it')
