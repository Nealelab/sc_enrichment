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

    parser.add_argument('--main-annot-genes',  help = 'Path to file with a list of genes to run as your annotation.This can also have an additional column for continuous annotations.')
    parser.add_argument('--condition-annot-genes', help = 'Path to file with a list of genes to run as a conditional annotation. You can specify multiple comma-separated files. These files have two format: 1) simple gene list 2) genelist + annotation columns. In the latter all the annotation columns will be used for conditional analysis.')
    parser.add_argument('--summary-stats-files', required=True,  help = 'File(s) (already processed with munge_sumstats.py) where to apply partition LDscore, files should end with .sumstats.gz. If multiple files are used, need a comma-separated list.')
    parser.add_argument('--prefix', required=True, help = 'Prefix for main-annot file.')
    parser.add_argument('--out', required=True, help = 'Path to save the results')
    parser.add_argument('--verbose', help="increase output verbosity",action="store_true")
    parser.add_argument('--quantiles', type=int, default=5,required=False, help='If using a continuous annotation,the number of quantiles to split it into for regression.')
    parser.add_argument('--cont-breaks',type=str,required=False,help='Specific boundary points to split your continuous annotation on, comma separated list e.g. 0.1,0.4,0.5,0.6. ATTENTION: if you use negative values add a space in the beginning e.g. <space>-0.1,-0.4,0.5,0.6')

    args = parser.parse_args()
    if not (args.main_annot_genes or args.summary_stats_files or args.prefix or args.out):
        parser.error("You have to specify --main_annot_genes and --summary-stats-files and --prefix and --out")

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
        


def prepare_magma_binary(args):

    """Download and prepare geneset file for MAGMA analysis for binary genelist"""

    with open("/mnt/data/"+ os.path.basename(args.main_annot_genes)) as input:
        content = input.read().splitlines() 
    content.insert(0,args.prefix)
    outlist = (" ".join(map(str, content)))
    with open('/mnt/data/gene_list_for_magma', 'w') as output:
        output.write(outlist)
    logging.info('Wrote geneset for MAGMA: /mnt/data/gene_list_for_magma')

    download_magma()



def prepare_magma_continuous(args):

    """Download and prepare geneset file for MAGMA analysis for continuous genelist"""

    df = pd.read_csv("/mnt/data/"+ os.path.basename(args.main_annot_genes), sep="\t", header=None)

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
        outlist = args.prefix + "_" + labs[ind] + " " + outlist
        with open('/mnt/data/gene_list_for_magma_'+str(ind), 'w') as output:
            output.write(outlist)
        logging.info('Wrote geneset for MAGMA: /mnt/data/gene_list_for_magma_'+str(ind))

    download_magma()


def process_conditional_genesets(cond_file,prefix_cond):

    """ Download, process and save conditional genesets """

    with open("/mnt/data/conditional_genesets/"+ os.path.basename(cond_file)) as input:
        content = input.read().splitlines()
    content.insert(0,prefix_cond)
    outlist = (" ".join(map(str, content)))
    with open('/mnt/data/conditional_gene_list_for_magma_'+prefix_cond, 'w') as output:
        output.write(outlist)


def combine_conditional_genesets():

    """Prepare conditional geneset file for MAGMA analysis"""

    all_gene_lists = glob.glob("/mnt/data/gene_list_for_magma*")
    all_cond_gene_lists = glob.glob("/mnt/data/conditional_gene_list_for_magma_*")
    to_cat = ' '.join(all_cond_gene_lists)
    for file in all_gene_lists:
        os.system('awk 1 ' + file + ' ' + to_cat + ' > ' +  os.path.dirname(file)+'/cond_'+ os.path.basename(file))

    

def run_magma(args,sumstat,phname,prefix_cond_string_dicot,prefix_cond_string_cont,ncol):

    """ Run MAGMA analysis for each sumistat and given the geneset """

    df = pd.read_csv(sumstat, compression='gzip', header=0, sep='\t')
    if 'Z' in list(df.columns.values) and 'P' not in list(df.columns.values):
        df["P"] = 2*st.norm.cdf(-abs(df.Z))
    elif not all(elem in ['SNP', 'P', 'N'] for elem in list(df.columns.values)):
         raise ValueError("Summmary statistics should have column SNP, P, N - or Z if P is not available")
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

        if args.condition_annot_genes and len(prefix_cond_string_dicot)>0 and len(prefix_cond_string_cont)>0:
            subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/cond_gene_list_for_magma',
                                'condition='+ prefix_cond_string_dicot,
                                '--gene-covar',prefix_cond_string_cont,
                                'condition=' + ncol_out,
                                '--out','/mnt/data/magma_results_0_' + phname])

        elif args.condition_annot_genes and len(prefix_cond_string_dicot)>0:
            subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/cond_gene_list_for_magma',
                                'condition='+ prefix_cond_string_dicot,
                                '--out','/mnt/data/magma_results_0_' + phname])

        elif args.condition_annot_genes and len(prefix_cond_string_cont)>0:
            subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/gene_list_for_magma',
                                '--gene-covar',prefix_cond_string_cont,
                                'condition=' + ncol_out,
                                '--out','/mnt/data/magma_results_0_' + phname])
        else:
            subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/gene_list_for_magma',
                                '--out','/mnt/data/magma_results_0_' + phname])

    elif n_magma_genefiles > 1:

        for quantvalue in range(n_magma_genefiles):

            if args.condition_annot_genes and len(prefix_cond_string_dicot)>0 and len(prefix_cond_string_cont)>0:
                subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/cond_gene_list_for_magma_' + str(quantvalue),
                                'condition='+ prefix_cond_string_dicot,
                                '--gene-covar',prefix_cond_string_cont,
                                'condition=' + ncol_out,
                                '--out','/mnt/data/magma_results_' + str(quantvalue) + "_" + phname])

            elif args.condition_annot_genes and len(prefix_cond_string_dicot)>0:
                subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/cond_gene_list_for_magma_' + str(quantvalue),
                                'condition='+ prefix_cond_string_dicot,
                                '--out','/mnt/data/magma_results_' + str(quantvalue) + "_" + phname])

            elif args.condition_annot_genes and len(prefix_cond_string_cont)>0:
                subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/gene_list_for_magma_' + str(quantvalue),
                                '--gene-covar',prefix_cond_string_cont,
                                'condition=' + ncol_out,
                                '--out','/mnt/data/magma_results_' + str(quantvalue) + "_" + phname])
            else:
                subprocess.call(['/home/magma',
                                '--gene-results','/mnt/data/tmp/genes_for_magma_'+ phname + '.genes.raw',
                                '--set-annot','/mnt/data/gene_list_for_magma_' + str(quantvalue),
                                '--out','/mnt/data/magma_results_' + str(quantvalue) + "_" + phname])

    logging.info('MAGMA file generated: '+ '/mnt/data/magma_results_' + phname)



if __name__ == "__main__":

    args = parse_args()
    main_file = args.main_annot_genes

    # Download main annotations
    logging.info('Downloading main annotation file(s):' + main_file)
    subprocess.call(['gsutil','cp',main_file,'/mnt/data/'])

    noun = type_of_file('/mnt/data/' + os.path.basename(main_file))
    logging.info('The type of file that will be used in the analysis: '+noun)


    # Download summary stats
    prefix = args.prefix
    ss_list = args.summary_stats_files.split(',')
    logging.info('The summary statistic(s) to download: ' + ':'.join(ss_list))

    logging.info('Downloading summary statistic(s):' + ':'.join(ss_list))
    subprocess.call(['mkdir','/mnt/data/tmp'])
    subprocess.call(['mkdir','/mnt/data/ss'])
    for ss in ss_list:
        subprocess.call(['gsutil','cp',ss,'/mnt/data/ss/'])

    # Summary statistics
    list_sumstats_file=glob.glob("/mnt/data/ss/*")

    #Prepare genes from main-annot-genes
    if noun=='binary':
        prepare_magma_binary(args)
    elif noun=='continuous':
        prepare_magma_continuous(args)


    # Download and prepare additional geneset for conditioning (if they are specified)
    # And attached them to the output from prepare_magma_*
    prefix_cond_string_dicot=[]
    prefix_cond_string_cont=[]
    ncol_out=None
    if args.condition_annot_genes:
        subprocess.call(['mkdir','/mnt/data/conditional_genesets'])
        cond_files = args.condition_annot_genes.split(',')
        counter = 0
        for k in cond_files:
            # Download
            subprocess.call(['gsutil','cp',k,'/mnt/data/conditional_genesets/'])
            # Get prefix
            prefix_cond = os.path.splitext(os.path.basename(k))[0]
            # Get if file is continuous or not
            noun_cond = type_of_file('/mnt/data/conditional_genesets/' + os.path.basename(k))
            if noun_cond == 'binary':
                process_conditional_genesets(k,prefix_cond)
                prefix_cond_string_dicot.append(prefix_cond)
            if noun_cond == 'continuous':
                counter = counter + 1
                if counter > 1:
                    raise ValueError("No more than 1 continous conditional annotation is allowed")    
                local_file_name='/mnt/data/conditional_genesets/' + os.path.basename(k)
                prefix_cond_string_cont=local_file_name
                ncol=pd.read_csv(local_file_name,delim_whitespace=True,header=None).shape[1]
                ncol_out=','.join([str(x+1) for x in range(ncol-1)])
                logging.info('Continous genesets for adjustment: ' + os.path.basename(prefix_cond_string_cont))

        if prefix_cond_string_dicot:
            prefix_cond_string_dicot = ','.join(prefix_cond_string_dicot)
            logging.info('Binary genesets for adjustment: ' + prefix_cond_string_dicot)
            combine_conditional_genesets()
        
        
 
    # Run MAGMA
    for sumstats in list_sumstats_file:
        phname = os.path.basename(sumstats).replace('.sumstats.gz','')
        run_magma(args,sumstats,phname,prefix_cond_string_dicot,prefix_cond_string_cont,ncol_out)

    # Writing the results
    subprocess.call(['gsutil','-m','cp','/mnt/data/magma_results_*',os.path.join(args.out,"")])

    logging.info('FINITO!')

