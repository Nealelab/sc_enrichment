#!/usr/bin/env python

from __future__ import print_function,division
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip
import json
import os
import logging


def make_annot_files(args):
    print('making gene set bed file')
    GeneSet = pd.read_csv(args.geneset_file, header = None,sep='\t')
    if GeneSet.shape[1] == 1:
	GeneSet.columns = [args.gene_col_name]
	binary=True
    else: 
        GeneSet.columns = [args.gene_col_name,'ANNOT']
        binary=False
    
    df_bim = pd.read_csv(args.bfile_chr + str(args.chrom) + '.bim',
	        delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    if 'rs' in GeneSet[args.gene_col_name].all():
        df_rs = pd.merge(GeneSet,df_bim,how='inner',left_on=args.gene_col_name,right_on='SNP')
        df = df_rs[['CHR','BP']]
	df = df.rename(columns={'BP':'START'})
	df['END'] = df['START']
        num_snps = len(df)
	if binary == False:
	    df['ANNOT'] = df_rs['ANNOT']
    else:
        num_snps=None
        all_genes = pd.read_csv(args.gene_annot, delim_whitespace = True)
        num_genes = len(all_genes)
        df = pd.merge(GeneSet, all_genes, on = args.gene_col_name, how = 'inner')
        num_final_genes = len(df)
        if (num_genes/num_final_genes) < 0.1:
            logging.warning("The number of genes in your annotation after merging is less than 10% of the original list.")
        df['START'] = np.maximum(0, df['START'] - args.windowsize)
        df['END'] = df['END'] + args.windowsize
    if binary == True:
        df = df.sort_values(by=['CHR','START'])
        iter_df = [['chr'+(str(x1)).lstrip('chr'), x2, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])]
        genesetbed = BedTool(iter_df).sort().merge()
    else:
        df = df.sort_values(by=['CHR','START'])
	iter_df = [['chr'+(str(int(x1)).lstrip('chr')), int(x2), int(x3),'annot',float(x4)] for (x1,x2,x3,x4) in np.array(df[['CHR', 'START', 'END','ANNOT']])]
        genesetbed = BedTool(iter_df).sort()
    print('making annot file')
    iter_bim = [['chr'+str(x1), x2, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    if binary == True:
        bimbed = BedTool(iter_bim)
        annotbed = bimbed.intersect(genesetbed)
        bp = [x.start for x in annotbed]
        df_int = pd.DataFrame({'BP': bp, 'ANNOT':1})
        df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
        df_annot.fillna(0, inplace=True)
        df_annot = df_annot[['ANNOT']].astype(int)
        num_snps_final = sum(df_annot.values)
        if num_snps is not None:
            if (num_snps/num_snps_final) < 0.1:
                logging.warning("The number of SNPs in your annotation after merging is less than 10% of the original list.")
    else:
        bimbed = BedTool(iter_bim).sort()
        annotbed = bimbed.map(genesetbed,c=5,o='mean',null=0).to_dataframe()
        bp = annotbed.start
        annot = annotbed.score
        df_int = pd.DataFrame({'BP': bp, 'ANNOT':annot})
        df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
        num_snps_final = df_annot.ANNOT.count()
        df_annot.fillna(0, inplace=True)
        df_annot = df_annot[['ANNOT']].astype(float)
        df_bim['ANNOT'] = df_annot[['ANNOT']]
        cont_annot = df_bim[['SNP','ANNOT']]
        cont_annot_file = args.ldscores_prefix+'.'+str(args.chrom)+'.cont_bin.gz'
        with gzip.open(cont_annot_file,'wb') as f:
            cont_annot.to_csv(f,sep="\t",index=False,header=None)
    
    annot_file = args.ldscores_prefix+'.'+str(args.chrom)+'.annot.gz' 
    with gzip.open(annot_file, 'wb') as f:
        df_annot.to_csv(f, sep = "\t", index = False)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--geneset-file', help = 'location of the Genset file')
    parser.add_argument('--gene-annot', help = 'location of the file mapping genes to positions')
    parser.add_argument('--bfile-chr', help = 'plink file for creating annot')
    parser.add_argument('--ldscores_prefix', help = 'path and prefix of the ldscore')
    parser.add_argument('--chrom', type=int)
    parser.add_argument('--windowsize', type=int, default=100000, help = 'size of the window around the gene')
    parser.add_argument('--dont-make-ldscores', action='store_true', default=False)
    parser.add_argument('--gene-col-name', default = 'GENENAME', help = 'which column to use as Gene Name')

    args = parser.parse_args()

    make_annot_files(args)

