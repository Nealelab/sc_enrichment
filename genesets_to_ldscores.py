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

def bed_to_bed(args):
    print('making gene set bed file')
    df = pd.read_csv(args.bed_file,sep='\t',header=None)
    if df.shape[1]==3:
        df.columns=['CHR','START','END']
        binary=True
    elif df.shape[1]>3:
        df.columns=['CHR','START','END','ANNOT']
        binary=False    
    return df, binary

def rsids_to_bed(args):
    print('making rsid list into bed file')
    df_bim = pd.read_csv(args.bfile_chr + str(args.chrom) + '.bim',
    delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    GeneSet = pd.read_csv(args.rsid_file, header = None,sep='\t')
    if GeneSet.shape[1] == 1:
        GeneSet.columns = [args.gene_col_name]
        binary=True
    else: 
        GeneSet.columns = [args.gene_col_name,'ANNOT']
        binary=False
    df_rs = pd.merge(GeneSet,df_bim,how='inner',left_on=args.gene_col_name,right_on='SNP')
    df = df_rs[['CHR','BP']]
    df = df.rename(columns={'BP':'START'})
    df['END'] = df['START']
    if binary == False:
        df['ANNOT'] = df_rs['ANNOT']
    return df, binary

def genes_to_bed(args):
    print('making gene set bed file')
    GeneSet = pd.read_csv(args.geneset_file, header = None,sep='\t')
    if GeneSet.shape[1] == 1:
        GeneSet.columns = [args.gene_col_name]
        binary=True
    else: 
        GeneSet.columns = [args.gene_col_name,'ANNOT']
        binary=False
    all_genes = pd.read_csv(args.gene_coord_file, delim_whitespace = True)
    df = pd.merge(GeneSet, all_genes, on = args.gene_col_name, how = 'inner')
    df['START'] = np.maximum(0, df['START'] - args.windowsize)
    df['END'] = df['END'] + args.windowsize
    
    return df, binary

def make_annot_files(args,df,binary):
    df = df.sort_values(by=['CHR','START'])
    if binary==True:
        iter_df = [['chr'+(str(x1)).lstrip('chr'), x2, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])]
        genesetbed = BedTool(iter_df).sort().merge()
    elif binary==False:
        iter_df = [['chr'+(str(x1).lstrip('chr')), int(x2), int(x3),'annot',str(x4)] for (x1,x2,x3,x4) in np.array(df[['CHR', 'START', 'END','ANNOT']])]
        genesetbed = BedTool(iter_df).sort()
    
    print('making annot file')
    df_bim = pd.read_csv(args.bfile_chr + str(args.chrom) + '.bim',
        delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = [['chr'+str(x1), x2, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim).sort()
    if binary == True:
        annotbed = bimbed.intersect(genesetbed)
        bp = [x.start for x in annotbed]
        df_int = pd.DataFrame({'BP': bp, 'ANNOT':1})
        df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
        df_annot.fillna(0, inplace=True)
        df_annot = df_annot[['ANNOT']].astype(int)
    else:
        annotbed = bimbed.map(genesetbed,c=5,o='mean',null=0).to_dataframe()
        bp = annotbed.start
        annot = annotbed.name
        df_int = pd.DataFrame({'BP': bp, 'ANNOT':annot})
        df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
        num_snps_final = df_annot.ANNOT.count()
        df_annot.fillna(0, inplace=True)
        df_annot = df_annot[['ANNOT']].astype(float)
        df_bim['ANNOT'] = df_annot[['ANNOT']]
        cont_annot = df_bim[['SNP','ANNOT']]
        cont_annot_file = args.prefix+'.'+str(args.chrom)+'.cont_bin.gz'
        with gzip.open(cont_annot_file,'wb') as f:
            cont_annot.to_csv(f,sep="\t",index=False,header=None)
    
    annot_file = args.prefix+'.'+str(args.chrom)+'.annot.gz' 
    with gzip.open(annot_file, 'wb') as f:
        df_annot.to_csv(f, sep = "\t", index = False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--geneset-file', help = 'location of the Genset file')
    parser.add_argument('--rsid-file',help= 'location of the rsids file')
    parser.add_argument('--bed-file',type=str, help='the UCSC bed file with the regions that make up your annotation')
    parser.add_argument('--gene-coord-file', help = 'location of the file mapping genes to positions')
    parser.add_argument('--bfile-chr', help = 'plink file for creating annot')
    parser.add_argument('--prefix', help = 'path and prefix of the ldscore')
    parser.add_argument('--chrom',type=int,help='chromosome')
    parser.add_argument('--windowsize', type=int, default=100000, help = 'size of the window around the gene')
    parser.add_argument('--dont-make-ldscores', action='store_true', default=False)
    parser.add_argument('--gene-col-name', default = 'GENENAME', help = 'which column to use as Gene Name')

    args = parser.parse_args()
    if args.geneset_file or args.rsid_file or args.bed_file is not None:
        if args.geneset_file:
            df, binary = genes_to_bed(args)
        if args.rsid_file:
            df, binary = rsids_to_bed(args) 
        if args.bed_file:
            df, binary = bed_to_bed(args)

        make_annot_files(args,df,binary)

