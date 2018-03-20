#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip
import json
import os

def make_annot_files(args):
    print('making gene set bed file')
    GeneSet = pd.read_csv(args.geneset_file, header = None, names = [args.gene_col_name])
    all_genes = pd.read_csv(args.gene_annot, delim_whitespace = True)
    df = pd.merge(GeneSet, all_genes, on = args.gene_col_name, how = 'inner')
    df['START'] = np.maximum(0, df['START'] - args.windowsize)
    df['END'] = df['END'] + args.windowsize
    iter_df = [['chr'+(str(x1).lstrip('chr')), x2, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])]
    genesetbed = BedTool(iter_df).sort().merge()

    print('making annot file')
    df_bim = pd.read_csv(args.bfile_chr + str(args.chrom) + '.bim',
            delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = [['chr'+str(x1), x2, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim)
    annotbed = bimbed.intersect(genesetbed)
    bp = [x.start for x in annotbed]
    df_int = pd.DataFrame({'BP': bp, 'ANNOT':1})
    df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
    df_annot.fillna(0, inplace=True)
    df_annot = df_annot[['ANNOT']].astype(int)
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

