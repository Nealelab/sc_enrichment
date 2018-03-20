#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
import subprocess
import glob
import sys
import logging
import os
from pybedtools import BedTool
from argparse import Namespace

sys.path.insert(0, '/home/ldscore/ldsc-master/')
sys.path.insert(0, '/home/mtag/mtag-master/')

import ldscore.ldscore as ldsc
import ldscore.sumstats as sumst
from mtag import Logger_to_Logging

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--main-annot-file', required=True, help = 'File containing the gene list to calculate partition h2')
    parser.add_argument('--summary-stats-files', required=True,  help = 'List of comma-separated files (already processed with munge_sumstats.py) where to apply partition LDscore, files should end with .sumstats.gz')
    parser.add_argument('--ldscores-prefix', required=True, help = 'Prefix for main annotation output')
    parser.add_argument('--out', required=True, help = 'Path to save the results')

    parser.add_argument('--no_baseline', action='store_false', default=True, help = 'condition on baseline annotations?')

    parser.add_argument('--condition-annot-file', help = 'file(s) containing the gene list for conditioning')
    parser.add_argument('--condition-annot-ldscores', help = 'comma-separated list of folder locations of ldscores to be used for conditioning')

    parser.add_argument('--export_ldscore_path', help = 'Path where to export the LDscores generated from --main-annot-file')

    parser.add_argument('--gene-col-name', default="GENENAME", help = 'Column name for the files containing a gene list')
    parser.add_argument('--windowsize', type=int, default=100000, help = 'size of the window around the gene')

    parser.add_argument('--snp-list-file', default="gs://ldscores/list.txt", help = 'Location of the file containing the list of SNPs to use for the generation of the LD-scores')
    parser.add_argument('--gene-anno-pos-file', default="gs://ldscores/GENENAME_gene_annot.txt", help = 'Location of the file containing start and end position for each gene')
    parser.add_argument('--tkg-weights-folder', default="gs://constraint_ukbb/1000G_Phase3_weights_hm3_no_MHC", help = 'Folder containing the chr-specific files with 1000 genomes weights for running LDscore regression')
    parser.add_argument('--tkg-plink-folder', default="gs://ldscores/plink_files", help = 'Folder containing the chr-specific plink files from 1000 genomes to be used to create LDscores')
    parser.add_argument('--baseline-ldscores-folder', default="gs://constraint_ukbb/baselineLD_v1.1", help = 'Folder containing the baseline chr-specific LDscores to be used for conditioning')

    parser.add_argument("--verbose", help="increase output verbosity",
                    action="store_true")


    args = parser.parse_args()
    if not (args.main_annot_file or args.summary_stats_files or args.ldscores_prefix or args.out):
        parser.error("You have to specify --main-annot-file and --summary-stats-files and --ldscores-prefix and --out")
    if (args.condition_annot_file and args.condition_annot_ldscores):
        parser.error("You cannto specify both --condition-annot-file and --condition-annot-ldscores, select one")

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    return args


def download_files(args,ss_list):

    """Download files for downstream analyses"""

    # Create folders
    logging.info('Creating folders')
    subprocess.call(['mkdir','/home/ss'])
    subprocess.call(['mkdir','/home/outld'])
    subprocess.call(['mkdir','/home/inld'])
    subprocess.call(['mkdir','/home/tmp'])

    # Download plink files
    logging.info('Downloading 1000 genomes plink files')
    subprocess.call(['gsutil','-m','cp','-r',args.tkg_plink_folder,'/home/'])

    # Downlad 1000 genome weights
    logging.info('Downloading 1000 genomes weights for ldscore')
    subprocess.call(['gsutil','-m','cp','-r',args.tkg_weights_folder,"/home/inld/"])

    # Download baseline
    if args.baseline:
        logging.info('Downloading baseline annotations')
        subprocess.call(['gsutil','-m','cp','-r',args.baseline_ldscores_folder,"/home/inld/"])

    if args.condition_annot_ldscores:
        logging.info('Downloading conditional ldscores annotations')
        subprocess.call(['mkdir','/home/cond_ldscores'])
        cond_ld_list = args.condition_annot_ldscores.split(',')
        for k in cond_ld_list:
            subprocess.call(['gsutil','-m','cp','-r',k,"/home/cond_ldscores/"])

    if args.condition_annot_file:
        subprocess.call(['mkdir','/home/outcondld'])
    
    # Dowload SNP-list for generating LD-scores
    logging.info('Downloading SNP list for LDscore')
    subprocess.call(['gsutil','cp',args.snp_list_file,'/home/list.txt'])

    # Download file mapping SNPs to positions
    logging.info('Downloading file to map genes to positions')
    subprocess.call(['gsutil','cp',args.gene_anno_pos_file,'/home/GENENAME_gene_annot.txt'])

    # Download main annotation file
    logging.info('Downloading main annotation file')
    subprocess.call(['gsutil','cp',args.main_annot_file,'/home/'])

    # Download summary stats
    for ss in ss_list:
        logging.info('Downloading summary statistics')
        subprocess.call(['gsutil','cp',ss,'/home/ss/'])



def prepare_annotations(args,gene_list,outldscore,plink_panel):

    """Prepare LDscores for analysis"""

    logging.info('Creating LDscores')

    for chrom in range(1, 23):

        logging.debug('Running genesets_to_ldscores.py for chr ' + str(chrom) )
        subprocess.call(['/home/sc_enrichement/sc_enrichement-master/genesets_to_ldscores.py',
                        '--geneset-file',gene_list,
                        '--gene-annot',"/home/GENENAME_gene_annot.txt",
                        '--bfile-chr',plink_panel,
                        '--ldscores_prefix','/home/tmp/temp_dscore',
                        '--windowsize',str(args.windowsize),
                        '--gene-col-name', str(args.gene_col_name),
                        '--chrom', str(chrom)])

        logging.debug('Running ldsc.py for chr ' + str(chrom) )
        subprocess.call(['/home/ldscore/ldsc-master/ldsc.py',
                        '--l2',
                        '--bfile',plink_panel + str(chrom),
                        '--ld-wind-cm', "1",
                        '--annot','/home/tmp/temp_dscore.' + str(chrom) + '.annot.gz',
                        '--thin-annot',
                        '--out', outldscore + "." + str(chrom),
                        '--print-snps',"/home/list.txt"])


def commonprefix(m):

    """Given a list of pathnames, returns the longest common leading component"""

    if not m: return ''
    s1 = min(m)
    s2 = max(m)
    for i, c in enumerate(s1):
        if c != s2[i]:
            return s1[:i]
    return s1

def prepare_params_file(args,params_file='/home/params.ldcts'):

    """ Save the parameter file containing the name of the ldscores to use for partitioning heritability """

    with open(params_file, 'a') as file:
        logging.debug('Save parameter file with this content: ' + '/home/outld/' + args.ldscores_prefix + '.')
        file.write(args.ldscores_prefix + "\t" +"/home/outld/" + args.ldscores_prefix + '.\n')



def write_report(report_name,sum_stat,base_panel,cond_panels,outfile):

    """ Write a report about which ldscores panels have been used etc.. """

    with open(report_name, 'a') as file:
        file.write("Summary statistic(s) used: " + sum_stat + '\n')
        file.write("Baseline panel used: " + base_panel + '\n')
        file.write("Conditional panel(s) used: " + cond_panels + '\n')
        file.write("Main output file(s): " + outfile + '\n')


def ldsc_h2(infile, phname, params_file, ld_ref_panel, ld_w_panel,outfile):

    """Perform partioning hertiability using version in MTAG package"""

    args_h2 =  Namespace(out=outfile, 
                         bfile=None,
                         l2=None,
                         extract=None,
                         keep=None,
                         ld_wind_snps=None,
                         ld_wind_kb=None,
                         ld_wind_cm=None,
                         print_snps=None,
                         annot=None,
                         thin_annot=False,
                         cts_bin=None,
                         cts_break=None,
                         cts_names=None,
                         per_allele=False,
                         pq_exp=None,
                         no_print_annot=False,
                         maf=0.05,
                         h2=None,
                         rg=None,
                         ref_ld=None,
                         ref_ld_chr=ld_ref_panel,
                         ref_ld_chr_cts=params_file,
                         w_ld=None,
                         w_ld_chr=ld_w_panel,
                         overlap_annot=False,
                         no_intercept=False,
                         intercept_h2=None,
                         intercept_gencov=None,
                         M=None,
                         two_step=None,
                         chisq_max=None,
                         print_cov=False,
                         print_delete_vals=False,
                         chunk_size=50,
                         pickle=False,
                         invert_anyway=False,
                         yes_really=False,
                         n_blocks=200,
                         not_M_5_50=False,
                         return_silly_things=False,
                         no_check_alleles=False,
                         print_coefficients=True,
                         samp_prev=None,
                         pop_prev=None,
                         frqfile=None,
                         h2_cts=infile,
                         frqfile_chr=None,
                         print_all_cts=False,
                         sumstats_frames=None,
                         rg_mat=False)

    logging.info('Running estimate_h2 on: ' + infile)

    sumst.cell_type_specific(args_h2, Logger_to_Logging())



if __name__ == "__main__":

    args = parse_args()

    ss_list = args.summary_stats_files.split(',')

    logging.info('The list of summary statistics downloaded is: ' + ':'.join(ss_list))

    ld_ref_panel = "No Baseline Panel"
    ld_cond_panel = "No Conditional Panel"

    # Set up the ennviroment
    download_files(args,ss_list)

    # 1000 genome files
    name_plink = os.path.split(args.tkg_plink_folder)
    name = glob.glob('/home/' + name_plink[-1] + "/*")
    plink_panel = commonprefix(name)
    logging.debug('plink_panel: ' + plink_panel)

    # Create annotations for main outcome 
    prepare_annotations(args,gene_list='/home/' + os.path.basename(args.main_annot_file), outldscore='/home/outld/' + args.ldscores_prefix, plink_panel=plink_panel)

    # If provided, prepare annotation for conditioning gene lists
    if args.condition_annot_file:
        cond_list = args.condition_annot_file.split(',')
        for k in cond_list:
            k_name = os.path.basename(k)
            print(k_name)
            subprocess.call(['mkdir','/home/outcondld/' + k_name])
            prepare_annotations(args,gene_list=k,outldscore='/home/outcondld/' + k_name + '/' + k_name, plink_panel=plink_panel)

    # Save parameter file
    prepare_params_file(args)

    # Weight panel
    name_w = os.path.split(args.tkg_weights_folder)
    name = glob.glob('/home/inld/' + name_w[-1] + "/*")
    ld_w_panel = commonprefix(name)
    logging.debug('ld_w_panel: ' + ld_w_panel)

    # LDscore baseline panel
    if args.baseline:
        name_ldref = os.path.split(args.baseline_ldscores_folder)
        name = glob.glob('/home/inld/' + name_ldref[-1] + "/*")
        ld_ref_panel = commonprefix(name)
        logging.debug('ld_ref_panel: ' + ld_ref_panel)

    # LDscore conditional panels
    if args.condition_annot_ldscores:
        name_ldcond = glob.glob('/home/cond_ldscores/*')
        ld_cond_panels_t = []
        for folder in name_ldcond:
            ld_cond_panels_t.append(commonprefix(glob.glob(folder + '/*')))
        logging.debug('ld_cond_panels_t: ' + ':'.join(ld_cond_panels_t))

    # LDscore conditional panels (created from files)
    if args.condition_annot_file:
        name_ldcond_file = glob.glob('/home/outcondld/*')
        ld_cond_panels_file_t = []
        for folder in name_ldcond_file:
            ld_cond_panels_file_t.append(commonprefix(glob.glob(folder + '/*')))
        logging.debug('ld_cond_panels_file_t: ' + ':'.join(ld_cond_panels_file_t))


    # Summary statistics
    list_sumstats_file=glob.glob("/home/ss/*")

    # Panels for conditioning
    if args.baseline:
         ld_cond_panel = ld_ref_panel
         if args.condition_annot_ldscores:
            ld_cond_panel = ','.join(ld_cond_panels_t + [ld_ref_panel])
         if args.condition_annot_file:
            ld_cond_panel = ','.join(ld_cond_panels_file_t + [ld_ref_panel])

    # Partitioning heritability
    outfiles_list = []
    for sumstats in list_sumstats_file:
        phname = os.path.basename(sumstats).replace('.sumstats.gz','')
        outfile = '/home/' + phname + '.' + args.ldscores_prefix + '.ldsc'
        outfiles_list.append('/home/' + phname + '.' + args.ldscores_prefix + '.ldsc.cell_type_results.txt')
        logging.info('Running partition LDscores for ' + phname)
        ldsc_results = ldsc_h2(infile=sumstats, phname=phname, params_file='/home/params.ldcts',ld_ref_panel=ld_cond_panel, ld_w_panel=ld_w_panel, outfile=outfile)

    # Writing report
    write_report(report_name=args.ldscores_prefix + '.report',sum_stat='\t'.join(ss_list),base_panel=ld_ref_panel, cond_panels=ld_cond_panel, outfile='\t'.join(outfiles_list))

    if args.export_ldscore_path:
        logging.info('LDscores copied to ' + str(args.export_ldscore_path))
        subprocess.call(['gsutil','-m','cp','/home/outld/*',args.export_ldscore_path])
    
    # Writing the results
    logging.info('Results copied to ' + str(args.export_ldscore_path))
    subprocess.call(['gsutil','cp','/home/*.ldsc.cell_type_results.txt',args.out])
    subprocess.call(['gsutil','cp',args.ldscores_prefix + '.report',args.out])

    logging.info('FINITO!')
