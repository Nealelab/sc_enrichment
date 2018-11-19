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

    parser.add_argument('--main-annot-genes', help = 'Path to file with a list of genes to run as your annotation.This can also have an additional column for continuous annotations.')
    parser.add_argument('--main-annot-rsids',help='Path to file with list of rsids to run as your annotation. This can also have an additional column for a continuous annotation.')
    parser.add_argument('--main-annot-ldscores',help='Path to folder with ldscores for the regression.')
    parser.add_argument('--main-annot-bed',help='Path to file in bed format to run as your annotation. This can also have an additional column for a continuous annotation.')
    
    parser.add_argument('--main-annot-ldcts',help='Path to file that has prefix for what you want your ldcsores to be named "\t" google bucket path to a geneset, one per line to run multiple genesets using --cts flag in ldsc software on one machine.')
    parser.add_argument('--main-annot-ldscores-ldcts',help='Path to file that has prefix of ldscores "\t" gs://path/to/ldscores/prefix.*')

    parser.add_argument('--condition-annot-genes', help = 'Path to file with a list of genes to run as a conditional annotation.This can also have an additional column for continuous annotations.')
    parser.add_argument('--condition-annot-rsids',help='Path to file with list of rsids to run as a conditional annotation. This can also have an additional column for a continuous annotation.')
    parser.add_argument('--condition-annot-ldscores',help='Path to a folder with ldscores to condition on for the regression.')
    parser.add_argument('--condition-annot-bed',help='Path to file in bed format to run as a conditional annotation. This can also have an additional column for a continuous annotation.')
    parser.add_argument('--just-ldscores', action='store_true', default=False, help='Use this flag if you only want to calculate LD-Scores and don\'t need to run a regression. Must be used with --export-ldscore-path')

    parser.add_argument('--summary-stats-files', required=False,  help = 'File(s) (already processed with munge_sumstats.py) where to apply partition LDscore, files should end with .sumstats.gz. If multiple files are used, need a comma-separated list.')
    parser.add_argument('--prefix', required=True, help = 'Prefix that will be used for the ldscore files and the regression output files. If using a --main-annot-*-ldcts flag, this should be a common descriptor of the analyses.')
    parser.add_argument('--out', required=False, help = 'Path to save the regression results.')
    parser.add_argument('--export-ldscore-path', help = 'Path to export the LDscores generated from --main-annot-rsids/genes/bed')
    parser.add_argument('--no-baseline', action='store_true', default=False, help = 'Do not condition on baseline annotations')
    parser.add_argument('--exclude-file', help = 'File in UCSC bed format of regions to exclude in regression')

    parser.add_argument('--windowsize', type=int, default=100000, help = 'size of the window around the gene')
    parser.add_argument('--snp-list-file', default="gs://singlecellldscore/list.txt", help = 'Path of the file containing the list of SNPs to use for the generation of the LD-scores')
    parser.add_argument('--full-report', help = 'Return a full report, including coefficients and enrichment for all annotations.',action="store_true", default=False)
    parser.add_argument('--gene-coord-file', default="gs://singlecellldscore/GENENAME_gene_annot.txt", help = 'Path of the file containing start and end position for each gene, default is ENTREZ')
    parser.add_argument('--gene-col-name', default="GENENAME", help = 'Gene column name in the file specified in --gene-coord-file')

    parser.add_argument('--tkg-weights-folder', default="gs://singlecellldscore/1000G_Phase3_weights_hm3_no_MHC", help = 'Folder containing the chr-specific files with 1000 genomes weights for running LDscore regression')
    parser.add_argument('--tkg-plink-folder', default="gs://singlecellldscore/plink_files", help = 'Folder containing the chr-specific plink files from 1000 genomes to be used to create LDscores')
    parser.add_argument('--tkg-freq-folder', default="gs://singlecellldscore/1000G_Phase3_frq", help = 'Folder containing the chr-specific plink files with 1000 genomes frequencies')
    parser.add_argument('--baseline-ldscores-folder', default="gs://singlecellldscore/baselineLD_v1.1", help = 'Folder containing the baseline chr-specific LDscores to be used for conditioning')
    parser.add_argument("--verbose", help="increase output verbosity",action="store_true")
    
    parser.add_argument('--quantiles', type=int, default=0,required=False, help='If using a continuous annotation,the number of quantiles to split it into for regression. Default is 0. Then the annotation is treated as continuous.')
    parser.add_argument('--cont-breaks',type=str,required=False,help='Specific boundary points to split your continuous annotation on, comma separated list e.g. 0.1,0.4,0.5,0.6. ATTENTION: if you use negative values add a space in the beginning e.g. <space>-0.1,-0.4,0.5,0.6')

    args = parser.parse_args()
    if not args.just_ldscores:
        if not ((args.main_annot_genes or args.main_annot_rsids or args.main_annot_ldscores or args.main_annot_bed or args.main_annot_ldcts or args.main_annot_ldscores_ldcts) or args.summary_stats_files or args.prefix or args.out):
            parser.error("You have to specify --main-annot-* and --summary-stats-files and --prefix and --out")
    else:
        if not ((args.main_annot_genes or args.main_annot_rsids or args.main_annot_ldscores or args.main_annot_bed) or args.prefix or args.export_ldscore_path):
            parser.error("You have to specify --main-annot-* and --prefix and --export-ldscore-path")

    if args.full_report:
        if not (args.main_annot_ldscores_ldcts or args.main_annot_ldcts):
            parser.error("--full-report can only used with --main-annot-ldscores-ldcts or --main-annot-ldcts")

    if (args.cont_breaks):
        args.quantiles = None

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    return args


def random_string(length):
    """ Generate a random string """
    return ''.join(random.choice(string.ascii_letters) for m in range(length))


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


def download_files(args,main_file,ss_list,prefix):

    """Download files for downstream analyses"""

    #Create folders
    logging.info('Creating folders')
    subprocess.call(['mkdir','/mnt/data/ss'])
    subprocess.call(['mkdir','/mnt/data/outld'])
    subprocess.call(['mkdir','/mnt/data/inld'])
    subprocess.call(['mkdir','/mnt/data/tmp'])
    subprocess.call(['mkdir','/mnt/data/genesets/'])

    # Download plink files
    logging.info('Downloading 1000 genomes plink files')
    subprocess.call(['gsutil','-m','cp','-r',args.tkg_plink_folder,'/mnt/data/'])

    # Downlad 1000 genome weights
    logging.info('Downloading 1000 genomes weights for ldscore')
    subprocess.call(['gsutil','-m','cp','-r',args.tkg_weights_folder,"/mnt/data/inld/"])

    # Downlad frequency files
    logging.info('Downloading 1000 genomes frequencies')
    subprocess.call(['gsutil','-m','cp','-r',args.tkg_freq_folder,"/mnt/data/"])


    # Download baseline
    if not args.no_baseline:
        logging.info('Downloading baseline annotation')
        subprocess.call(['gsutil','-m','cp','-r',args.baseline_ldscores_folder,"/mnt/data/inld/"])


    # Download main annotations
    if args.main_annot_ldscores:  
        logging.info('Downloading main annotation LDscores(s):' + main_file)
        subprocess.call(['mkdir','/mnt/data/outld'])
        if '*' in main_file:
            subprocess.call(['gsutil','-m','cp','-r',main_file,'/mnt/data/outld/'])
        else:
            subprocess.call(['gsutil','-m','cp','-r',os.path.join(main_file, "") + '*' ,'/mnt/data/outld/'])
    elif (args.main_annot_genes or args.main_annot_rsids or args.main_annot_bed):
        logging.info('Downloading main annotation file(s):' + main_file)
        subprocess.call(['gsutil','cp',main_file,'/mnt/data/'])
    elif args.main_annot_ldcts:
        logging.info('Downloading main annotation files from list of files provided.')
        subprocess.call(['gsutil','cp',main_file,'/mnt/data/file.ldcts'])
        with open('/mnt/data/file.ldcts','r') as ldcts_file:
            for line in ldcts_file:
                subprocess.call(['gsutil','cp',line,'/mnt/data/genesets/'])
    elif args.main_annot_ldscores_ldcts:
        logging.info('Downloading main annotation files from list of files provided.')
        subprocess.call(['gsutil','cp',main_file,'/mnt/data/file.ldcts'])
        with open('/mnt/data/file.ldcts','r') as ldcts_file:
            for line in ldcts_file:
                path = line.split()[1]
                if '*' in path:
                    subprocess.call(['gsutil','-m','cp','-r',path,'/mnt/data/outld/'])
                else:
                    subprocess.call(['gsutil','-m','cp','-r',os.path.join(path, "") + '*' ,'/mnt/data/outld/'])

    # Download conditional annotations
    if (args.condition_annot_ldscores or args.condition_annot_genes or args.condition_annot_rsids or args.condition_annot_bed):
        if args.condition_annot_bed:
            cond_files = args.condition_annot_bed.split(',')
        if args.condition_annot_genes:
            cond_files = args.condition_annot_genes.split(',')
        if args.condition_annot_rsids:
            cond_files = args.condition_annot_rsids.split(',')
        if args.condition_annot_ldscores:
            cond_files = args.condition_annot_ldscores.split(',') 
        if args.condition_annot_ldscores:
            logging.info('Downloading conditional ldscores annotation(s)')
            subprocess.call(['mkdir','/mnt/data/cond_ldscores'])
            for k in cond_files:
                ts = os.path.join(random_string(7),"")
                subprocess.call(['mkdir','/mnt/data/cond_ldscores/' + ts])
                subprocess.call(['gsutil','-m','cp','-r',os.path.join(k, "") + '*' ,'/mnt/data/cond_ldscores/' + ts])
        else:
            logging.info('Downloading file(s) containing conditional annotations')
            subprocess.call(['mkdir','/mnt/data/outcondld'])
            for k in cond_files:
                subprocess.call(['gsutil','cp',k,"/mnt/data/"])
	    
    # Dowload SNP-list for generating LD-scores
    logging.info('Downloading SNP list for LDscore')
    subprocess.call(['gsutil','cp',args.snp_list_file,'/mnt/data/list.txt'])

    if args.exclude_file:
        logging.info('Downloading file to exclude in regression')
        subprocess.call(['gsutil','cp',args.exclude_file,'/mnt/data/exclude.bed'])
    # Download file mapping SNPs to positions
    logging.info('Downloading file to map genes to positions')
    subprocess.call(['gsutil','cp',args.gene_coord_file,'/mnt/data/GENENAME_gene_annot.txt'])

    # Download summary stats
    if not args.just_ldscores:
        logging.info('Downloading summary statistic(s):' + ':'.join(ss_list))
        for ss in ss_list:
            subprocess.call(['gsutil','cp',ss,'/mnt/data/ss/'])

def prepare_annotations_bed(args,bed_file,outldscore,plink_panel):

    """Prepare LDscores for analysis"""
    logging.info('Creating LDscores')

    for chrom in range(1, 23):
        logging.debug('Running genesets_to_ldscores.py for chr ' + str(chrom) + ' and bed-file ' + str(gene_list))
        subprocess.call(['/home/sc_enrichment/sc_enrichment-master/genesets_to_ldscores.py',
                        '--bed-file',bed_file,
                        '--gene-coord-file',"/mnt/data/GENENAME_gene_annot.txt",
                        '--bfile-chr',plink_panel,
                        '--prefix',outldscore,
                        '--windowsize',str(args.windowsize),
                        '--gene-col-name', str(args.gene_col_name),
                        '--chrom', str(chrom)])

def prepare_annotations_genes(args,gene_list,outldscore,plink_panel):
    """Prepare LDscores for analysis"""
    logging.info('Creating LDscores')

    for chrom in range(1, 23):
        logging.debug('Running genesets_to_ldscores.py for chr ' + str(chrom) + ' and geneset-file ' + str(gene_list))
        subprocess.call(['/home/sc_enrichment/sc_enrichment-master/genesets_to_ldscores.py',
                        '--geneset-file',gene_list,
                        '--chrom',str(chrom),
                        '--gene-coord-file',"/mnt/data/GENENAME_gene_annot.txt",
                        '--bfile-chr',plink_panel,
                        '--prefix',outldscore,
                        '--windowsize',str(args.windowsize),
                        '--gene-col-name', str(args.gene_col_name)])

def prepare_annotations_genes_ldcts(args,gene_list,outldscore,plink_panel,local_prefix):
    """Prepare LDscores for analysis"""
    logging.info('Creating LDscores')

    for chrom in range(1, 23):
        logging.debug('Running genesets_to_ldscores.py for chr ' + str(chrom) + ' and geneset-file ' + str(gene_list))
        subprocess.call(['/home/sc_enrichment/sc_enrichment-master/genesets_to_ldscores.py',
                        '--geneset-file',gene_list,
                        '--chrom',str(chrom),
                        '--gene-coord-file',"/mnt/data/GENENAME_gene_annot.txt",
                        '--bfile-chr',plink_panel,
                        '--prefix',outldscore+local_prefix,
                        '--windowsize',str(args.windowsize),
                        '--gene-col-name', str(args.gene_col_name)])

def prepare_annotations_rsids(args,gene_list,outldscore,plink_panel):
    """Prepare LDscores for analysis"""
    logging.info('Creating LDscores')

    for chrom in range(1, 23):

        logging.debug('Running genesets_to_ldscores.py for chr ' + str(chrom) + ' and rsid-file ' + str(gene_list))
        subprocess.call(['/home/sc_enrichment/sc_enrichment-master/genesets_to_ldscores.py',
                        '--rsid-file',gene_list,
                        '--gene-coord-file',"/mnt/data/GENENAME_gene_annot.txt",
                        '--bfile-chr',plink_panel,
                        '--prefix',outldscore,
                        '--windowsize',str(args.windowsize),
                        '--gene-col-name', str(args.gene_col_name),
                        '--chrom', str(chrom)])      

def calculate_ldscores(args,outldscore,plink_panel,noun):
    for chrom in range(1,23):
        if 'binary' in noun:
            logging.debug('Running ldsc.py for chr ' + str(chrom) )
            subprocess.call(['/home/ldscore/ldsc-kt_exclude_files/ldsc.py',
                            '--l2',
                            '--bfile',plink_panel + str(chrom),
                            '--ld-wind-cm', "1",
                            '--annot',outldscore + '.' + str(chrom) + '.annot.gz',
                            '--thin-annot',
                            '--out', outldscore + "." + str(chrom),
                            '--print-snps',"/mnt/data/list.txt"])
        elif ('continuous' in noun and args.quantiles==0):
            logging.debug('Running ldsc.py for chr ' + str(chrom) )
            subprocess.call(['/home/ldscore/ldsc-kt_exclude_files/ldsc.py',
                            '--l2',
                            '--bfile',plink_panel + str(chrom),
                            '--ld-wind-cm', "1",
                            '--annot',outldscore + '.' + str(chrom) + '.annot.gz',
                            '--thin-annot',
                            '--out', outldscore + "." + str(chrom),
                            '--print-snps',"/mnt/data/list.txt"])
        elif (('continuous' in noun) and args.quantiles):
            try:
                logging.debug('Running ldsc.py for chr ' + str(chrom) )
                subprocess.call(['/home/ldscore/ldsc-kt_exclude_files/ldsc.py',
                                '--l2',
                                '--bfile',plink_panel + str(chrom),
                                '--ld-wind-cm', "1",
                                '--cont-bin',outldscore + '.' + str(chrom) + '.cont_bin.gz',
                                '--cont-quantiles',str(args.quantiles),
                                '--thin-annot',
                                '--out', outldscore + "." + str(chrom)])
            except ValueError:
                sys.exit("The continuous annotation you've entered has non-unique quantile bin edges. Please use --cont-breaks flag instead with user specified bins.")    
        elif (('continuous' in noun) and args.cont_breaks):
            logging.debug('Running ldsc.py for chr ' + str(chrom) )
            subprocess.call(['/home/ldscore/ldsc-kt_exclude_files/ldsc.py',
                            '--l2',
                            '--bfile',plink_panel + str(chrom),
                            '--ld-wind-cm', "1",
                            '--cont-bin',outldscore + '.' + str(chrom) + '.cont_bin.gz',
                            '--cont-breaks',args.cont_breaks,
                            '--thin-annot',
                            '--out', outldscore + "." + str(chrom)])

def calculate_ldscores_ldcts(args,outldscore,plink_panel,local_prefix):
    for chrom in range(1,23):
        logging.debug('Running ldsc.py for chr ' + str(chrom) )
        subprocess.call(['/home/ldscore/ldsc-kt_exclude_files/ldsc.py',
                        '--l2',
                        '--bfile',plink_panel + str(chrom),
                        '--ld-wind-cm', "1",
                        '--annot',outldscore + local_prefix + '.' + str(chrom) + '.annot.gz',
                        '--thin-annot',
                        '--out', outldscore + local_prefix + '.' + str(chrom),
                        '--print-snps',"/mnt/data/list.txt"])    

def commonprefix(m):

    """Given a list of pathnames, returns the longest common leading component"""

    if not m: return ''
    s1 = min(m)
    s2 = max(m)
    for i, c in enumerate(s1):
        if c != s2[i]:
            return s1[:i]
    return s1

def prepare_params_file(args,prefix,name_main_ldscore,params_file='/mnt/data/params.ldcts'):

    """ Save the parameter file containing the name of the ldscores to use for partitioning heritability """
    with open(params_file, 'w') as file:
        logging.debug('Save parameter file with prefix: ' + prefix + ' and ldscore: /mnt/data/outld/' + name_main_ldscore)
        file.write(prefix + "\t" + '/mnt/data/outld/' + name_main_ldscore + '\n')

def prepare_params_file_ldcts(args,main_file,params_file='/mnt/data/params.ldcts'):
    with open(params_file,'w') as file:
        with open('/mnt/data/file.ldcts','r') as ldcts_file:
            for line in ldcts_file:
                local_prefix = line.split()[0]
                file.write(local_prefix + "\t" + '/mnt/data/outld/'+local_prefix+'.'+"\n")


def write_report(report_name,sum_stat,main_panel,cond_panels,outfile):

    """ Write a report about which ldscores panels have been used etc.. """

    with open(report_name, 'a') as file:
        file.write("Summary statistic(s) used: " + sum_stat + '\n')
        file.write("Main panel(s) used: " + main_panel + '\n')
        file.write("Conditional panel(s) used: " + cond_panels + '\n')
        file.write("Main output file(s): " + outfile + '\n')

def ldsc_h2_exclude(infile, params_file, ld_ref_panel, ld_w_panel, tg_f_panel,outfile,exclude_file):

    """Perform partioning hertiability """
    subprocess.call(['/home/ldscore/ldsc-kt_exclude_files/ldsc.py',
                                '--h2-cts',infile,
                                '--ref-ld-chr',ld_ref_panel,
                                '--ref-ld-chr-cts',params_file,
                                '--w-ld-chr',ld_w_panel,
                                '--frqfile-chr',tg_f_panel,
                                '--overlap-annot',
                                '--exclude-file',exclude_file,
                                '--print-all-cts',
                                '--print-coefficients',
                                '--out',outfile])

    logging.info('Running estimate_h2 on: ' + infile)

def ldsc_h2(infile, params_file, ld_ref_panel, ld_w_panel, tg_f_panel,outfile):

    """Perform partioning hertiability """
    subprocess.call(['/home/ldscore/ldsc-kt_exclude_files/ldsc.py',
                                '--h2-cts',infile,
                                '--ref-ld-chr',ld_ref_panel,
                                '--ref-ld-chr-cts',params_file,
                                '--w-ld-chr',ld_w_panel,
                                '--frqfile-chr',tg_f_panel,
                                '--overlap-annot',
                                '--print-all-cts',
                                '--print-coefficients',
                                '--out',outfile])

    logging.info('Running estimate_h2 on: ' + infile)



def ldsc_h2_full(infile, ld_ref_panel, ld_w_panel, tg_f_panel,outfile):

    """Perform partioning hertiability - full report"""
    subprocess.call(['/home/ldscore/ldsc-kt_exclude_files/ldsc.py',
                                '--h2',infile,
                                '--ref-ld-chr',ld_ref_panel,
                                '--w-ld-chr',ld_w_panel,
                                '--frqfile-chr',tg_f_panel,
                                '--overlap-annot',
                                '--print-coefficients',
                                '--out',outfile])

    logging.info('Running estimate_h2 - full report on: ' + infile)





if __name__ == "__main__":

    args = parse_args()

    
    if args.main_annot_bed:
        main_file = args.main_annot_bed
    if args.main_annot_genes:
        main_file = args.main_annot_genes 
    if args.main_annot_rsids:
        main_file = args.main_annot_rsids
    if args.main_annot_ldscores:
        main_file = args.main_annot_ldscores 
    if args.main_annot_ldcts:
        main_file = args.main_annot_ldcts
    if args.main_annot_ldscores_ldcts:
        main_file = args.main_annot_ldscores_ldcts    


    prefix = args.prefix
    if not args.just_ldscores:
        ss_list = args.summary_stats_files.split(',')
    else:
        ss_list=None    

    logging.info('The main annotation file(s) or LDscore(s) to Download: '+ main_file)
    if not args.just_ldscores:
        logging.info('The summary statistic(s) to download: ' + ':'.join(ss_list))

    ld_ref_panel = "No Baseline Panel"
    ld_cond_panel = "No Conditional Panel"

    # Set up the ennviroment
    download_files(args,main_file,ss_list,prefix)
    
    # 1000 genome files
    name_plink = os.path.split(args.tkg_plink_folder)
    name = glob.glob('/mnt/data/' + name_plink[-1] + "/*")
    plink_panel = commonprefix(name)
    logging.debug('plink_panel: ' + plink_panel)

    #Create annotations for main outcome (put each annotation in a different folder)
    #If it is an LDscore put it in a folder and get the name of the LDscore
    if (args.main_annot_rsids or args.main_annot_genes or args.main_annot_bed):
        noun = type_of_file('/mnt/data/' + os.path.basename(main_file))
        logging.info('The type of file that will be used in the analysis: '+noun)
        outldscore='/mnt/data/outld/' + prefix
        if args.main_annot_bed:
            prepare_annotations_bed(args,bed_file='/mnt/data/' + os.path.basename(main_file), plink_panel=plink_panel,outldscore=outldscore)
        elif args.main_annot_genes:
            prepare_annotations_genes(args,gene_list='/mnt/data/' + os.path.basename(main_file), plink_panel=plink_panel,outldscore=outldscore)
        elif args.main_annot_rsids:
            prepare_annotations_rsids(args,gene_list='/mnt/data/' + os.path.basename(main_file), plink_panel=plink_panel,outldscore=outldscore)
        calculate_ldscores(args,outldscore=outldscore,plink_panel=plink_panel,noun=noun)
        name_main_ldscore = prefix + '.'   
    elif (args.main_annot_ldscores):
        temp_name_list =  [os.path.basename(x) for x in glob.glob('/mnt/data/outld/*')]
        name_main_ldscore = commonprefix(temp_name_list)
    elif (args.main_annot_ldcts):
        with open('/mnt/data/file.ldcts','r') as ldcts_file:
            for line,geneset in zip(ldcts_file,os.listdir('/mnt/data/genesets/')):
                local_prefix = line.split()[0]
                prepare_annotations_genes_ldcts(args,gene_list='/mnt/data/genesets/' + geneset,outldscore='/mnt/data/outld/',plink_panel=plink_panel,local_prefix=local_prefix)
                calculate_ldscores_ldcts(args,outldscore='/mnt/data/outld/',plink_panel=plink_panel,local_prefix=local_prefix)

	    
    # If provided, prepare annotation for conditioning gene lists
    if (args.condition_annot_rsids or args.condition_annot_genes or args.condition_annot_bed):
        if args.condition_annot_bed:
            cond_files = args.condition_annot_bed.split(',')
        if args.condition_annot_genes:
            cond_files = args.condition_annot_genes.split(',')
        if args.condition_annot_rsids:
            cond_files = args.condition_annot_rsids.split(',')
        for k in cond_files:
            k_name = os.path.basename(k)
            noun = type_of_file('/mnt/data/' + k_name)
            subprocess.call(['mkdir','/mnt/data/outcondld/' + k_name])
            if args.condition_annot_bed:
                prepare_annotations_bed(args,bed_file='/mnt/data/' + k_name, plink_panel=plink_panel)
            elif args.condition_annot_genes:
                prepare_annotations_genes(args,gene_list='/mnt/data/' + k_name, plink_panel=plink_panel)
            elif args.condition_annot_rsids:
                prepare_annotations_rsids(args,gene_list='/mnt/data/' + k_name, plink_panel=plink_panel)
            calculate_ldscores(args,outldscore='/mnt/data/outcondld/' + k_name + '/' + k_name,plink_panel=plink_panel,noun=noun)   
    
    # Save parameter file
    if not (args.main_annot_ldcts or args.main_annot_ldscores_ldcts):
        prepare_params_file(args,prefix,name_main_ldscore)
    else:
        prepare_params_file_ldcts(args,main_file)

    # Weight panel
    name_w = os.path.split(args.tkg_weights_folder)
    name = glob.glob('/mnt/data/inld/' + name_w[-1] + "/*")
    ld_w_panel = commonprefix(name)
    logging.debug('ld_w_panel: ' + ld_w_panel)


    # Frequency panel
    name_f = os.path.split(args.tkg_freq_folder)
    name = glob.glob('/mnt/data/' + name_f[-1] + "/*")
    tg_f_panel = commonprefix(name)
    logging.debug('tg_f_panel: ' + tg_f_panel)

    # LDscore baseline panel
    if not args.no_baseline:
        name_ldref = os.path.split(args.baseline_ldscores_folder)
        name = glob.glob('/mnt/data/inld/' + name_ldref[-1] + "/*")
        ld_ref_panel = commonprefix(name)
        logging.debug('ld_ref_panel: ' + ld_ref_panel)

    # LDscore conditional panels
    if args.condition_annot_ldscores:
        name_ldcond = glob.glob('/mnt/data/cond_ldscores/*')
        ld_cond_panels_t = []
        for folder in name_ldcond:
            ld_cond_panels_t.append(commonprefix(glob.glob(folder + '/*')))
        logging.debug('ld_cond_panels_t: ' + ':'.join(ld_cond_panels_t))

    # LDscore conditional panels (created from files)
    if (args.condition_annot_rsids or args.condition_annot_genes or args.condition_annot_bed):
        name_ldcond_file = glob.glob('/mnt/data/outcondld/*')
        ld_cond_panels_file_t = []
        for folder in name_ldcond_file:
            ld_cond_panels_file_t.append(commonprefix(glob.glob(folder + '/*')))
        logging.debug('ld_cond_panels_file_t: ' + ':'.join(ld_cond_panels_file_t))

    # Summary statistics
    if not args.just_ldscores:
        list_sumstats_file=glob.glob("/mnt/data/ss/*")

    # Panels for conditioning
    if not args.no_baseline:
         ld_cond_panel = ld_ref_panel
         if args.condition_annot_ldscores:
            ld_cond_panel = ','.join(ld_cond_panels_t + [ld_ref_panel])
         elif (args.condition_annot_rsids or args.condition_annot_genes or args.condition_annot_bed):
            ld_cond_panel = ','.join(ld_cond_panels_file_t + [ld_ref_panel])
    elif (args.no_baseline and (args.condition_annot_ldscores or args.condition_annot_genes or args.condition_annot_rsids or args.condition_annot_bed)):
        if args.condition_annot_ldscores:
            ld_cond_panel = ','.join(ld_cond_panels_t)
        else:
            ld_cond_panel = ','.join(ld_cond_panels_file_t)
    else:
        sys.exit("No baseline panel or conditional panel specified - Interrupting")

    logging.info('The following panel(s) will be used for conditioning: ' + ':'.join([ld_cond_panel]))
    
    if args.just_ldscores:
        logging.info('LDscores copied to ' + str(args.export_ldscore_path))
        subprocess.call(['gsutil','-m','cp','-r','/mnt/data/outld/*',os.path.join(args.export_ldscore_path,"")])
        # Writing report
    

    else:
        # Partitioning heritability
        outfiles_list = []
        for sumstats in list_sumstats_file:
            phname = os.path.basename(sumstats).replace('.sumstats.gz','')
            logging.info('Running partition LDscores for ' + phname)
             # If full report, then run  LDscore for each panel
            if args.full_report:
                with open('/mnt/data/params.ldcts','r') as f:
                    for x in f:
                        x = x.strip().split("\t")
                        ld_cond_panel_full=ld_cond_panel+","+x[1]
                        ld_cond_panel_full=ld_cond_panel_full.replace(" ", "")
                        print('ld_cond_panel_full: ' + ld_cond_panel_full)
                        outfile_full = '/mnt/data/' + phname + '.' + prefix + '.' + x[0] + '.ldsc_full'
                        ldsc_h2_full(infile=sumstats, ld_ref_panel=ld_cond_panel_full, ld_w_panel=ld_w_panel,tg_f_panel=tg_f_panel,outfile=outfile_full)
                        outfiles_list.append('/mnt/data/' + phname + '.' + prefix + '.' + x[0] + '.ldsc_full.cell_type_results.txt')
            else:
                outfiles_list.append('/mnt/data/' + phname + '.' + prefix + '.ldsc.cell_type_results.txt')
                outfile = '/mnt/data/' + phname + '.' + prefix + '.ldsc'
                if not args.exclude_file:
                    ldsc_results = ldsc_h2(infile=sumstats, params_file='/mnt/data/params.ldcts',ld_ref_panel=ld_cond_panel, ld_w_panel=ld_w_panel,tg_f_panel=tg_f_panel,outfile=outfile)
                else:
                    ldsc_results = ldsc_h2_exclude(infile=sumstats, params_file='/mnt/data/params.ldcts',ld_ref_panel=ld_cond_panel, ld_w_panel=ld_w_panel,tg_f_panel=tg_f_panel,outfile=outfile,exclude_file='/mnt/data/exclude.bed')



        # Writing report
        write_report(report_name=prefix + '.report',sum_stat='\t'.join(ss_list),main_panel=main_file, cond_panels=ld_cond_panel, outfile='\t'.join(outfiles_list))

        if args.export_ldscore_path:
            logging.info('LDscores copied to ' + str(args.export_ldscore_path))
            subprocess.call(['gsutil','-m','cp','-r','/mnt/data/outld/*',os.path.join(args.export_ldscore_path,"")])
    
    # Writing the results
        logging.info('Results copied to ' + str(args.export_ldscore_path))
        subprocess.call(['gsutil','cp','/mnt/data/*.ldsc*.cell_type_results.txt',os.path.join(args.out,"")])
        subprocess.call(['gsutil','cp',prefix + '.report',os.path.join(args.out,"")])

    logging.info('FINITO!')
