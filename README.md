# sc_enrichement

This is a cloud-based pipeline that uses the job submission tool [dsub](https://github.com/DataBiosphere/dsub) to run [stratified LD-Score Regression](https://www.biorxiv.org/content/early/2015/01/23/014241) and [MAGMA](http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004219) in a parallelized way. The pipeline is split into two scripts, main_ldscore.py and main_magma.py. 

Flags for ```main_ldscore.py```:

```
--main-annot-genes 
This flag accepts as input a list of genes over which you want to partition heritability.
This gene list will be converted into a per snp annotation. The file can have a single 
column indicating it will be a binary annotation or contain a secondary column that 
is a quantitative annotation for each gene. 
```
```
--main-annot-rsids
This flag accepts as input a list of rsids over which you want to partition heritability.
The file can have a single column indicating it will be a binary annotation or contain a
secondary column that is a quantitative annotation for each SNP.
```
```
--main-annot-ldscores
This flag accepts a path to a folder that contains pre-calculated ldscores for an 
annotation. This ldscores will be run directly into the regression portion of the
pipeline to produce partitioned heritability results.
```
```
--main-annot-bed
This flag accepts a file that is in UCSC bed file format for regions over which 
you want to partition heritability. There can be a 4th column that is a continuous
annotation.
```
```
--main-annot-ldcts
This flag accepts a file that has two columns, the first is the prefix for your ldscores for a geneset,
and the second is the google bucket path to the corresponding geneset. One geneset per line. This allows
the user to take advantage of the --cts flags within LDSC software to run many genesets on one VM.
```
```
--main-annot-ldscores-ldcts
This flag accepts a file that has two columns, the first is the prefix of the ldscores for a geneset,
and the second the google bucket path to the corresponding ldscores. On set of ldscores per line.
This allows the user to take advantage of the --cts flags within LDSC when you already have
ldscores calculated. E.g. test_analyses.GeneSet1 gs://test_analyses/ldscores/test_analyses.GeneSet1.*
```
```
--condition-annot-genes/--condition-annot-rsids/--condition-annot-ldscores/--condition-annot-bed
These flags work the same as the --main-annot-* flags but are used when you want 
to condition the regression on another annotation.
```
```
--just-ldscores
This flag allows you to just calculate ldscores for a particular annotation. 
If given, this flag will prevent any regression from being run.
```
```
--sumary-stats-files
A comma separated list of summary statistics files ending in .sumstats.gz that 
have already been processed using munge_sumstats.py.
```
```
--prefix
Prefix that will be used when naming ldscore files and regression output files.
```
```
--out
Path to folder to save regression results to.
```
```
--export-ldscore-path
Path to folder to save ldscores to. If given this flag will copy the ldscores to 
the path, if not ldscore files will not be written out.
```
```
--gene-coord-file
Path to file that has gene coordinates. Format is GENE CHR START END including the header.
If not using the default (ENSGID based) file, you need to include --gene-col-name flag 
to indicate what the first column of your --gene-coord-file is called. 
e.g if your --gene-coord-file is headed as such: ENTREZ CHR START END you would indicate
 --gene-col-name ENTREZ
```

Steps to run the pipeline:

1. Prepare a tab-separated file containing the inputs for the `dsub` command. See an example in `/example/submit_list_example.tsv`. These environmental variables are then read in by the script called by `dsub` as explained below.
Depending on your analysis, these fields can change but below is an example:
```
--env INPUT_MAIN - This will provide your main annotation, can be path to gene list, rsids,
                   ldscore folder, bed file or ldcts file with genesets or ldscores.
		   In this example it is a gene list.
--env INPUT_SUMSTAT - List of comma-separated files (already processed with munge_sumstats.py) 
                      where to apply partition LDscore.
--env PREFIX - Prefix for the ldscores files that will be created and the results file 
               from the regression.
--env OUT - Path to save the regression results
```

Example:
```
--env INPUT_MAIN    --env INPUT_SUMSTAT    --env PREFIX    -env OUT
gs://singlecellldscore/example/example.geneset    gs://singlecellldscore/example/asd_summary_stats.sumstats.gz,gs://singlecellldscore/example/scz_summary_stats.sumstats.gz    example    gs://singlecellldscore/example/    
```
2. Build a `.py` command to run the analysis. One example is provided in `example/run_sc_enrichment_example.py`

The code should look something like this:

  a) Assign the enviromental variables defined in the file created in step 2.

  ```
  INPUT_MAIN = os.environ['INPUT_MAIN']
  INPUT_SUMSTAT = os.environ['INPUT_SUMSTAT']
  PREFIX = os.environ['PREFIX']
  OUT = os.environ['OUT']
  ```
  b) Call the `main_ldscore.py` script, for example:
  ```
  subprocess.call(['/home/sc_enrichement/sc_enrichement-master/main_ldscore.py',
                    '--main-annot-genes',INPUT_MAIN,
                    '--summary-stats-files',INPUT_SUMSTAT,
                    '--prefix',PREFIX,
                    '--out',OUT,
                    '--verbose'])
  ```

There are many other options that can be used. For another example check `example/run_sc_enrichment_example_advanced.py`.


3. Install `dsub` if you have not done yet
```pip install dsub```

4. Run `dsub` command, similar to this:

```
dsub \
	--provider google \
	--project ldscore-data \
	--zones "us-central1-*" \
	--min-ram 4 \
	--min-cores 4 \
	--logging gs://singlecellldscore/example/log/ \
	--disk-size 100 \
	--image gcr.io/ldscore-data/ldscore \
	--tasks example/submit_list_example.tsv \
	--script example/run_sc_enrichment_example.py
```

sometime `dsub` does not recognize the google cloud credential, then you have to `export GOOGLE_APPLICATION_CREDENTIALS="your_google_cloud_service_account_key_file.json"`

Check the `example/` folder for other examples of submissions programs. 
