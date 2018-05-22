# sc_enrichement

How to use it

Input Possibilities:


1. Prepare a tab-separated file containing the inputs for the `dsub` command. See an example in `/example/submit_list_example.tsv`.
These fields are mandatories:
```
--env INPUT_GENELIST - Path to either a gene list or list of rsids which can have a second column with a continuous annotation, in which case ldscores will be calculated, OR to a folder with ldscores to use in the regression.
--env INPUT_SUMSTAT - List of comma-separated files (already processed with munge_sumstats.py) where to apply partition LDscore, files should end with .sumstats.gz
--env PREFIX - Prefix for the ldscores files that will be created and the results file from the regression.
--env OUT - Path to save the regression results
```

2. Build a `.py` command to run the analysis. One example is provided in `example/run_sc_enrichment_example.py`

The code should look something like this:

  a) Assign the enviromental variables defined in the file created in step 2.

  ```
  INPUT_GENELIST = os.environ['INPUT_GENELIST']
  INPUT_SUMSTAT = os.environ['INPUT_SUMSTAT']
  PREFIX = os.environ['PREFIX']
  OUT = os.environ['OUT']
  ```
  b) Call the `main_ldscore.py` script, for example:
  ```
  subprocess.call(['/home/sc_enrichement/sc_enrichement-master/main_ldscore.py',
                    '--main-annot',INPUT_GENELIST,
                    '--summary-stats-files',INPUT_SUMSTAT,
                    '--ldscores-prefix',PREFIX,
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
	--logging gs://singlecellldscore/example/log/ \
	--disk-size 200 \
	--image gcr.io/ldscore-data/ldscore \
	--tasks example/submit_list_example.tsv \
	--script example/run_sc_enrichment_example.py
```

sometime `dsub` does not recognize the google cloud credential, then you have to `export GOOGLE_APPLICATION_CREDENTIALS="your_google_cloud_service_account_key_file.json"`

Check the `example/` folder for other examples of submissions programs. 
