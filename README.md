# sc_enrichement
Cloud-based single-cell enrichment analysis

1. Build a docker, this has already been done and available at `gcr.io/ukbb-gay-was/ldscore`
```
docker build -t gcr.io/ukbb-gay-was/ldscore .
gcloud docker -- push gcr.io/ukbb-gay-was/ldscore
```

2. Prepare a tab-separated file containing the inputs for the `dsub` command. See an example in the `/example/` folder.
This fields are mandatories:
```
--input INPUT_SUMSTAT - File containing the gene list to calculate partition h2
--env INPUT_GENELIST - List of comma-separated files (already processed with munge_sumstats.py) where to apply partition LDscore, files should end with .sumstats.gz
--env PREFIX - Prefix for main annotation output
--env OUT - Path to save the results
```

3. Install `dsub` if you have not done yet
```pip install dsub```

4. Run `dsub` command, similar to this:

```
dsub \
	--provider google \
	--project ukbb-gay-was \
	--zones "us-central1-*" \
	--min-ram 4 \
	--logging gs://ldscores/example/log/ \
	--disk-size 200 \
	--image gcr.io/ukbb-gay-was/ldscore \
	--tasks example/submit_list_example.tsv \
	--script run_sc_enrichment.py
```

sometime `dsub` does not recognize the google cloud credential, then you have to `export GOOGLE_APPLICATION_CREDENTIALS="your_google_cloud_service_account_key_file.json"`
