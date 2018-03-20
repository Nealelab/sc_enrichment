# sc_enrichement
Cloud-based single-cell enrichment analysis

1. Build a docker, this has already been done and available at `gcr.io/ukbb-gay-was/ldscore`
```
docker build -t gcr.io/ukbb-gay-was/ldscore .
gcloud docker -- push gcr.io/ukbb-gay-was/ldscore
```