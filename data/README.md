# Public NHANES inputs

Run `Rscript scripts/download_data.R --download` from the repository root.
The downloader obtains the 32 official files listed in `data_manifest.csv` and
checks their SHA-256 hashes. Files stay in this ignored directory by cycle.
No participant data are distributed in the repository.
