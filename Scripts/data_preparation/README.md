# Data preparation

These scripts document how the four input RData files in `Data/` were built.

They write to `Data/regenerated/` and never overwrite the files in `Data/`.
In this folder they are run by `run_all.R` (step "prepare"), which reads the raw data in `Data/raw_data_notclosed/` (and the non-public files in `Data/raw_data_closed/`, when present) and writes the four files to `<run_dir>/Data/`.
Run them from the repository root.

| Step | Script 
|---|---|---|---|
| 1a | `prep_01a_aging_expression.R` | expression matrices of the aging datasets 
| 1b | `prep_01b_aging_metadata_Alldata_20Sep.R` 
| 2 | `prep_02_Psychedelics_PFC.R` 
| 3 | `prep_03_ExerciseAndAlcohol.R` 
| 4 | `prep_04_Dementia_alldata_meta.R` 

`ensembl_release_map.R` (sourced by steps 1a, 1b, 2 and 3) gives the Ensembl
gene ID -> symbol table of a fixed Ensembl release (see below).
`Scripts/00_LSD_DEGs.R` builds `DE_GSE179379` from the `GSE179379` counts in
`Psychedelics_PFC.RData`.

## Ensembl releases

The original scripts mapped Ensembl gene IDs to gene symbols with biomaRt on
the live Ensembl site (`useMart("ensembl")`), so the symbols depend on the
release that was current at the time. 

| Data | Ensembl release | Original tables checked |
|---|---|---|
| Aging datasets with Ensembl IDs (GSE101521, GSE59288_H, GSE64810, GSE102741, GSE174409) | 105 (Dec 2021) | `anno`, `ids` in `Aging_alldata27Apr.RData` (52,090 and 60,524 pairs) |
| GTEx_PFC (`00_GTEx_filt.R`, Sep 2022) | 107 (Jul 2022) | 2022 `GTEx_RPM_PFC_filt.RData` (rebuilt identically) |
| Psychedelics_PFC (mouse, rat) | 107 | `ids` (mouse, 52,332), `geni_coord` (rat, 33,606) |
| ExerciseAndAlcohol (GSE111273) | 107 | `ids` (43,926) |

`ensembl_release_map(release, species)` returns the same table as
`getBM(attributes = c("ensembl_gene_id", "<db>_symbol"))` on that release. It
reads the BioMart database dump of the release from the Ensembl FTP
(`https://ftp.ensembl.org/pub/release-<N>/mysql/ensembl_mart_<N>/`), which is
permanent, instead of the BioMart archive sites (`dec2021.archive.ensembl.org`,
`jul2022.archive.ensembl.org`), which are often unreachable. The equivalent
biomaRt call is `useEnsembl("ensembl", dataset = ..., version = 105)` (or 107).

`Results/RData/GTEx_RPM_PFC_filt.RData` in this repository was regenerated in
September 2025 with a newer Ensembl release (31,756 genes instead of 30,888);
the 2022 file used for `Alldata_20Sep.RData` is in
`Repository_HT/AgingVsLSD/Results/RData/`. Step 1b therefore rebuilds GTEx_PFC
from the GTEx counts with Ensembl 107 instead of loading that file.

## Verification against the original files (November 2026)

- **Identical** (differences at floating-point level, < 1e-6):
  - all expression matrices of `Alldata_20Sep.RData`, including the Ensembl-based
    datasets and GTEx_PFC;
  - the whole `metadata` table (same rows, order and values);
  - `Dementia_alldata_meta.RData`;
  - the 7 objects of `ExerciseAndAlcohol.RData`;
  - the 14 objects of `Psychedelics_PFC.RData`.
- **Objects saved**: only those read by the analysis scripts.
  - Dropped from `Alldata_20Sep.RData`: 7 expression matrices (their metadata
    is kept).
  - Dropped from `Psychedelics_PFC.RData`: `DE_GSE14720`, and `DE_GSE179379`
    (recomputed by `00_LSD_DEGs.R`).
  - Dropped from `ExerciseAndAlcohol.RData`: GSE105453, which is not in the
    paper.
- **Package versions**:
  - `DE_MDMA` reproduces exactly only with `org.Mm.eg.db` 3.14.0
    (Bioconductor 3.14), which provides the probe symbols of `mouse4302.db`.
  - Recent limma versions need `eBayes(legacy = TRUE)` for GSE60676; this is
    handled in prep_03.
- **GEO files**: microarray gene symbols come from the GPL annotation files
  downloaded in 2022 (in `GEO_DIR` / `RAW_DIR`); newer GEO annotation files
  give different symbols.
- `GSE125681` (column names are the idat file paths) is not used: it has no
  healthy samples.
