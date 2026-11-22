# AgingVsLSD: analysis pipeline

This folder contains everything needed to rerun the analyses of the paper, from the
raw data to the supplementary tables and figures, without an internet connection:
- the input data (`Data/`), except the non-public raw data (section 4);
- the four input RData files used for the paper (`Data/original_RData_notclosed/`), without
  the raw tables of the non-public datasets but with their DE tables (section 4);
- the code (`Scripts/`);
- the R packages (`Data/packages/`) and the scripts that install and run them.

The analysis always runs on the deposited input RData of `Data/original_RData_notclosed/`,
which also contain the DE tables of the non-public datasets, so every analysis script is
run. The MSBB gene correlations are read from `Data/closed_data_intermediate/` (section 4).
The RData rebuilt from the raw data are only compared with the deposited ones.

Contents:
1. How to run
2. Folder layout
3. Input data
4. Non-public data
5. Data preparation: the four input RData files
6. Analysis scripts
7. Intermediate files
8. R environments
9. Random numbers
10. Numerical precision

## 1. How to run

Requirements: Windows, **R 4.5**, about 10 GB of free disk space.

Open a terminal in this folder and run:

```
Rscript Scripts/Reproducibility/install_environment.R
Rscript run_all.R C:/AgingVsLSD_run
```

1. `Scripts/Reproducibility/install_environment.R` installs every R package from the files in `Data/packages/`
   into `C:/Users/<user>/AgingVsLSD_Rlib` (or into the folder given as argument).
   - It takes about 10 minutes and is needed only once.
   - The user's own R library is not changed.
2. `run_all.R` does everything else, in about 8 hours, and writes only into the
   folder given as argument (here `C:/AgingVsLSD_run`; use a folder outside Dropbox,
   about 8 GB). It runs three steps:

   | Step | What it does | Output in `C:/AgingVsLSD_run` |
   |---|---|---|
   | `prepare` | rebuilds the four input RData files from the raw data in `Data/raw_data_notclosed/` (`Scripts/data_preparation`, section 5); `Psychedelics_PFC.RData` without the DE tables of the non-public datasets (section 4) | `rebuilt_RData/*.RData` |
   | `verify` | compares the rebuilt input RData with the deposited ones (`Data/original_RData_notclosed/`), object by object | `comparison/data_preparation_verification.csv` |
   | `analysis` | runs the analysis scripts (section 6), each in a fresh R session with its own R library (section 8), on the deposited input RData (`Data/original_RData_notclosed/`), never on the rebuilt ones | `Results/`, `figures_NES/`, `figures_dementia_cohen/` |

   The steps can also be run one at a time, for example
   `Rscript run_all.R C:/AgingVsLSD_run analysis` (about 5 h).
   The option `--rdata=notclosed` (the default here) chooses the files `verify` compares
   with; `--rdata=full` selects `Data/original_RData/`, the complete original files, which
   are not distributed (section 4).

All R processes get a non-working proxy, so nothing can be downloaded: every input
is in this folder. Set `ALLOW_NETWORK=1` to allow connections.

Logs, `sessionInfo()` and the objects left in the workspace of each script are in
`C:/AgingVsLSD_run/logs/`.

**Verification of the input RData** (`C:/AgingVsLSD_run/comparison/data_preparation_verification.csv`):
one row per object of the rebuilt input RData files, with the dimensions, whether row and
column names are the same, and the largest absolute difference between the rebuilt and
the deposited object (section 10).

## 2. Folder layout

| Folder or file | Content |
|---|---|
| `run_all.R`, `Scripts/Reproducibility/install_environment.R` | the two commands of section 1 |
| `Data/original_RData_notclosed/` | the four input RData files used for the paper, deposited as intermediate files: `Alldata_20Sep.RData`, `ExerciseAndAlcohol.RData`, `Dementia_alldata_meta.RData` unchanged, `Psychedelics_PFC.RData` without the raw tables of the non-public datasets (section 4); the `verify` step compares the rebuilt files with them. Also the other distributed intermediate files (section 7) |
| `Data/raw_data_notclosed/` | public data: raw data of the data preparation (four subfolders) and the other inputs of the analysis scripts (section 3) |
| `Data/intermediate_RData/` | outputs of the main pipeline read by the threshold-sensitivity analyses (section 7) |
| `Data/homologs/` | human-rat and human-mouse homolog tables |
| `Data/ensembl/` | Ensembl ID-symbol tables of releases 105, 107 and 115 |
| `Data/closed_data_intermediate/` | `MSBB_gene_correlations.RData`, read by `13_AD_validation.R` without the MSBB data (section 4) |
| `Scripts/` | the analysis scripts and `utilities.R` |
| `Scripts/data_preparation/` | the scripts that build the four input RData files |
| `Scripts/Reproducibility/` | the pipeline runner (`run_pipeline.R`, `run_script.R`), the installation of the R libraries (`install_environment.R`), the verification of the input RData (`verify_data_preparation.R`), the script that made `Data/original_RData_notclosed/` (`make_notclosed_RData.R`, not part of the pipeline), the public release (`make_public_release.R`), the Ensembl release check (`check_ensembl_releases.R`) |
| `Data/envs/` | the R library of each script (`environments.csv`) and the package versions of each library (`*_versions.csv`) |
| `Data/packages/` | all R packages, for the offline installation |

The scripts read these paths directly (for example
`Data/original_RData_notclosed/Alldata_20Sep.RData`, `Data/homologs/Human rat homologs.txt`,
`Data/raw_data_notclosed/Invivo RNAseq/`). `run_all.R` links the folders of `Data/` (except
`original_RData/` and `packages/`) into `<run_dir>/Data/`.

## 3. Input data

**Raw data of the data preparation** (subfolders of `Data/raw_data_notclosed/`):

| Folder | Content |
|---|---|
| `aging_GEO/` | human aging datasets: GEO series matrices, CEL files, supplementary files and GPL annotation files |
| `psychedelics_GEO/` | psychedelic and antidepressant datasets (GEO) |
| `exercise_alcohol_GEO/` | exercise, enriched environment and alcohol datasets (GEO) |
| `dementia_GEO/` | human dementia datasets (GEO) |

The GEO files were downloaded in 2022-2024. GEO updates the platform annotation
files: newer downloads give slightly different gene symbols.

**Other inputs of the analysis scripts** (`Data/raw_data_notclosed/`, unless stated):
- homolog tables (`Data/homologs/`: `Human rat homologs.txt`, `Human mouse homologs.txt`);
- GTEx v8 (dbGaP phs000424.v8, release 2017-06-05): open-access files of the GTEx
  portal, `gene_reads_2017-06-05_v8_brain_frontal_cortex_ba9.gct.gz` and the two
  `GTEx_Analysis_v8_Annotations_*` files;
- `Invivo RNAseq/`: read counts of the LSD validation experiment (in vivo, mouse PFC), as
  deposited in GEO (GSE263662): 15 samples, two sequencing lanes each (30 files,
  `ASA_<n>_S<sample>_L00<lane>.txt`), the samples used in `05_DEvsAging_validation.R`, listed in `Samples_LSD.xlsx`;
- `In vitro RNA-seq_all/`: in vitro RNA-seq;
- GSE75772 (rat aging), GSE25219 and GSE36192 (validation datasets): series matrices
  and platform tables (`GPL5175.soft.gz`, `GPL6947.annot.gz`);
- `ACEL-17-e12819-s002.csv`;
- `cmap/` (in `Data/original_RData_notclosed/`): CMap (LINCS) input of
  `08_2_random_drugs.R`, prepared from the LINCS data with the clue.io API
  (`CMAP_NPC_mean.RData`, `compounds_NPC.RData`).

**Non-public data**: MDMA, DMT, pharmahuasca, harmaline and MSBB (section 4).

**Ensembl releases.** Ensembl gene IDs are mapped to gene symbols with the tables in
`Data/ensembl/`, taken from the permanent BioMart dumps on the Ensembl FTP site
(`Scripts/data_preparation/ensembl_release_map.R`). The original analyses used the
live Ensembl BioMart; the releases below were identified by comparing the ID-symbol
tables of the original analyses with each release (identical pair by pair;
`Scripts/Reproducibility/check_ensembl_releases.R`, tables in `Data/ensembl/original_workspace_tables/`).

| Data | Ensembl release |
|---|---|
| aging datasets with Ensembl IDs (GSE101521, GSE59288_H, GSE64810, GSE102741, GSE174409) | 105 (Dec 2021) |
| GTEx_PFC, psychedelics (mouse, rat), exercise (GSE111273) | 107 (Jul 2022) |
| analysis scripts `00_GTEx_filt` (Fig S16), `05` (mouse), `09` (rat GSE75772), `13` (MSBB) | 115 (Sep 2025) |

**Probe annotation.** Microarray gene symbols come from the GPL annotation files in
`Data/raw_data_notclosed/`. The MDMA arrays use `org.Mm.eg.db` 3.14.0 with `mouse4302.db` 3.13.0
(later versions rename a few hundred genes).

## 4. Non-public data

Three datasets used in the paper are not distributed with this folder: they do not
belong to the authors of this work, or are controlled-access data. The pipeline does
not need them: it uses their DE tables (`Data/original_RData_notclosed/Psychedelics_PFC.RData`)
and the MSBB gene correlations (`Data/closed_data_intermediate/`), computed from them
for the paper.

| Data | Study | How to obtain it | Where to put it |
|---|---|---|---|
| MDMA, mouse PFC, Affymetrix Mouse 430 2.0 (`PF1.CEL`-`PF9.CEL`) | Fernàndez-Castillo et al. 2012, *Genes Brain Behav* 11, 38-51 (Table S5) | not deposited in public repositories: request it from the authors of the study | `Data/raw_data_closed/psychedelics/MDMA data/PREFRONTAL CORTEX/` |
| DMT, pharmahuasca, harmaline, rat PFC, RNA-seq read counts (one table, 16 samples) | Kelley et al. 2022, *ACS Chem. Neurosci.* 13, 257-274 (Table S5) | not deposited in public repositories: request it from the authors of the study | `Data/raw_data_closed/psychedelics/readcount copy 2.xlsx` |
| MSBB (Mount Sinai Brain Bank) RNA-seq: `MSBB_Normalized_counts_(CQN).tsv`, `RNAseq_Harmonization_MSBB_combined_metadata.csv` | AD Knowledge Portal, RNAseq Harmonization study (Synapse `syn21241740`): `MSBB_Normalized_counts_(CQN).tsv` = `syn27068756` version 1 (MD5 `4665daf23145048d0daecf113d6f9fdb`), `RNAseq_Harmonization_MSBB_combined_metadata.csv` = `syn27000243` version 2 (MD5 `e225b3396859528f1a9a6eea2927b8a7`, same content as version 1); downloaded on 9 January 2023, MD5 checked against Synapse | controlled access: request access through the AD Knowledge Portal (https://adknowledgeportal.synapse.org) | `Data/AD portal/` |

The processing of these data is described in the Methods of the paper (MDMA, DMT,
pharmahuasca, harmaline: "Comparing the aging-signature-opposition/-concordance of
psychoplastogen with that of aging-concordant and aging-reverting comparators" and
Table S5; MSBB: "External validations") and in the scripts
`Scripts/data_preparation/prep_02_Psychedelics_PFC.R` and `Scripts/13_AD_validation.R`.

**What changes without them.** The pipeline checks which data are present and runs
accordingly:

| Missing data | What is used instead | Not recomputed from the raw data |
|---|---|---|
| MDMA or DMT/pharmahuasca/harmaline | nothing changes for the analysis, which reads `DE_MDMA`, `DE_DMT`, `DE_pharm`, `DE_harm` from the deposited `Data/original_RData_notclosed/Psychedelics_PFC.RData`. `prep_02` builds `Psychedelics_PFC.RData` without them | the four DE tables (the `verify` step lists them as not rebuilt) |
| MSBB | `Data/closed_data_intermediate/MSBB_gene_correlations.RData`: gene-wise correlations of expression with age and with the Braak score (`cor_age`, `cor_braak`; no sample data), read by `13_AD_validation.R`. Fig S11B is computed from them (same GSEA as with the MSBB files) | the two correlation vectors |

**Input RData.** `Data/original_RData_notclosed/` holds the four input RData files used for
the paper. `Alldata_20Sep.RData`, `ExerciseAndAlcohol.RData` and
`Dementia_alldata_meta.RData` are the original files. `Psychedelics_PFC.RData` is the
original file without the two tables of raw non-public data (`MDMA`, the RMA expression
matrix of the MDMA arrays; `DMT`, the read counts of DMT, pharmahuasca and harmaline);
every other object, the DE tables included, is unchanged
(`Scripts/Reproducibility/make_notclosed_RData.R`). The analysis reads only the DE tables, so
its results are the same as with the complete original file.

With the non-public raw data in the folders above, `prepare` also rebuilds the four DE
tables and `verify` compares them with the deposited ones.

## 5. Data preparation: the four input RData files

| File | Content | Script (`Scripts/data_preparation/`) |
|---|---|---|
| `Alldata_20Sep.RData` | human aging datasets (expression matrices, incl. `GTEx_PFC`) + `metadata` | `prep_01a_aging_expression.R`, `prep_01b_aging_metadata_Alldata_20Sep.R` |
| `Psychedelics_PFC.RData` | DE tables of psychedelics and antidepressants in PFC + LSD counts (MDMA, DMT, pharmahuasca, harmaline only with the non-public raw data, section 4) | `prep_02_Psychedelics_PFC.R` |
| `ExerciseAndAlcohol.RData` | DE tables of exercise, enriched environment and alcohol | `prep_03_ExerciseAndAlcohol.R` |
| `Dementia_alldata_meta.RData` | human dementia datasets + `metadata` | `prep_04_Dementia_alldata_meta.R` |

The source and processing of each dataset are described in the header of each
script. `GTEx_PFC` is built from the GTEx read counts: Ensembl IDs to HGNC symbols
(Ensembl 107; for duplicated symbols the gene with the highest mean count is kept);
log2(RPM + 1); removal of the 10% of samples with the highest mitochondrial read
counts (209 -> 188 samples); removal of genes with expression > 0 in fewer than 10
samples (30,888 genes).

The original files were produced in 2022 by interactive scripts (RStudio workspaces);
the scripts in `Scripts/data_preparation/` rewrite that code with the same
processing. The rebuilt files contain only the objects read by the analysis scripts.
The `verify` step compares them with the deposited files of `Data/original_RData_notclosed/`
(section 10). These lack only the raw MDMA and DMT tables, which are never rebuilt, so
every rebuilt object is compared.

## 6. Analysis scripts

`Scripts/Reproducibility/run_pipeline.R` runs the scripts in three layers; the scripts of a
layer run in parallel, each in its own R session (`Scripts/Reproducibility/run_script.R`),
with the R library given in `Data/envs/environments.csv` (section 8).

| Layer | Script | Main results |
|---|---|---|
| 0 | `00_LSD_DEGs.R` | DE of chronic LSD (`Results/RData/DE_GSE179379.RData`) |
| 0 | `00_GTEx_filt.R` | Fig S16 |
| 0 | `16_invitro_RNAseq.R` | Figs 3I, S13, S15 |
| 1 | `01_aging_metanalysis.R` | Figs S1, S3 |
| 1 | `03_Age_reversal.R` | Fig 1B |
| 1 | `04_integrated_dataset.R` | Figs 1A, S2, S4 |
| 1 | `05_DEvsAging_validation.R` | Fig 2B; Tables S15, S16 |
| 1 | `06_GO_aging.R` | gseGO of the aging datasets (input of Tables S17, S18) |
| 1 | `09_age_reversal_otherpsychedelics.R` (*) | Fig 1E; Tables S4, S7 |
| 1 | `11_validation_datasets.R` | Fig S11A |
| 1 | `12_dementia_cohen.R` (*) | Fig 1F; Tables S12, S13 |
| 1 | `13_AD_validation.R` | Fig S11B |
| 1 | `save_mean_logFC_selected_genes.R` (*) | Tables S8, S10 |
| 2 | `08_2_random_drugs.R` | Figs 1C-D, S7, S8 |
| 2 | `14_Signature_agreement.R` | Fig S6 |
| 2 | `17_AD_DEGs.R` | Fig S14 |
| 2 | `18_GO_allLSD.R` | Figs 2C-D; Tables S17, S18 |
| 2 | `19_GO_revert_dementia.R` | Figs 3J-K, S12; Tables S19, S20, S21 |
| 2 | `20_NES_figures.R` | Fig S5 |
| 2 | `21_cohen_figures.R` | Fig S10 |
| 2 | `Cumulative_logFC.R` | Fig S9; Table S9 |

(*) Without the non-public raw data, these scripts read the DE tables of MDMA, DMT,
pharmahuasca and harmaline from `Data/original_RData_notclosed/Psychedelics_PFC.RData`
(section 4).

Not run by the pipeline: `08_1_cMAP_data.R`, which prepares `Data/original_RData_notclosed/cmap/` from the
LINCS data and needs a personal clue.io API key.

The scripts start with `setwd("workdir")`: `run_script.R` turns it into a no-op and
runs each script in the run folder, which contains `Data/` (links) and `Results/`.

## 7. Intermediate files

**Written by the pipeline** (in the run folder):
- `rebuilt_RData/*.RData`: the four input RData files rebuilt by `prepare` (section 5),
  used only by `verify`;
- `Results/RData/`: objects passed between scripts, for example `DE_GSE179379.RData`
  (from `00_LSD_DEGs`), `fgsea_res*.RData` (GSEA results of `03`, `09`, `08_2`),
  `alldata_resid.RData`, `age_all.RData`, `dat_names.RData` (`04`), `gsea_GO*.RData`
  (`06`, `19`), `GTEx_RPM_PFC_filt.RData` (`00_GTEx_filt`);
- `Results/*.xlsx`: the tables.

**Distributed with the data** (`Data/original_RData_notclosed/`, unless stated):
- `cmap/`: the CMap input (section 3), copied to `Results/RData/` by the pipeline;
- `Data/intermediate_RData/`: outputs of the main pipeline read by
  `Threshold_sensitivity/` and `Threshold_sensitivity_invivo/`, which run from the
  repository root without the main pipeline (same content as the files of a full run):
  `DE_GSE179379.RData` (`00_LSD_DEGs.R`), `gsea_GO.RData` and `genes_cor.RData`
  (`06_GO_aging.R`);
- `seeded_intermediate/fgsea_res_drugs_RAW.RData`: the CMap GSEA of
  `08_2_random_drugs.R` with `set.seed(1234)` (section 9), as a reference copy;
- `clusterProfiler_cache_2022.rds`: the annotation caches of clusterProfiler used by
  `19_GO_revert_dementia.R` (section 8);
- `Data/original_RData_notclosed/rng_states.rds`: the random states read by the scripts
  (section 9);
- `Data/closed_data_intermediate/MSBB_gene_correlations.RData`: the MSBB gene correlations,
  read by `13_AD_validation.R` (section 4);
- the four input RData files of the paper (section 4).

## 8. R environments

R 4.5 with Bioconductor 3.21 software (`Data/packages/offline_repo_R4.5_windows`). Main
versions:

| Package | Version |
|---|---|
| fgsea | 1.34.2 |
| clusterProfiler | 4.16.0 |
| DESeq2 | 1.48.2 |
| limma | 3.64.3 |
| sva | 3.56.0 |
| GSVA | 2.2.1 |
| metap | 1.12 |
| ggplot2 | 4.0.2 |

Full lists: `Data/envs/*_versions.csv`.

The analyses were run at different times, with different versions of a few packages.
`install_environment.R` therefore builds one base library and a few small libraries
("overlays") that are put in front of it for the scripts that need them
(`Data/envs/environments.csv`):

| Library | Content | Scripts | Why |
|---|---|---|---|
| `base` | Bioconductor 3.21 software and annotations (`org.*.eg.db`, `GO.db` 3.21.0) | all others | |
| `prep` | `org.Mm.eg.db` 3.14.0, `mouse4302.db` 3.13.0, array-reading packages | data preparation | MDMA arrays (section 3) |
| `ann322` | `AnnotationDbi` 1.72.0 + `org.Hs/Mm/Rn.eg.db`, `GO.db` 3.22.0 | `save_mean_logFC_selected_genes.R`, `16_invitro_RNAseq.R` | gene annotation of Bioconductor 3.22 |
| `ann322_orgMm321` | as `ann322`, with `org.Mm.eg.db` 3.21.0 | `19_GO_revert_dementia.R` | as above; mouse annotation of 3.21 |
| `ggplot2_35` | `ggplot2` 3.5.2 + `ggrepel` 0.9.6 built against ggplot2 3.5 | `01_aging_metanalysis.R`, `05_DEvsAging_validation.R`, `12_dementia_cohen.R` | figures drawn with ggplot2 3.5 (theme, histogram bins) |
| `factominer212` | `FactoMineR` 2.12 | `04_integrated_dataset.R` | random numbers (section 9) |

Notes:
- The 3.22 annotation packages need `AnnotationDbi` >= 1.71.1, which is why the
  `ann322` libraries also contain `AnnotationDbi` 1.72.0.
- ggplot2 4.0 changes the default theme (axis text black instead of grey) and places
  the bins of `geom_histogram()` differently. The `ggrepel` 0.9.6 binary of the
  offline repository was built against ggplot2 4 and stops with ggplot2 3.5 ("could
  not find function gg_par"), so `ggplot2_35` also holds a `ggrepel` 0.9.6 built
  against ggplot2 3.5.
- `FactoMineR` 2.18 computes `PCA()` with `irlba`, which draws random numbers. In
  `04` the PCA runs between the random state set after `load()` and the GSEA of
  Fig 1A, so with 2.18 the GSEA p-values and NES change. The PCA coordinates are the
  same with both versions.
- **clusterProfiler caches.** clusterProfiler keeps the gene-GO table and the GO
  term table it builds in `.GlobalEnv` (`.GO_clusterProfiler_Env`, `.GOTERM_Env`,
  `.Anno_clusterProfiler_Env`) and reuses them instead of `org.Rn.eg.db` and
  `GO.db`. The original `ExerciseAndAlcohol.RData` workspace (21/09/2022) contained
  the three caches for rat (SYMBOL, BP), so the chronic-LSD `enrichGO` of `19` used
  the 2022 annotations. `19` restores them from `Data/original_RData_notclosed/clusterProfiler_cache_2022.rds`
  (extracted from that workspace) after its `load()` calls.

## 9. Random numbers

**Random states.** The original input RData files were RStudio workspace images that
contained `.Random.seed`, and each `load()` reset the generator to the state saved in
the file. These states are in `Data/original_RData_notclosed/rng_states.rds` (a list of four
626-integer Mersenne-Twister states; they cannot be written as `set.seed(n)`). The
scripts set them explicitly right after the `load()`:

```r
load("Data/Alldata_20Sep.RData")
assign(".Random.seed", readRDS("Data/original_RData_notclosed/rng_states.rds")$Alldata_20Sep, envir = .GlobalEnv)
```

- `$Alldata_20Sep`: `03`, `04`, `05`, `06`, `09`, `13`, `14`;
- `$Dementia_alldata_meta`: `17`, `19`.

The rebuilt input RData therefore do not need `.Random.seed`.

**Other fixed states.**
- `09`: the rat-aging GSEA (GSE75772) starts from the `Alldata_20Sep` state advanced
  by 1559 uniform draws (`runif(1559)`), the state of the run that produced it. The
  other interventions keep the state set after `load()`.
- `08_2_random_drugs.R`: `set.seed(1234)` right after the `load()`. The result is also
  distributed as `Data/original_RData_notclosed/seeded_intermediate/fgsea_res_drugs_RAW.RData`.
- The other scripts start with `set.seed(1234)`. Where a state is set after `load()`,
  that state is the one used.

**Rules for editing the scripts.**
- Do not add calls that draw random numbers between the `load()` and the line that
  sets the state, or between that line and the GSEA.
- Attaching `clusterProfiler` draws one random number: in `19` it is attached before
  the `load()` calls, in `14` after the line that sets the state.
- The results are reproduced exactly only with the package versions of section 8.
  A newer version of a package can draw random numbers where the old one did not, and
  then every later GSEA changes. Example: `FactoMineR` 2.18 draws random numbers in
  `PCA()` and changes the p-values and NES of Fig 1A; `04` therefore runs with
  `FactoMineR` 2.12 (library `factominer212`), which reproduces them.

**Graphics.** Some figures use `geom_jitter()`, `position_jitter()` or
`geom_text_repel()` without a seed (`05`, `09`, `12`, `16`, `20`, `21`, `Cumulative_logFC`): the positions of
the jittered points and repelled labels change from run to run; the plotted data do
not.

## 10. Numerical precision

The input RData rebuilt from the raw data agree with the original ones within
numerical precision:
- `Alldata_20Sep`: expression matrices within 7e-14 (RMA normalization and
  arithmetic of the 2025 package versions); `metadata` identical (row names differ);
- `Psychedelics_PFC`: within 1.4e-5 (DESeq2/limma; the DE tables of the non-public
  datasets were checked with the non-public data);
- `ExerciseAndAlcohol`: within 2.7e-4 (`DE_GSE111273_Y`, DESeq2 dispersion fit);
- `Dementia_alldata_meta`: expression matrices within 2e-15; `metadata` has 69 rows
  fewer (GSE110731, GSE53697, GSE64810: datasets with no expression matrix in the
  file, not used by the analysis).

The analysis always uses the deposited files, from which the results of the paper are
reproduced exactly; run on the rebuilt files, a few values would differ at the level of
the precision above.
