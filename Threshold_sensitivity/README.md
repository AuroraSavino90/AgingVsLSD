# Threshold sensitivity (APR and CRB, Fig. 1E–F)

In the original pipeline the intervention genes are selected with mixed thresholds:
**LSD (chronic, GSE179379) uses FDR < 0.05**, every other intervention dataset uses the
**nominal p < 0.05**. This folder re-runs the APR and CRB analyses with the **same
threshold for every dataset**, in two versions:

| Mode     | Gene selection, all datasets                         | Output folder      |
|----------|------------------------------------------------------|--------------------|
| `pvalue` | nominal p < 0.05 (`P.Value` limma, `pvalue` DESeq2)  | `Results_pvalue/`  |
| `fdr`    | BH FDR < 0.05 (`adj.P.Val` limma, `padj` DESeq2)     | `Results_FDR/`     |

The original `Scripts/` are never modified. The scripts here read the four complete input
RData in `Data/original_RData/`, `Data/intermediate_RData/DE_GSE179379.RData` (the
chronic LSD DE table of `00_LSD_DEGs.R`, distributed as an intermediate file),
`Data/raw_data_notclosed/` (GSE75772), `Data/homologs/`, `Data/ensembl/`, `Scripts/utilities.R`,
`Scripts/data_preparation/ensembl_release_map.R` and `Data/original_RData_notclosed/rng_states.rds`.
The main pipeline does not need to be run first.

## How to run

From the repository root (`AgingVsLSD/`), not from inside this folder, with the R
libraries of `install_environment.R`:

```
cd AgingVsLSD
REPRO_LIBS=C:/Users/<user>/AgingVsLSD_Rlib Rscript Threshold_sensitivity/run_all.R          # both modes
Rscript Threshold_sensitivity/run_all.R pvalue   # one mode
Rscript Threshold_sensitivity/run_all.R fdr
```

With `REPRO_LIBS` set, each script uses the library of the main-pipeline script it copies
(`Data/envs/environments.csv`).

About 1 h per mode with all data, mostly for the fgsea runs in `09`. Per-script logs are in `Results_*/logs/`.

**Data:** for now this analysis runs only with the complete RData, including the non-public
data (MDMA, DMT, pharmahuasca, harmaline; main README, section 4). `run_all.R` stops if
they are missing from `Data/original_RData/Psychedelics_PFC.RData`. A version for the public data only
will be added once the strategy for the public release is defined.

## Outputs

`Results_<mode>/Main_figures/`:

- `APR_aging_<mode>.pdf`, `APR_aging_scores_<mode>.xlsx`: APR vs human aging
- `APR_dementia_<mode>.pdf`, `APR_dementia_scores_<mode>.xlsx`: APR vs dementia
- `CRB_panel_E_aging_<mode>.pdf`, `CRB_panel_F_dementia_<mode>.pdf`,
  `CRB_values_panels_E_F_<mode>.xlsx`: cumulative reversal burden, Fig. 1E–F

All intermediate files (fgsea results, per-intervention GSEA plots, full tables) are in
`Results_<mode>/`, `Results_<mode>/RData/` and `Results_<mode>/Figures/`.

## Scripts (`Scripts/`)

Copies of the original scripts, run in this order by `run_all.R`:

1. `03_Age_reversal.R`: LSD (chronic) GSEA vs aging
2. `09_age_reversal_otherpsychedelics.R`: other interventions and rat aging, APR aging
3. `12_dementia_cohen.R`: APR dementia
4. `save_mean_logFC_selected_genes.R`: mean logFC tables (input to CRB)
5. `Cumulative_logFC.R`: CRB panels E–F

Differences from the originals (check with `diff Scripts/X.R Threshold_sensitivity/Scripts/X.R`):

- `setwd("workdir")` → `source("Threshold_sensitivity/Scripts/sensitivity_setup.R")`
- gene selection `$P.Value` / `$pvalue` / `$padj < 0.05` → `[[P_ARRAY]]` / `[[P_SEQ]]`,
  set by `SELECTION_MODE` for **all** datasets, LSD (chronic) included
- output paths `Results/...` → `Results_<mode>/...`
- `09`, rat aging (GSE75772):
  - Ensembl → symbol mapping from the pinned Ensembl release 115, as in `Scripts/09`
    (`Scripts/data_preparation/ensembl_release_map.R`, files in `Data/ensembl/`).
  - Before the rat-aging GSEA, the random state of the published run is set: the state of
    `Alldata_20Sep.RData` (`Data/original_RData_notclosed/rng_states.rds`) advanced by 1559
    `runif()` draws. `Scripts/09` uses `set.seed(1234)` here, which does not give the
    published NES/p. With this state the 15 rat-aging GSEA results are identical to
    `Results/RData/fgsea_res_GSE75772.RData` (size, ES, NES, p).

`sensitivity_setup.R` also defines small wrappers (`fgsea`, `plotEnrichment`, `two2one`,
`sumlog`, `gsva`, `t.test`, `print`). With FDR some interventions have no selected genes
(MDMA, DOI 2h, psilocybin 3h/4w, exercise 30d) or only 1–3 genes (LSD single,
alcohol IP, EE old). The
original code assumes non-empty UP and DN sets, so the wrappers return NA for an empty
set instead of stopping. Such interventions are missing from the FDR figures. With
non-empty sets the wrappers return exactly the original results.

## Package versions

The results in `Results_pvalue/` and `Results_FDR/` were produced with R 4.5.2 (Windows),
Bioconductor 3.22:

| Package | Version |
|---|---|
| fgsea | 1.36.2 |
| GSVA | 2.4.4 |
| limma | 3.66.0 |
| metap | 1.13 |
| AnnotationDbi | 1.72.0 |
| org.Rn.eg.db | 3.22.0 |
| data.table | 1.18.0 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| ggplot2 | 4.0.3 |
| ggrepel | 0.9.6 |
| pheatmap | 1.0.13 |
| openxlsx | 4.2.8.1 |
| readxl | 1.4.5 |

Rat aging gene symbols come from Ensembl release 115. fgsea 1.34 and 1.36 give identical
results from the same random state (see the main README, "Random numbers").
