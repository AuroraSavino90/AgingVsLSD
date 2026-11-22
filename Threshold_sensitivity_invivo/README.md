# Threshold sensitivity, in vivo LSD (Fig. 2B–D, Tables S15–S18)

The paper text says that the genes of each acute LSD dose (mouse PFC, 50–500 µg/kg) were
selected with **FDR < 0.05**. The code (`Scripts/05_DEvsAging_validation.R`) selects them
with the **nominal p < 0.05** (`pvalue` column of DESeq2). The same gene sets feed the GO
analysis of `Scripts/18_GO_allLSD.R` (Fig. 2C–D, Tables S17–S18). Chronic LSD (GSE179379)
is selected with FDR < 0.05 in the code and in both modes here.

This folder re-runs 05 and 18 in two versions:

| Mode     | Gene selection, in vivo doses | Output folder     |
|----------|-------------------------------|-------------------|
| `pvalue` | nominal p < 0.05 (paper code) | `Results_pvalue/` |
| `fdr`    | BH FDR < 0.05 (`padj`)        | `Results_FDR/`    |

The original `Scripts/` and `Results/` are never modified. Run from the repository root
(`AgingVsLSD/`). The scripts read `Data/original_RData/Alldata_20Sep.RData`,
`Data/intermediate_RData/DE_GSE179379.RData`, `Data/raw_data_notclosed/` (`Invivo RNAseq/`, with
`Samples_LSD.xlsx`), `Data/homologs/`, `Data/ensembl/` and the main-pipeline outputs distributed in
`Data/intermediate_RData/`: `gsea_GO.RData`, `genes_cor.RData` (06). The main pipeline does
not need to be run first.

## How to run

Run it from the repository root, with the R libraries of `install_environment.R`
(05 with `ggplot2_35`, 18 with `base`, as in the main pipeline):

```
cd AgingVsLSD
REPRO_LIBS=C:/Users/<user>/AgingVsLSD_Rlib Rscript Threshold_sensitivity_invivo/run_all.R          # both modes
Rscript Threshold_sensitivity_invivo/run_all.R pvalue   # one mode
Rscript Threshold_sensitivity_invivo/run_all.R fdr
```

About 15 min per mode. Per-script logs are in `Results_*/logs/`. For now this analysis is
run with the complete RData; a version for the public release will be added once that
strategy is defined.

## Outputs

`Results_<mode>/Main_figures/`:

- `Fig2B_APR_invivo_<mode>.pdf`, `TableS16_APR_invivo_<mode>.xlsx`: APR per dose
- `TableS15_GSEA_invivo_<mode>.xlsx`: GSEA per aging dataset, dose and direction
- `Fig2C_left_GOup_LSD_<mode>.pdf`, `Fig2C_right_GOup_aging_<mode>.pdf`,
  `TableS17_GO_up_<mode>.xlsx`: GO terms up-regulated by LSD and down-regulated with aging
- `Fig2D_left_GOdn_LSD_<mode>.pdf`, `Fig2D_right_GOdn_aging_<mode>.pdf`,
  `TableS18_GO_dn_<mode>.xlsx`: GO terms down-regulated by LSD and up-regulated with aging

## Scripts (`Scripts/`)

Copies of the original scripts, run in this order by `run_all.R`:

1. `05_DEvsAging_validation.R`: in vivo DEGs, GSEA vs aging, APR, GO per dose
2. `18_GO_allLSD.R`: GO terms shared by the LSD regimens and aging

Differences from the originals (check with `diff Scripts/X.R Threshold_sensitivity_invivo/Scripts/X.R`):

- `setwd("workdir")` → `source("Threshold_sensitivity_invivo/Scripts/invivo_setup.R")`
- in vivo gene selection `$pvalue < 0.05` → `[[P_SEQ]] < 0.05`, set by `SELECTION_MODE`
- output paths `Results/...` and the working directory → `Results_<mode>/...`

`invivo_setup.R` also defines small wrappers (`fgsea`, `plotEnrichment`, `two2one`,
`sumlog`, `enrichGO`, `print`). With FDR the 50 µg/kg dose has no selected genes and
100 µg/kg only 3. The original code assumes non-empty gene sets, so the wrappers return NA
(GSEA) or an empty table (GO) for an empty set instead of stopping. Without empty sets they
call the original functions with the original arguments.

## Number of selected genes (up / down)

| Dose (µg/kg) | nominal p < 0.05 | FDR < 0.05 |
|---|---|---|
| 50  | 111 / 211 | 0 / 0 |
| 100 | 121 / 352 | 1 / 2 |
| 200 | 541 / 956 | 36 / 170 |
| 500 | 270 / 471 | 3 / 30 |

## Results with FDR < 0.05

Fig. 2B / Table S16, average APR:

| Dose (µg/kg) | nominal p | FDR |
|---|---|---|
| 500 | 22.1 | 1.0 |
| 200 | 264.8 | 63.6 |
| 100 | 55.6 | 0.5 |
| 50  | 57.8 | not computable (no genes) |
| mean over doses | 100.1 | 21.7 |

Only 200 µg/kg keeps a clear aging-opposition signal.

GO (Tables S17/S18):

- up-regulated by LSD: 163 terms instead of 184 (nominal p), 84 in common;
- down-regulated by LSD: 173 terms instead of 480 (nominal p), 153 in common.

Top 20 terms of Fig. 2C (up): out with FDR are learning or memory, memory, learning,
cognition, glutamate receptor signaling pathway, proton transmembrane transport, signal
release; in are vesicle localization, dendrite development, neurotransmitter secretion,
signal release from synapse, postsynapse organization, establishment of organelle
localization, positive regulation of synaptic transmission.

Top 20 terms of Fig. 2D (down): out with FDR are blood vessel morphogenesis, embryonic organ
development, wound healing, positive regulation of cell adhesion, glial cell development; in
are leukocyte migration, gland development, small molecule catabolic process, cell
chemotaxis, extracellular matrix organization.

## Known issue in the original code (kept here)

In `18_GO_allLSD.R` the GO results of the doses are in the order 500, 200, 100, 50
(`DEGs_seq` of 05) but the columns are labelled 50, 100, 200, 500. The dose columns of
Fig. 2C–D (left) and of Tables S17–S18 are therefore reversed, in the paper and here.
The corrected pipeline (`AgingVsLSD_corrected`) fixes the labels.

## Package versions

Run with R 4.5.2 (Windows) and the `bioc321` library (Bioconductor 3.21; full list in
`Data/envs/bioc321_versions.csv`), Ensembl release 115 for the mouse gene symbols.
