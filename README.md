# AgingVsLSD

Code to reproduce the results described in the manuscript Title, Authors (link).

Input data can be downloaded at (link) and placed in the project folder.

- 00_GTEx_filt: filters GTEx data removing the 10% of samples with the highest expression of MT- genes (not equally distributed across age groups)
- 00_LSD_DEGs: computes the differentially expressed genes in the chronic LSD dataset
- 01_aging_metanalysis: produces distribution plots of metadata in the brain aging resource, identifies outlier datasets and finds genes coherently up- or down-regulated with aging
- 02_GO_LSD: computes the gene ontology enrichment for genes differentially expressed upon chronic LSD
- 03_Age_reversal: tests the transcriptional age reversal of chronic LSD for each human prefrontal cortex dataset separately
- 04_integrated_dataset: integrates human aging datasets of dorsolateral prefrontal cortex and tests the transcriptional age reversal of chronic LSD
- 05_DE_vsAging_validation: processes LSD validation data and tests their transcriptional age reversal potential
- 06_GO_aging: computes the gene ontology enrichment for genes correlated with age across the human brain aging datasets
- 07_GO_allLSD: computes the gene ontology enrichment for genes differentially expressed upon acute LSD administration
- 08_1_cMAP_data: downloads connectivity map data
- 08_2_random_drugs: tests the age-reversal potential of random drugs
- 09_age_reversal_otherpsychedelics: tests the age-reversal potential of other psychedelics or psychoplastogens with available post-perturbational transcriptomic data for each human prefrontal cortex dataset separately
- 10_GSEA_integrated_otherpsychedelics: tests the age-reversal potential of other psychedelics or psychoplastogens with available post-perturbational transcriptomic data using the integrated DLPFC aging dataset
- 11_validation_datasets: tests chronic LSD age-reversal potential in two aging datasets not included in the aging database
- 12_dementia_cohen: tests the potential psychedelics or psychoplastogens in reverting gene expression signatures of dementia
- 13_AD_validation: employs an external Alzheimer disease dataset to validate the age-reversal effect of chronic LSD
- 14_Signature_agreement: comparison of aging gene expression signatures obtained in the current study with the gene expression signature obtained by Donertas et al, 2018.
- 15_Gender_correct: tests chronic LSD age-reversal effect regressing out gender from gene expression
