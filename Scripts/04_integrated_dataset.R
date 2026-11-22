suppressMessages(library(DESeq2))  # DESeqResults objects need DESeq2/S4Vectors loaded (rownames() returns NULL otherwise)
rm(list=ls())
setwd("workdir")#workdir = working directory
set.seed(1234)  # fixed seed; where a random state is set after load() below, that state is used instead (README.md, "Random numbers")

library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
load("Data/original_RData_notclosed/Alldata_20Sep.RData")
assign(".Random.seed", readRDS("Reproducibility/rng_states.rds")$Alldata_20Sep, envir = .GlobalEnv)  # random state of the published runs (the .Random.seed saved in Data/original_RData_notclosed/Alldata_20Sep.RData), set explicitly
homologs<-read.csv("Data/homologs/Human rat homologs.txt")

region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE5388", "GSE102741"))#outlier

dat_names<-c()
age_all<-c()
dataset_all<-c()
for(dd in dat){

  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  age<-as.numeric(age)
  if(length(which(age>=7300))>=10){
    sample<-sample[which(age>=7300)]
    data<-data[, sample]
    age<-age[which(age>=7300)]
   dat_names<-c(dat_names, dd)
   age_all<-c(age_all, age)
   dataset_all<-c(dataset_all, rep(dd, ncol(data)))
  }
}

save(age_all, file="Results/RData/age_all.RData")


all<-c()
for(dd in dat_names){
  data<-get(dd)
  all<-c(all, rownames(data))
}

sum(table(all)>=11)
selgenes<-names(which(table(all)>=11))

############################ALLDATA
data<-get(dat_names[1])
sample<-metadata$Sample[which(metadata$Dataset==dat_names[1] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
data<-data[, sample]
age<-metadata$Age[which(metadata$Dataset==dat_names[1] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
age<-as.numeric(age)
if(length(which(age>=7300))>=10){
  sample<-sample[which(age>=7300)]
  data<-data[, sample]
alldata<-data[selgenes,]
}

for(i in 2:length(dat_names)){
  data<-get(dat_names[i])
  sample<-metadata$Sample[which(metadata$Dataset==dat_names[i] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dat_names[i] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  age<-as.numeric(age)
  if(length(which(age>=7300))>=10){
    sample<-sample[which(age>=7300)]
    data<-data[, sample]
    age<-age[which(age>=7300)]


  alldata<-cbind(alldata, data[selgenes,])
  }
}

library(FactoMineR)
pca<-PCA(t(alldata))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], age=age_all/365, dataset=dataset_all)

pdf("Results/Figures/PCA_DLPF_dataset.pdf",7,5)
ggplot(df, aes(PC2, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()
pdf("Results/Figures/PCA_DLPF_age.pdf",6,5)
ggplot(df, aes(PC2, PC1, colour=age))+geom_point(size=1)+ scale_color_gradient(low = "orange", high = "blue")+theme_classic()
dev.off()
pdf("Results/Figures/PC1_age.pdf",7,5)
ggplot(df, aes(age, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()

df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2],PC3=pca$ind$coord[,3],PC4=pca$ind$coord[,4], age=age_all/365, dataset=dataset_all)
ggplot(df, aes(PC1, PC4, colour=dataset))+geom_point()+theme_classic()
ggplot(df, aes(PC1, PC4, colour=age))+geom_point()+theme_classic()


batch = dataset_all

# parametric adjustment
library(sva)
mod = model.matrix(~age_all, data=data.frame(alldata))
mod0 = model.matrix(~1,data=data.frame(t(alldata)))
n.sv = 8  # the published alldata_resid used 8 SVs (num.sv(..., method="leek") now returns 28)
svobj = sva(as.matrix(alldata),mod,mod0,n.sv=n.sv)
plot(age_all, svobj[[1]][,1], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,2], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,3], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,4], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,5], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,6], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,7], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,8], col=as.factor(dataset_all))

alldata_resid<-apply(alldata,1, function(x){residuals(lm(unlist(x)~svobj[[1]][,1:n.sv]))})
#alldata_resid<-apply(alldata,1, function(x){residuals(lm(unlist(x)~svobj[[1]][,1]+svobj[[1]][,2]+svobj[[1]][,3]+svobj[[1]][,4]+svobj[[1]][,5]+svobj[[1]][,6]+svobj[[1]][,7]))})
alldata_resid<-t(alldata_resid)

pca<-PCA(t(alldata_resid))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2],PC3=pca$ind$coord[,3],PC4=pca$ind$coord[,4], age=age_all/365, dataset=dataset_all)
pdf("Results/Figures/PCA_SVA_dataset.pdf",7,5)
ggplot(df, aes(PC2, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()
pdf("Results/Figures/PCA_SVA_age.pdf",6,5)
ggplot(df, aes(PC2, PC1, colour=age))+geom_point()+ scale_color_gradient(low = "orange", high = "blue")+theme_classic()
dev.off()

ggplot(df, aes(PC3, PC4, colour=age))+geom_point()

pdf("Results/Figures/PC1_SVA_age.pdf",7,5)
ggplot(df, aes(age, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()
cor.test(df$age, df$PC1)

save(alldata_resid, file="Results/RData/alldata_resid.RData")

#################################################
#########GSEA with the batch corrected dataset
##################################################

library(fgsea)
library(ggplot2)

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

genes_cor<-cor(t(alldata_resid), as.numeric(age_all))

rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

  forgesea<-unlist(genes_cor)
  names(forgesea)<-rownames(genes_cor)
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_res<-fgsea(list(UP=genes_up, DN=genes_dn), sort(forgesea))
  
  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)
  
  df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data))), 
                       data.frame(p2$data, dir=rep("dn", nrow(p2$data))))
  
  p<-ggplot(df, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))

pdf("Results/Figures/GSEA_all_batchcorrect.pdf", 5,4)
ggplot(df, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))
dev.off()


#########################Integration performance

############################################################
### KNN prediction of dataset of origin
### Before vs after SVA batch correction
############################################################


library(class)
library(ggplot2)

set.seed(1234)


############################################################
### 1. Prepare data
############################################################

# samples x genes
X_raw <- t(alldata)
X_sva <- t(alldata_resid)

# dataset of origin
y <- factor(dataset_all)

# Basic checks
stopifnot(nrow(X_raw) == length(y))
stopifnot(nrow(X_sva) == length(y))
stopifnot(ncol(X_raw) == ncol(X_sva))

cat("Number of samples:", length(y), "\n")
cat("Number of datasets:", nlevels(y), "\n")
cat("Number of genes:", ncol(X_raw), "\n")

print(table(y))


############################################################
### 2. Stratified cross-validation folds
############################################################

make_stratified_folds <- function(y, K = 5) {
  
  folds <- vector("list", K)
  
  for (lev in levels(y)) {
    
    idx <- which(y == lev)
    idx <- sample(idx)
    
    fold_id <- rep(seq_len(K), length.out = length(idx))
    
    for (k in seq_len(K)) {
      folds[[k]] <- c(
        folds[[k]],
        idx[fold_id == k]
      )
    }
  }
  
  return(folds)
}


############################################################
### 3. Multiclass performance metrics
############################################################

multiclass_metrics <- function(obs, pred, all_levels) {
  
  obs <- factor(obs, levels = all_levels)
  pred <- factor(pred, levels = all_levels)
  
  cm <- table(obs, pred)
  
  # Overall accuracy
  accuracy <- sum(diag(cm)) / sum(cm)
  
  # Recall / sensitivity for each dataset
  recall <- diag(cm) / rowSums(cm)
  
  # Precision for each dataset
  precision <- diag(cm) / colSums(cm)
  
  # F1 for each dataset
  f1 <- 2 * precision * recall / (precision + recall)
  
  recall[!is.finite(recall)] <- NA
  precision[!is.finite(precision)] <- NA
  f1[!is.finite(f1)] <- NA
  
  balanced_accuracy <- mean(recall, na.rm = TRUE)
  macro_f1 <- mean(f1, na.rm = TRUE)
  
  return(
    c(
      Accuracy = accuracy,
      BalancedAccuracy = balanced_accuracy,
      MacroF1 = macro_f1
    )
  )
}


############################################################
### 4. KNN for one train/test split
###
### Important:
### PCA is fitted ONLY using training samples.
### Test samples are projected into the training PCA space.
############################################################

knn_one_split <- function(X,
                          y,
                          test_idx,
                          k = 5,
                          nPC = 20) {
  
  train_idx <- setdiff(
    seq_len(nrow(X)),
    test_idx
  )
  
  X_train <- X[train_idx, , drop = FALSE]
  X_test  <- X[test_idx, , drop = FALSE]
  
  y_train <- y[train_idx]
  
  ########################################################
  # Remove genes with zero / undefined variance
  # in the training set
  ########################################################
  
  sd_train <- apply(
    X_train,
    2,
    sd,
    na.rm = TRUE
  )
  
  keep <- which(
    is.finite(sd_train) &
      sd_train > 0
  )
  
  X_train <- X_train[, keep, drop = FALSE]
  X_test  <- X_test[, keep, drop = FALSE]
  
  
  ########################################################
  # PCA fitted only on training samples
  ########################################################
  
  npc_use <- min(
    nPC,
    nrow(X_train) - 1,
    ncol(X_train)
  )
  
  pca_train <- prcomp(
    X_train,
    center = TRUE,
    scale. = TRUE,
    rank. = npc_use
  )
  
  train_pc <- pca_train$x[
    ,
    seq_len(npc_use),
    drop = FALSE
  ]
  
  test_pc <- predict(
    pca_train,
    newdata = X_test
  )[
    ,
    seq_len(npc_use),
    drop = FALSE
  ]
  
  
  ########################################################
  # KNN
  ########################################################
  
  pred <- class::knn(
    train = train_pc,
    test = test_pc,
    cl = y_train,
    k = k
  )
  
  pred <- factor(
    pred,
    levels = levels(y)
  )
  
  return(pred)
}


############################################################
### 5. Repeated stratified cross-validation
###
### The SAME train/test splits are used for
### raw and SVA-corrected expression.
############################################################

run_knn_cv <- function(X_raw,
                       X_sva,
                       y,
                       K = 5,
                       n_repeats = 20,
                       k = 5,
                       nPC = 20,
                       seed = 1234) {
  
  set.seed(seed)
  
  perf <- list()
  predictions <- list()
  
  perf_i <- 0
  pred_i <- 0
  
  all_levels <- levels(y)
  
  
  for (r in seq_len(n_repeats)) {
    
    cat(
      "Repeat",
      r,
      "of",
      n_repeats,
      "\n"
    )
    
    # Generate folds once:
    # these are shared by raw and corrected data
    folds <- make_stratified_folds(
      y = y,
      K = K
    )
    
    
    for (f in seq_len(K)) {
      
      test_idx <- folds[[f]]
      
      
      ################################################
      # RAW DATA
      ################################################
      
      pred_raw <- knn_one_split(
        X = X_raw,
        y = y,
        test_idx = test_idx,
        k = k,
        nPC = nPC
      )
      
      met_raw <- multiclass_metrics(
        obs = y[test_idx],
        pred = pred_raw,
        all_levels = all_levels
      )
      
      
      ################################################
      # SVA-CORRECTED DATA
      ################################################
      
      pred_sva <- knn_one_split(
        X = X_sva,
        y = y,
        test_idx = test_idx,
        k = k,
        nPC = nPC
      )
      
      met_sva <- multiclass_metrics(
        obs = y[test_idx],
        pred = pred_sva,
        all_levels = all_levels
      )
      
      
      ################################################
      # Save fold-level performances
      ################################################
      
      perf_i <- perf_i + 1
      
      perf[[perf_i]] <- data.frame(
        rep_id = r,
        fold = f,
        data_type = "Before correction",
        Accuracy = unname(met_raw["Accuracy"]),
        BalancedAccuracy = unname(
          met_raw["BalancedAccuracy"]
        ),
        MacroF1 = unname(met_raw["MacroF1"])
      )
      
      
      perf_i <- perf_i + 1
      
      perf[[perf_i]] <- data.frame(
        rep_id = r,
        fold = f,
        data_type = "After SVA",
        Accuracy = unname(met_sva["Accuracy"]),
        BalancedAccuracy = unname(
          met_sva["BalancedAccuracy"]
        ),
        MacroF1 = unname(met_sva["MacroF1"])
      )
      
      
      ################################################
      # Save individual predictions
      ################################################
      
      pred_i <- pred_i + 1
      
      predictions[[pred_i]] <- data.frame(
        sample = rownames(X_raw)[test_idx],
        dataset = as.character(y[test_idx]),
        prediction = as.character(pred_raw),
        correct = as.character(pred_raw) ==
          as.character(y[test_idx]),
        rep_id = r,
        fold = f,
        data_type = "Before correction"
      )
      
      
      pred_i <- pred_i + 1
      
      predictions[[pred_i]] <- data.frame(
        sample = rownames(X_sva)[test_idx],
        dataset = as.character(y[test_idx]),
        prediction = as.character(pred_sva),
        correct = as.character(pred_sva) ==
          as.character(y[test_idx]),
        rep_id = r,
        fold = f,
        data_type = "After SVA"
      )
    }
  }
  
  
  perf <- do.call(
    rbind,
    perf
  )
  
  predictions <- do.call(
    rbind,
    predictions
  )
  
  
  return(
    list(
      performance = perf,
      predictions = predictions
    )
  )
}


############################################################
### 6. Run analysis
############################################################

res_knn <- run_knn_cv(
  X_raw = X_raw,
  X_sva = X_sva,
  y = y,
  K = 5,
  n_repeats = 20,
  k = 5,
  nPC = 20,
  seed = 1234
)

knn_perf <- res_knn$performance
knn_predictions <- res_knn$predictions


############################################################
### 7. Save results
############################################################

save(
  knn_perf,
  knn_predictions,
  file = "Results/RData/KNN_dataset_prediction.RData"
)

write.csv(
  knn_perf,
  "Results/KNN_dataset_prediction_performance.csv",
  row.names = FALSE
)

write.csv(
  knn_predictions,
  "Results/KNN_dataset_prediction_predictions.csv",
  row.names = FALSE
)


############################################################
### 8. Overall performance summary
############################################################

performance_summary <- aggregate(
  cbind(
    Accuracy,
    BalancedAccuracy,
    MacroF1
  ) ~ data_type,
  data = knn_perf,
  FUN = function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE)
    )
  }
)

print(performance_summary)


############################################################
### 9. Chance / null performance
############################################################

# Accuracy obtained by always predicting
# the largest dataset
majority_baseline <- max(
  prop.table(table(y))
)

# Expected balanced accuracy under random classification
chance_balanced_accuracy <- 1 / nlevels(y)

cat(
  "Majority-class accuracy baseline:",
  majority_baseline,
  "\n"
)

cat(
  "Balanced accuracy chance level:",
  chance_balanced_accuracy,
  "\n"
)


############################################################
### 10. Paired decrease in performance
############################################################

raw_perf <- knn_perf[
  knn_perf$data_type == "Before correction",
]

sva_perf <- knn_perf[
  knn_perf$data_type == "After SVA",
]

raw_perf <- raw_perf[
  order(
    raw_perf$rep_id,
    raw_perf$fold
  ),
]

sva_perf <- sva_perf[
  order(
    sva_perf$rep_id,
    sva_perf$fold
  ),
]

stopifnot(
  identical(
    raw_perf$rep_id,
    sva_perf$rep_id
  )
)

stopifnot(
  identical(
    raw_perf$fold,
    sva_perf$fold
  )
)


delta <- data.frame(
  rep_id = raw_perf$rep_id,
  fold = raw_perf$fold,
  
  delta_Accuracy =
    raw_perf$Accuracy -
    sva_perf$Accuracy,
  
  delta_BalancedAccuracy =
    raw_perf$BalancedAccuracy -
    sva_perf$BalancedAccuracy,
  
  delta_MacroF1 =
    raw_perf$MacroF1 -
    sva_perf$MacroF1
)

print(summary(delta))


cat(
  "\nMean reduction in accuracy:",
  mean(delta$delta_Accuracy),
  "\n"
)

cat(
  "Mean reduction in balanced accuracy:",
  mean(delta$delta_BalancedAccuracy),
  "\n"
)

cat(
  "Mean reduction in Macro-F1:",
  mean(delta$delta_MacroF1),
  "\n"
)


############################################################
### 11. Better statistical unit:
### average folds within each CV repeat
############################################################

perf_repeat <- aggregate(
  cbind(
    Accuracy,
    BalancedAccuracy,
    MacroF1
  ) ~ rep_id + data_type,
  data = knn_perf,
  FUN = mean
)

raw_repeat <- perf_repeat[
  perf_repeat$data_type == "Before correction",
]

sva_repeat <- perf_repeat[
  perf_repeat$data_type == "After SVA",
]

raw_repeat <- raw_repeat[
  order(raw_repeat$rep_id),
]

sva_repeat <- sva_repeat[
  order(sva_repeat$rep_id),
]


############################################################
### 12. Paired statistical tests
############################################################

test_accuracy <- wilcox.test(
  raw_repeat$Accuracy,
  sva_repeat$Accuracy,
  paired = TRUE,
  exact = FALSE
)

test_balanced <- wilcox.test(
  raw_repeat$BalancedAccuracy,
  sva_repeat$BalancedAccuracy,
  paired = TRUE,
  exact = FALSE
)

test_macroF1 <- wilcox.test(
  raw_repeat$MacroF1,
  sva_repeat$MacroF1,
  paired = TRUE,
  exact = FALSE
)

print(test_accuracy)
print(test_balanced)
print(test_macroF1)


############################################################
### 13. Main figure:
### dataset prediction balanced accuracy
############################################################

pdf(
  "Results/Figures/KNN_dataset_prediction_before_after_SVA.pdf",
  width = 4,
  height = 4
)

print(
  ggplot(
    perf_repeat,
    aes(
      x = data_type,
      y = BalancedAccuracy
    )
  ) +
    geom_hline(
      yintercept = chance_balanced_accuracy,
      linetype = "dashed"
    ) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA
    ) +
    geom_point(
      position = position_jitter(
        width = 0.07,
        height = 0
      ),
      size = 1.5,
      alpha = 0.5
    ) +
    theme_classic() +
    labs(
      x = NULL,
      y = "Dataset prediction\nbalanced accuracy"
    )
)

dev.off()


############################################################
### 14. Paired representation
############################################################

pdf(
  "Results/Figures/KNN_dataset_prediction_paired.pdf",
  width = 4,
  height = 4
)

print(
  ggplot(
    perf_repeat,
    aes(
      x = data_type,
      y = BalancedAccuracy,
      group = rep_id
    )
  ) +
    geom_hline(
      yintercept = chance_balanced_accuracy,
      linetype = "dashed"
    ) +
    geom_line(
      alpha = 0.25
    ) +
    geom_point(
      size = 1.5
    ) +
    theme_classic() +
    labs(
      x = NULL,
      y = "Dataset prediction\nbalanced accuracy"
    )
)

dev.off()


############################################################
### 15. Dataset-specific accuracy
############################################################

dataset_perf <- aggregate(
  correct ~ dataset + data_type,
  data = knn_predictions,
  FUN = mean
)

colnames(dataset_perf)[3] <- "Accuracy"

print(dataset_perf)

write.csv(
  dataset_perf,
  "Results/KNN_dataset_specific_accuracy.csv",
  row.names = FALSE
)


############################################################
### 16. Dataset-specific figure
############################################################

pdf(
  "Results/Figures/KNN_dataset_specific_accuracy.pdf",
  width = 6,
  height = max(
    5,
    0.3 * nlevels(y)
  )
)

print(
  ggplot(
    dataset_perf,
    aes(
      x = Accuracy,
      y = reorder(
        dataset,
        Accuracy
      ),
      shape = data_type
    )
  ) +
    geom_point(
      size = 2.5
    ) +
    theme_classic() +
    labs(
      x = "Dataset-specific prediction accuracy",
      y = NULL,
      shape = NULL
    )
)

dev.off()


############################################################
### 17. Confusion matrices
### aggregated across repeated CV
############################################################

cm_raw <- table(
  True = knn_predictions$dataset[
    knn_predictions$data_type == "Before correction"
  ],
  Predicted = knn_predictions$prediction[
    knn_predictions$data_type == "Before correction"
  ]
)

cm_sva <- table(
  True = knn_predictions$dataset[
    knn_predictions$data_type == "After SVA"
  ],
  Predicted = knn_predictions$prediction[
    knn_predictions$data_type == "After SVA"
  ]
)

print(cm_raw)
print(cm_sva)

write.csv(
  cm_raw,
  "Results/KNN_confusion_matrix_before_correction.csv"
)

write.csv(
  cm_sva,
  "Results/KNN_confusion_matrix_after_SVA.csv"
)


############################################################
### 18. Final concise summary
############################################################

cat("\n---------------------------------------------\n")
cat("KNN DATASET-OF-ORIGIN CLASSIFICATION\n")
cat("---------------------------------------------\n")

cat(
  "Balanced accuracy before correction:",
  mean(
    raw_repeat$BalancedAccuracy,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Balanced accuracy after SVA:",
  mean(
    sva_repeat$BalancedAccuracy,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Absolute reduction:",
  mean(
    raw_repeat$BalancedAccuracy -
      sva_repeat$BalancedAccuracy,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Relative reduction (%):",
  100 *
    mean(
      raw_repeat$BalancedAccuracy -
        sva_repeat$BalancedAccuracy,
      na.rm = TRUE
    ) /
    mean(
      raw_repeat$BalancedAccuracy,
      na.rm = TRUE
    ),
  "\n"
)

cat(
  "Chance balanced accuracy:",
  chance_balanced_accuracy,
  "\n"
)

cat(
  "Paired Wilcoxon p-value:",
  test_balanced$p.value,
  "\n"
)

cat("---------------------------------------------\n")


perf_repeat <- aggregate(
  cbind(
    Accuracy,
    BalancedAccuracy,
    MacroF1
  ) ~ rep_id + data_type,
  data = knn_perf,
  FUN = mean
)

ggplot(
  perf_repeat,
  aes(x = data_type, y = BalancedAccuracy)
) +
  geom_hline(
    yintercept = 1 / nlevels(y),
    linetype = "dashed"
  ) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitter(width = 0.06),
    alpha = 0.6,
    size = 1.5
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Dataset prediction\nbalanced accuracy"
  )

ggplot(
  perf_repeat,
  aes(
    x = data_type,
    y = BalancedAccuracy,
    group = rep_id
  )
) +
  geom_hline(
    yintercept = 1 / nlevels(y),
    linetype = "dashed"
  ) +
  geom_line(alpha = 0.3) +
  geom_point(size = 1.7) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Dataset prediction\nbalanced accuracy"
  )


############################################################
### Variance explained by Dataset vs Age
### Before and after SVA correction
###
### No variancePartition package required
### Base R + ggplot2 only
############################################################



############################################################
### 1. Prepare expression matrices
############################################################

expr_raw <- as.matrix(alldata)
expr_sva <- as.matrix(alldata_resid)

# Keep exactly the same genes
common_genes <- intersect(
  rownames(expr_raw),
  rownames(expr_sva)
)

expr_raw <- expr_raw[
  common_genes,
  ,
  drop = FALSE
]

expr_sva <- expr_sva[
  common_genes,
  ,
  drop = FALSE
]


# Check/reorder samples if necessary
if (!identical(
  colnames(expr_raw),
  colnames(expr_sva)
)) {
  
  if (setequal(
    colnames(expr_raw),
    colnames(expr_sva)
  )) {
    
    expr_sva <- expr_sva[
      ,
      colnames(expr_raw),
      drop = FALSE
    ]
    
  } else {
    
    stop(
      "Raw and SVA matrices do not contain the same samples."
    )
  }
}


############################################################
### 2. Metadata
############################################################

stopifnot(
  ncol(expr_raw) == length(age_all)
)

stopifnot(
  ncol(expr_raw) == length(dataset_all)
)


meta <- data.frame(
  Age = as.numeric(age_all) / 365,
  Dataset = factor(dataset_all)
)

rownames(meta) <- colnames(expr_raw)

# Standardization is not needed for R2,
# but makes the age coefficient more stable/readable
meta$Age_scaled <- as.numeric(
  scale(meta$Age)
)


cat(
  "Samples:",
  nrow(meta),
  "\n"
)

cat(
  "Genes:",
  nrow(expr_raw),
  "\n"
)

cat(
  "Datasets:",
  nlevels(meta$Dataset),
  "\n\n"
)

print(table(meta$Dataset))


############################################################
### 3. Model matrices
############################################################

# Null model
X_null <- model.matrix(
  ~ 1,
  data = meta
)

# Age only
X_age <- model.matrix(
  ~ Age_scaled,
  data = meta
)

# Dataset only
X_dataset <- model.matrix(
  ~ Dataset,
  data = meta
)

# Full model
X_full <- model.matrix(
  ~ Age_scaled + Dataset,
  data = meta
)


############################################################
### 4. Function to calculate RSS for ALL genes simultaneously
###
### This is much faster than fitting lm() separately
### for thousands of genes.
############################################################

get_rss <- function(expr, design) {
  
  # lm.fit expects:
  # rows = samples
  # columns of response = genes
  
  Y <- t(expr)
  
  fit <- lm.fit(
    x = design,
    y = Y
  )
  
  residuals <- fit$residuals
  
  if (is.null(dim(residuals))) {
    
    rss <- sum(
      residuals^2,
      na.rm = TRUE
    )
    
  } else {
    
    rss <- colSums(
      residuals^2,
      na.rm = TRUE
    )
  }
  
  return(rss)
}


############################################################
### 5. Function for variance decomposition
############################################################

variance_decomposition <- function(expr,
                                   X_null,
                                   X_age,
                                   X_dataset,
                                   X_full) {
  
  cat("  Null model...\n")
  
  rss_null <- get_rss(
    expr,
    X_null
  )
  
  
  cat("  Age model...\n")
  
  rss_age <- get_rss(
    expr,
    X_age
  )
  
  
  cat("  Dataset model...\n")
  
  rss_dataset <- get_rss(
    expr,
    X_dataset
  )
  
  
  cat("  Full model...\n")
  
  rss_full <- get_rss(
    expr,
    X_full
  )
  
  
  ########################################################
  ### Total variance
  ########################################################
  
  tss <- rss_null
  
  
  ########################################################
  ### Full model R2
  ########################################################
  
  full_R2 <- 1 - (
    rss_full / tss
  )
  
  
  ########################################################
  ### UNIQUE effect of AGE
  ###
  ### Compare:
  ### Dataset
  ### vs
  ### Dataset + Age
  ########################################################
  
  age_unique <- (
    rss_dataset -
      rss_full
  ) / tss
  
  
  ########################################################
  ### UNIQUE effect of DATASET
  ###
  ### Compare:
  ### Age
  ### vs
  ### Age + Dataset
  ########################################################
  
  dataset_unique <- (
    rss_age -
      rss_full
  ) / tss
  
  
  ########################################################
  ### Residual variance
  ########################################################
  
  residual <- rss_full / tss
  
  
  ########################################################
  ### Shared Age/Dataset component
  ###
  ### This captures variance that cannot be uniquely
  ### assigned because Age and Dataset are associated.
  ########################################################
  
  shared <- full_R2 -
    age_unique -
    dataset_unique
  
  
  ########################################################
  ### Numerical cleanup
  ########################################################
  
  age_unique[
    age_unique < 0 &
      age_unique > -1e-10
  ] <- 0
  
  dataset_unique[
    dataset_unique < 0 &
      dataset_unique > -1e-10
  ] <- 0
  
  
  ########################################################
  ### Output
  ########################################################
  
  result <- data.frame(
    Gene = rownames(expr),
    Age_unique = age_unique,
    Dataset_unique = dataset_unique,
    Shared_Age_Dataset = shared,
    Residual = residual,
    Full_R2 = full_R2,
    TSS = tss,
    RSS_full = rss_full
  )
  
  return(result)
}


############################################################
### 6. RAW DATA
############################################################

cat(
  "\n====================================\n",
  "RAW EXPRESSION\n",
  "====================================\n"
)

vp_raw <- variance_decomposition(
  expr = expr_raw,
  X_null = X_null,
  X_age = X_age,
  X_dataset = X_dataset,
  X_full = X_full
)


############################################################
### 7. SVA-CORRECTED DATA
############################################################

cat(
  "\n====================================\n",
  "SVA-CORRECTED EXPRESSION\n",
  "====================================\n"
)

vp_sva <- variance_decomposition(
  expr = expr_sva,
  X_null = X_null,
  X_age = X_age,
  X_dataset = X_dataset,
  X_full = X_full
)


############################################################
### 8. Basic checks
############################################################

summary(vp_raw$Dataset_unique)
summary(vp_sva$Dataset_unique)

summary(vp_raw$Age_unique)
summary(vp_sva$Age_unique)


############################################################
### 9. Gene-wise summary
############################################################

summary_gene <- data.frame(
  
  Component = c(
    "Dataset",
    "Age"
  ),
  
  Before_mean = c(
    mean(
      vp_raw$Dataset_unique,
      na.rm = TRUE
    ),
    mean(
      vp_raw$Age_unique,
      na.rm = TRUE
    )
  ),
  
  After_mean = c(
    mean(
      vp_sva$Dataset_unique,
      na.rm = TRUE
    ),
    mean(
      vp_sva$Age_unique,
      na.rm = TRUE
    )
  ),
  
  Before_median = c(
    median(
      vp_raw$Dataset_unique,
      na.rm = TRUE
    ),
    median(
      vp_raw$Age_unique,
      na.rm = TRUE
    )
  ),
  
  After_median = c(
    median(
      vp_sva$Dataset_unique,
      na.rm = TRUE
    ),
    median(
      vp_sva$Age_unique,
      na.rm = TRUE
    )
  )
)


summary_gene$Absolute_change_median <-
  summary_gene$After_median -
  summary_gene$Before_median


summary_gene$Relative_change_percent <-
  100 *
  (
    summary_gene$After_median -
      summary_gene$Before_median
  ) /
  summary_gene$Before_median


cat(
  "\nGENE-WISE SUMMARY\n"
)

print(summary_gene)


############################################################
### 10. Transcriptome-wide GLOBAL variance decomposition
###
### Instead of averaging genes equally, here all
### squared deviations across genes and samples are pooled.
############################################################

global_decomposition <- function(vp) {
  
  total_tss <- sum(
    vp$TSS,
    na.rm = TRUE
  )
  
  residual_ss <- sum(
    vp$RSS_full,
    na.rm = TRUE
  )
  
  age_ss <- sum(
    vp$Age_unique *
      vp$TSS,
    na.rm = TRUE
  )
  
  dataset_ss <- sum(
    vp$Dataset_unique *
      vp$TSS,
    na.rm = TRUE
  )
  
  full_explained_ss <-
    total_tss -
    residual_ss
  
  shared_ss <-
    full_explained_ss -
    age_ss -
    dataset_ss
  
  
  data.frame(
    Age_unique =
      age_ss / total_tss,
    
    Dataset_unique =
      dataset_ss / total_tss,
    
    Shared_Age_Dataset =
      shared_ss / total_tss,
    
    Residual =
      residual_ss / total_tss,
    
    Full_R2 =
      full_explained_ss / total_tss
  )
}


global_raw <- global_decomposition(
  vp_raw
)

global_sva <- global_decomposition(
  vp_sva
)


global_summary <- rbind(
  
  data.frame(
    Correction = "Before correction",
    global_raw
  ),
  
  data.frame(
    Correction = "After SVA",
    global_sva
  )
)


cat(
  "\n====================================\n",
  "GLOBAL TRANSCRIPTOME-WIDE RESULTS\n",
  "====================================\n"
)

print(global_summary)


############################################################
### 11. Percentage version for easier interpretation
############################################################

global_percent <- global_summary

num_cols <- setdiff(
  colnames(global_percent),
  "Correction"
)

global_percent[
  ,
  num_cols
] <- global_percent[
  ,
  num_cols
] * 100


cat(
  "\nGLOBAL VARIANCE (%):\n"
)

print(global_percent)


############################################################
### 12. Explicit before/after changes
############################################################

dataset_before <-
  global_summary$Dataset_unique[
    global_summary$Correction ==
      "Before correction"
  ]

dataset_after <-
  global_summary$Dataset_unique[
    global_summary$Correction ==
      "After SVA"
  ]


age_before <-
  global_summary$Age_unique[
    global_summary$Correction ==
      "Before correction"
  ]

age_after <-
  global_summary$Age_unique[
    global_summary$Correction ==
      "After SVA"
  ]


dataset_reduction <-
  dataset_before -
  dataset_after


dataset_reduction_percent <-
  100 *
  dataset_reduction /
  dataset_before


age_change <-
  age_after -
  age_before


age_change_percent <-
  100 *
  age_change /
  age_before


cat(
  "\n----------------------------------------\n"
)

cat(
  "Dataset unique variance BEFORE:",
  round(
    100 * dataset_before,
    3
  ),
  "%\n"
)

cat(
  "Dataset unique variance AFTER:",
  round(
    100 * dataset_after,
    3
  ),
  "%\n"
)

cat(
  "Dataset variance reduction:",
  round(
    100 * dataset_reduction,
    3
  ),
  "percentage points\n"
)

cat(
  "Relative dataset variance reduction:",
  round(
    dataset_reduction_percent,
    1
  ),
  "%\n\n"
)


cat(
  "Age unique variance BEFORE:",
  round(
    100 * age_before,
    3
  ),
  "%\n"
)

cat(
  "Age unique variance AFTER:",
  round(
    100 * age_after,
    3
  ),
  "%\n"
)

cat(
  "Age variance change:",
  round(
    100 * age_change,
    3
  ),
  "percentage points\n"
)

cat(
  "Relative age variance change:",
  round(
    age_change_percent,
    1
  ),
  "%\n"
)

cat(
  "----------------------------------------\n"
)


############################################################
### 13. Save results
############################################################

write.csv(
  vp_raw,
  "Results/variance_components_per_gene_before_SVA.csv",
  row.names = FALSE
)

write.csv(
  vp_sva,
  "Results/variance_components_per_gene_after_SVA.csv",
  row.names = FALSE
)

write.csv(
  summary_gene,
  "Results/variance_components_gene_summary.csv",
  row.names = FALSE
)

write.csv(
  global_summary,
  "Results/variance_components_global_summary.csv",
  row.names = FALSE
)

save(
  vp_raw,
  vp_sva,
  summary_gene,
  global_summary,
  file =
    "Results/RData/variance_partition_age_dataset.RData"
)


############################################################
### 14. Main plotting dataframe
############################################################

plot_df <- rbind(
  
  data.frame(
    Gene = vp_raw$Gene,
    Correction = "Before correction",
    Component = "Dataset",
    Variance = vp_raw$Dataset_unique
  ),
  
  data.frame(
    Gene = vp_sva$Gene,
    Correction = "After SVA",
    Component = "Dataset",
    Variance = vp_sva$Dataset_unique
  ),
  
  data.frame(
    Gene = vp_raw$Gene,
    Correction = "Before correction",
    Component = "Age",
    Variance = vp_raw$Age_unique
  ),
  
  data.frame(
    Gene = vp_sva$Gene,
    Correction = "After SVA",
    Component = "Age",
    Variance = vp_sva$Age_unique
  )
)


plot_df$Correction <- factor(
  plot_df$Correction,
  levels = c(
    "Before correction",
    "After SVA"
  )
)

plot_df$Component <- factor(
  plot_df$Component,
  levels = c(
    "Dataset",
    "Age"
  )
)


############################################################
### 15. Main figure:
### gene-wise variance distributions
############################################################

pdf(
  "Results/Figures/variance_dataset_age_before_after_SVA.pdf",
  width = 7,
  height = 4
)

print(
  ggplot(
    plot_df,
    aes(
      x = Correction,
      y = Variance * 100
    )
  ) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA
    ) +
    facet_wrap(
      ~ Component,
      scales = "free_y"
    ) +
    theme_classic() +
    labs(
      x = NULL,
      y = "Unique variance explained (%)"
    )
)

dev.off()


############################################################
### 16. Dataset variance:
### each gene before vs after SVA
############################################################

dataset_scatter <- data.frame(
  Gene = vp_raw$Gene,
  Before = vp_raw$Dataset_unique,
  After = vp_sva$Dataset_unique
)


pdf(
  "Results/Figures/variance_dataset_gene_before_vs_after_SVA.pdf",
  width = 5,
  height = 5
)

print(
  ggplot(
    dataset_scatter,
    aes(
      x = Before * 100,
      y = After * 100
    )
  ) +
    geom_point(
      size = 0.5,
      alpha = 0.15
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    theme_classic() +
    labs(
      x =
        "Dataset-associated variance before SVA (%)",
      y =
        "Dataset-associated variance after SVA (%)"
    )
)

dev.off()


############################################################
### 17. Age variance:
### each gene before vs after SVA
############################################################

age_scatter <- data.frame(
  Gene = vp_raw$Gene,
  Before = vp_raw$Age_unique,
  After = vp_sva$Age_unique
)


pdf(
  "Results/Figures/variance_age_gene_before_vs_after_SVA.pdf",
  width = 5,
  height = 5
)

print(
  ggplot(
    age_scatter,
    aes(
      x = Before * 100,
      y = After * 100
    )
  ) +
    geom_point(
      size = 0.5,
      alpha = 0.15
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    theme_classic() +
    labs(
      x =
        "Age-associated variance before SVA (%)",
      y =
        "Age-associated variance after SVA (%)"
    )
)

dev.off()


############################################################
### 18. Compact global figure
############################################################

global_plot <- rbind(
  
  data.frame(
    Correction = "Before correction",
    Component = "Dataset",
    Variance =
      global_raw$Dataset_unique * 100
  ),
  
  data.frame(
    Correction = "After SVA",
    Component = "Dataset",
    Variance =
      global_sva$Dataset_unique * 100
  ),
  
  data.frame(
    Correction = "Before correction",
    Component = "Age",
    Variance =
      global_raw$Age_unique * 100
  ),
  
  data.frame(
    Correction = "After SVA",
    Component = "Age",
    Variance =
      global_sva$Age_unique * 100
  )
)


global_plot$Correction <- factor(
  global_plot$Correction,
  levels = c(
    "Before correction",
    "After SVA"
  )
)


pdf(
  "Results/Figures/global_variance_dataset_age_before_after_SVA.pdf",
  width = 5,
  height = 4
)

print(
  ggplot(
    global_plot,
    aes(
      x = Correction,
      y = Variance,
      group = Component,
      shape = Component
    )
  ) +
    geom_line() +
    geom_point(
      size = 3
    ) +
    theme_classic() +
    labs(
      x = NULL,
      y = "Transcriptome-wide variance explained (%)",
      shape = NULL
    )
)

dev.off()

