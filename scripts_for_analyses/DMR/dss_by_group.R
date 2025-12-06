library(bsseq)
library(DSS)
######################
#Rscript dss_human_celltype.R number_of_randomly_chosen_samples

## take input RDS file as first argument
args <- commandArgs(trailingOnly=TRUE)
rds_file <- args[1]
bs.combine <- readRDS(rds_file)

## take input covariates file as second argument
#covariates_file <- args[2]
#covariates_R <- read.table(covariates_file,header=T)
#row.names(covariates_R) <- covariates_R$Sample

sample_names <- sampleNames(bs.combine)
#covariates_R_subset <- covariates_R[sample_names,]
#covariates_R_subset <- droplevels(covariates_R_subset)
#covariates_R_subset$State <- factor(covariates_R_subset$State, levels=c("Ref","CKD"))

#ind.REF <- rownames(covariates_R_subset[covariates_R_subset$State == "control",])
#ind.CKD <- rownames(covariates_R_subset[covariates_R_subset$State == "disease",])

min.cov <- 5
#ind.num <- 2
bs.cov <- getCoverage(bs.combine,type="Cov")
#bs.good.loci <- which(rowSums(bs.cov[,ind.REF] >= min.cov) >= round(length(ind.neun)*ind.per,0) & rowSums(bs.cov[,ind.olig] >= min.cov) >= round(length(ind.olig)*ind.per,0))
bs.good.loci <- which(bs.cov[,"control"] >= min.cov  & bs.cov[,"disease"] >= min.cov)
bs.combine.filtered <- bs.combine[bs.good.loci,]
bs.combine.filtered <- sort(bs.combine.filtered)
bs.combined.filtered.REF.CKD <- bs.combine.filtered[,c("control","disease")]


##########################
#GT_ID   UTSW_ID Species Sex     Age_Class       Pmi     Conversion_rates
#X3590_Control_NeuN      X3590_Control_NeuN      Control M       3       11.5    0.9963784915977140
#pData(bs.combine.filtered)$CellType <- covariates_R_order$CellType
#pData(bs.combine.filtered)$Species <- covariates_R_order$Species
#pData(bs.combine.filtered)$Sex <- covariates_R_order$Sex
#pData(bs.combine.filtered)$Age_Class <- covariates_R_order$Age_Class
#pData(bs.combine.filtered)$Conversion_rates <- covariates_R_order$Conversion_rates
#pData(bs.combine.filtered)<- droplevels(pData(bs.combine.filtered))

###################

#model.matrix(~Species+CellType+Sex+Age_Class+Conversion_rates+Species:CellType,pData(bs.combine.filtered))
#model.matrix(~Celltype + Age + Gender + Disease_level, covariates_R_subset)

dmlTest.sm = DMLtest(bs.combined.filtered.REF.CKD, group1="control", group2="disease",smoothing=TRUE, ncores=4)
dmls = callDML(dmlTest.sm, delta=0.05)
dmrs = callDMR(dmlTest.sm, delta=0.05, p.threshold=0.005,minCG=3)
dmls$pos <- format(dmls$pos, scientific = FALSE)
dmrs$start <- format(dmrs$start, scientific = FALSE)
dmrs$end <- format(dmrs$end, scientific = FALSE)
write.table(dmls,file=args[2],sep='\t',row.names=F,quote=F)
write.table(dmrs,file=args[3],sep='\t',row.names=F,quote=F)
