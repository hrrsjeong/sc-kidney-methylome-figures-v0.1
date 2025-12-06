library(bsseq)
library(DSS)
######################
#Rscript dss_human_celltype.R number_of_randomly_chosen_samples

## take input RDS file as first argument
args <- commandArgs(trailingOnly=TRUE)
##print args
print(args)
print(args[1])
print(args[2])
print(args[3])
print(args[4])

rds_file <- args[1]
bs.combine <- readRDS(rds_file)

## take input covariates file as second argument
covariates_file <- args[2]
covariates_R <- read.table(covariates_file,header=T)
row.names(covariates_R) <- covariates_R$Sample

sample_names <- sampleNames(bs.combine)
covariates_R_subset <- covariates_R[sample_names,]
covariates_R_subset <- droplevels(covariates_R_subset)
#covariates_R_subset$State <- factor(covariates_R_subset$State, levels=c("control","disease"))
covariates_R_subset$Cell_state <- factor(covariates_R_subset$Cell_state, levels=c("healthy","altered"))

ind.REF <- rownames(covariates_R_subset[grep("healthy",covariates_R_subset$Cell_state),])
ind.CKD <- rownames(covariates_R_subset[grep("altered",covariates_R_subset$Cell_state),])

min.cov <- 2
ind.num <- 6
bs.cov <- getCoverage(bs.combine,type="Cov")
#bs.good.loci <- which(rowSums(bs.cov[,ind.REF] >= min.cov) >= round(length(ind.neun)*ind.per,0) & rowSums(bs.cov[,ind.olig] >= min.cov) >= round(length(ind.olig)*ind.per,0))
bs.good.loci <- which(rowSums(bs.cov[,ind.REF] >= min.cov) >= 3 & rowSums(bs.cov[,ind.CKD] >= min.cov) >= 3)
bs.combine.filtered <- bs.combine[bs.good.loci,]
bs.combine.filtered <- sort(bs.combine.filtered)
bs.combined.filtered.REF.CKD <- bs.combine.filtered[,c(ind.REF,ind.CKD)]


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

DMLfit = DMLfit.multiFactor(bs.combine.filtered,covariates_R_subset, ~Cell_state+Sample_id)
#saveRDS(DMLfit,"dss_human_celltype_DMLfit_own_ref_20.rds")
#DMLfit <- readRDS("dss_human_celltype_DMLfit_own_ref.rds")
DMLtest_state = DMLtest.multiFactor(DMLfit, coef="Cell_statealtered")
#fdr_max <- DMLtest_Age[DMLtest_Age$fdrs<0.1,][which.max(DMLtest_Age[DMLtest_Age$fdrs<0.2,]$pvals),]
#ix=sort.int(DMLtest_Age[,"pvals"], index.return=TRUE,na.last=TRUE)$ix
write.table(DMLtest_state,file=args[3],sep="\t",row.names=F,quote=F)
DMRtest_state <- callDMR(DMLtest_state, p.threshold=0.05,minCG=3)
write.table(DMRtest_state,file=args[4],sep="\t",row.names=F,quote=F)

