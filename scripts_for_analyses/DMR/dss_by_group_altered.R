library(bsseq)
library(DSS)
######################
#Rscript dss_human_celltype.R number_of_randomly_chosen_samples

## take input RDS file as first argument
args <- commandArgs(trailingOnly=TRUE)
rds_file <- args[1]
bs.combine <- readRDS(rds_file)

sample_names <- sampleNames(bs.combine)

min.cov <- 5
#ind.num <- 2
bs.cov <- getCoverage(bs.combine,type="Cov")
#bs.good.loci <- which(rowSums(bs.cov[,ind.REF] >= min.cov) >= round(length(ind.neun)*ind.per,0) & rowSums(bs.cov[,ind.olig] >= min.cov) >= round(length(ind.olig)*ind.per,0))
bs.good.loci <- which(bs.cov[,"healthy"] >= min.cov  & bs.cov[,"altered"] >= min.cov)
bs.combine.filtered <- bs.combine[bs.good.loci,]
bs.combine.filtered <- sort(bs.combine.filtered)
bs.combined.filtered.REF.CKD <- bs.combine.filtered[,c("healthy","altered")]


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
bs.cov.filtered <- getCoverage(bs.combined.filtered.REF.CKD)
write.table(bs.cov.filtered,file=args[4],sep='\t',row.names=F,quote=F)

dmlTest.sm = DMLtest(bs.combined.filtered.REF.CKD, group1="healthy", group2="altered",smoothing=TRUE, ncores=8)
dmls = callDML(dmlTest.sm, delta=0.001,p.threshold=1)
dmrs = callDMR(dmlTest.sm, delta=0.1, p.threshold=0.05,minCG=3)
dmls$pos <- format(dmls$pos, scientific = FALSE)
dmrs$start <- format(dmrs$start, scientific = FALSE)
dmrs$end <- format(dmrs$end, scientific = FALSE)
write.table(dmls,file=args[2],sep='\t',row.names=F,quote=F)
write.table(dmrs,file=args[3],sep='\t',row.names=F,quote=F)
