#Load libraries

suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(DMwR))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(edgeR))
library(ggrepel)
library("ggsci")
library("ggplot2")
library("gridExtra")
library("ggpubr")
library(ggjoy)
library("scatterplot3d")
library(stringr)
#library(tidyverse)
library(ggstats)
library(ggplot2)
library(chameleon)

df_all2 <- read.table("combinedAll.stat.txt",sep='\t',header=T,na.strings = "NaN")
df_all2 <- df_all2[df_all2$Sample %in% c("SDKZ0030","SDKZ0034","SDKZ0036","SDKZ0037","SDKZ0043","SDKZ0057"), ]
head(df_all2)
#
df_all2_PT <- read.table("combinedAll_PT_DMR.stat.txt",sep='\t',header=T,na.strings = "NaN")
df_all2_PT <- df_all2_PT[df_all2_PT$Sample %in% c("SDKZ0030","SDKZ0034","SDKZ0036","SDKZ0037","SDKZ0043","SDKZ0057") ,]
#
#df_all2_PT <- df_all2_PT[df_all2_PT$nReads > 10,]
head(df_all2_PT)
#rownames(df_all2) <- df_all2$loci
#rownames(df_all2) <- df_all2$loci
head(df_all2)
df_all3 <- df_all2[!(df_all2$loci %in% df_all2_PT$loci),]
df_all2_PT$label <- "PT-DMR"
df_all3$label <- "non-PT-DMR"
dim(df_all2_PT)
dim(df_all3)
#dim(df_all2_PT)[1]
df_all4 <- rbind(df_all2_PT,df_all3[1:dim(df_all2_PT)[1],])
#df_all2$consensus_annotation3 <- factor(df_all2$consensus_annotation3,levels=c("PEC","PT","DCT","TAL","CNT","IC","PC","POD","EC","FIB","T","Myeloid","B"))
ggbox <- ggplot(df_all4, aes(x=PT_state,y=Entropy,fill =PT_state))+#,fill=condition))+
  #geom_point(size=0.005,aes(color=consensus_annotation3))+
  #geom_bar()+#outlier.shape = NA,aes(fill=cell_class))+
  #geom_bar(stat="identity",color="black", position=position_dodge())+#, stat="identity")+
  geom_boxplot()+
  stat_compare_means(
    aes(label = paste0("p = ", after_stat(p.format)))
  )+
  #stat_compare_means(
  #  aes(label = paste0("p = ", after_stat(p.format)))
  #)
  #geom_point(size=0.005,aes(color=L1))+
  #geom_text(size=0.5)+
  #xlim(c(0,100))+
  #ylim(c(0,100))+
  #coord_cartesian(ylim=c(0, 0.015))+
  #coord_flip(ylim=c(0, 0.015))+
  #coord_flip()+#ylim=c(0, 0.015))+
  
  
  #facet_wrap(~label)+#, scales = "free")+
#geom_abline(intercept = 0,slope=1,linetype="dashed")+  #geom_boxplot() +
#geom_bar(position="fill")+#, stat="identity")+
#geom_text(stat = "prop", position = position_fill(.5),size=3.5)+
#stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE)+
#xlab("Sample")+ 
#ylab("Cell type proportion")+
#scale_x_continuous(trans='log10')+
#scale_y_continuous(trans='log10')+
#ggtitle("cell-type deconvolution")+
facet_wrap(~label, nrow=1)+#scales = "free",nrow=2)+
#scale_fill_chameleon(minimal_saturation = 10,
#                      minimal_lightness = 30,
#                     maximal_lightness = 70, name="") +
#scale_color_continuous_qualitative()+
theme_classic()+
  #scale_fill_manual(values=c("black","firebrick2"))+
  #scale_color_manual(values=c("black","firebrick2"))+ 
  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")+
  #scale_color_npg()+
  #scale_x_discrete(limits = rev)+
  scale_fill_npg()+
  #guides(color=guide_legend(override.aes=list(size=3))) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  #theme(legend.text=element_text(size=12), legend.key.height=unit(14, 'pt'))+
  theme_bw() + theme(
    plot.title = element_text(size=16),
    #legend.position='in',
    legend.position='none', 
    #legend.position = c(0.80, 0.6),
    #legend.justification='left',
    #legend.direction='vertical',
    plot.background = element_blank(),
    axis.text.x = element_text(angle = 45, color = "black",size=14,hjust = 1),
    #axis.text.x=element_text(color = "black", size = 14),
    axis.text.y=element_text(color = "black", size = 14),
    strip.text = element_text(size = 14),
    text = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    #axis.title.y=element_blank(),
    #axis.title.y = element_text(size = 14),
    legend.title=element_text(size=12),
    legend.text=element_text(size=14),
    #legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggbox


