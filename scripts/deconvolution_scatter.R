
#####################

# Load libraries

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

setwd("/Users/hjeong/Projects/sciMetv2/deconvolution/")
df <- read.table("deconv_LOO_major_celltype.prediction.hypo_only.csv",sep='\t',header=T)
df$CellType2 <- df$CellType
df <- df[df$CellType2 != "IMMU",]
#df <- df[df$CellType2 != "ATL",]
df[df$CellType2 == "Lymphoid",]$CellType <- "T"
df2 <- df[df$CellType %in% c("PT-S1","PT-S2","PT-S3","PT-alt-1","PT-alt-2","PT-alt-3","PT-alt-4","TAL","TAL-alt-1","TAL-alt-2","TAL-alt-3"),]
head(df)
unique(df$CellType)
df$CellType <- factor(df$CellType,levels=c("CNT" , "DCT", "IC" ,"PC", "TAL", "FIB", "EC"  , "B", "T", "Myeloid" , "PEC",  "POD" , "PT" ))
df$Sample <- factor(df$Sample)
ggscatter <- ggplot(df, aes(x=sciMET_celltype_ratio,y=predicted_celltype_ratio,label=CellType))+#,fill=condition))+
  geom_point(size=1.5,aes(color=CellType,shape=Sample))+
  #geom_text_repel(max.overlaps=5,size=4)+
  #geom_text(size=0.5)+
  #xlim(c(0,0.45))+
  #ylim(c(0,0.45))+
  #facet_wrap(~Sample, scales = "free",nrow=2)+
  geom_smooth(method=lm,linewidth=0.5, se=TRUE, fullrange=TRUE)+
  stat_cor(method = "pearson",aes(label = ..r.label..),size=5)+#, label.x = 3, label.y = 30)+
  #geom_abline(intercept = 0,slope=1,linetype="dashed")+  #geom_boxplot() +
  #geom_bar(position="fill")+#, stat="identity")+
  #geom_text(stat = "prop", position = position_fill(.5),size=3.5)+
  #stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE)+
  xlab("Observed")+ 
  ylab("Predicted")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent)+
  scale_shape_manual(values=1:nlevels(df$Sample)) +
  #scale_x_continuous(trans='log10')+
#scale_y_continuous(trans='log10')+
#ggtitle("cell-type deconvolution")+
#scale_color_chameleon(minimal_saturation = 25,
#                      minimal_lightness = 20,
#                      maximal_lightness = 80) +
  scale_color_chameleon(minimal_saturation = 15,
                       minimal_lightness = 20,
                       maximal_lightness = 80, name="") +
  
  #scale_color_continuous_qualitative()+
  theme_classic()+
  #scale_fill_manual(values=c("black","firebrick2"))+
  #scale_color_manual(values=c("black","firebrick2"))+ 
  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")+
  #scale_color_npg()+
  #scale_fill_npg()+
  guides(color=guide_legend(override.aes=list(size=3))) +
  #theme(legend.text=element_text(size=12), legend.key.height=unit(14, 'pt'))+
  theme_bw() + theme(
    plot.title = element_text(size=12),
    #legend.position='in',
    #legend.position='none', 
    #legend.position = c(0.80, 0.6),
    #legend.justification='left',
    #legend.direction='vertical',
    plot.background = element_blank(),
    axis.text.x = element_text(color = "black",size=12),
    #axis.text.x=element_text(color = "black", size = 16),
    axis.text.y=element_text(color = "black", size = 12),
    strip.text = element_text(size = 12),
    text = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    #axis.title.x=element_blank(),
    axis.title.y = element_text(size = 12),
    legend.title=element_text(size=10),
    legend.text=element_text(size=10),
    #legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggscatter

##########################################

df <- read.table("DNAmAge_SampleLevel.csv",sep=',',header=T)
head(df)
colnames(df)
ggboxplot <- ggplot(df, aes(x=Age,y=PCGrimAge))+#,fill=condition))+
  #geom_histogram(position="identity", alpha=0.5)+
  geom_point(size=2,aes(color=state,shape=sex))+
  xlim(c(0,100))+
  ylim(c(0,100))+
  facet_wrap(~downsample_cov)+#, scales = "free")+
  geom_abline(intercept = 0,slope=1,linetype="dashed")+  #geom_boxplot() +
  #geom_bar(position="fill")+#, stat="identity")+
  #geom_text(stat = "prop", position = position_fill(.5),size=3.5)+
  #stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE)+
  #xlab("Sample")+ 
  #ylab("Cell type proportion")+
  #scale_x_continuous(trans='log10')+
  #scale_y_continuous(trans='log10')+
  #ggtitle("cell-type deconvolution")+
  #facet_wrap(~state, nrow=1)+#scales = "free",nrow=2)+
  theme_classic()+
  #scale_fill_manual(values=c("black","firebrick2"))+
  #scale_color_manual(values=c("black","firebrick2"))+ 
  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")+
  #scale_color_npg()+
  #scale_fill_npg()+
  
  theme_bw() + theme(
    plot.title = element_text(size=12),
    #legend.position='in',
    legend.position='right', 
    #legend.position = c(0.80, 0.6),
    #legend.justification='left',
    #legend.direction='vertical',
    plot.background = element_blank(),
    axis.text.x = element_text(color = "black",size=12),
    #axis.text.x=element_text(color = "black", size = 16),
    axis.text.y=element_text(color = "black", size = 12),
    strip.text = element_text(size = 12),
    text = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    #axis.title.x=element_blank(),
    axis.title.y = element_text(size = 12),
    legend.title=element_text(size=12),
    legend.text=element_text(size=12),
    #legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggboxplot

