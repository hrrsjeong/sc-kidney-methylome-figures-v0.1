library(tidyverse)
library(methylCIPHER)
library(prcPhenoAge)
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
library(ggpubr)
#install.packages("ggpubr")
source("~/DNAmClock_sample_code/functions.R")
setwd("/Users/hjeong/Projects/sciMetv2/methylation_clock/")


#clock
#[1] "a"                  "b"                  "file_base"          "sample_id"          "cell_type"          "state"             
#[7] "Age"                "sex"                "Basename"           "Female"             "total_missing"      "horvath1_missing"  
#[13] "horvath2_missing"   "nonPRC_PhenoAge"    "PRC_PhenoAge"       "Alcohol_McCartney"  "BMI_McCartney"      "Bocklandt"         
#[19] "PhenoAge...18"      "PhenoAge...19"      "DNAmTL"             "DunedinPoAm38"      "EpiTOC"             "epiTOC2"           
#[25] "Garagnani"          "Hannum"             "Horvath1"           "Horvath2"           "HRSInChPhenoAge"    "hypoClock"         
#[31] "Knight"             "LeeControl"         "LeeRefinedRobust"   "LeeRobust"          "Lin"                "Mayne"             
#[37] "MiAge"              "PEDBE"              "PhenoAge...38"      "Smoking_McCartney"  "VidalBralo"         "Weidner"           
#[43] "Zhang"              "Zhang2019"          "UniversalAge2_450K" "PCHorvath1"         "PCHorvath2"         "PCHannum"          
#[49] "PCPhenoAge"         "PCDNAmTL"           "PCPACKYRS"          "PCADM"              "PCB2M"              "PCCystatinC"       
#[55] "PCGDF15"            "PCLeptin"           "PCPAI1"             "PCTIMP1"            "PCGrimAge"          "Y.pred"            
#[61] "ENCen40"            "ENCen100"           "DNAmGDF15"          "DNAmB2M"            "DNAmCystatinC"      "DNAmTIMP1"         
#[67] "DNAmADM"            "DNAmPAI1"           "DNAmLeptin"         "DNAmPACKYRS"        "DNAmlogCRP"         "DNAmlogA1C"        
#[73] "DNAmGrimAge2"       "AgeAccelGrim2"     

clock3 <- read.table(file="clock_coefficient.major_celltype.txt",sep='\t')

#clock4 <- clock2[clock2$state == "control",]
clock3$group <- "epithelia"
clock3[clock3$cell_type == "B",]$group <- "non-epithelia"
clock3[clock3$cell_type == "EC",]$group <- "non-epithelia"
clock3[clock3$cell_type == "FIB",]$group <- "non-epithelia"
clock3[clock3$cell_type == "Lymphoid",]$group <- "non-epithelia"
clock3[clock3$cell_type == "Myeloid",]$group <- "non-epithelia"
clock3 <- clock3[clock3$cell_type != "PEC",]

clock4 <- clock3[clock3$group == "epithelia",]
clock5 <- clock3[clock3$group == "non-epithelia",]
clock4$AgeAccel <- residuals(lm(PCHorvath1~true_age, data=clock4, na.action = na.exclude))
clock5$AgeAccel <- residuals(lm(PCHorvath1~true_age, data=clock5, na.action = na.exclude))
clock6 <- rbind(clock4,clock5)
ggboxplot3 <- ggplot(clock3, aes(x=state,y=AgeAccel))+#,fill=condition))+
  #geom_smooth(method='lm',linewidth=0.5,aes(color=state))+
  #geom_histogram(position="identity", alpha=0.5)+
  #geom_bar(stat="identity", aes(fill=state))+
  geom_boxplot(aes(fill=state))+
  stat_compare_means(aes(group = state), label = "p.format")+
  #geom_text(mapping=aes(x = 2, y = 60,label = true_age), vjust = -1)+
  geom_point(aes(fill=state),color="black",shape=21,position = position_jitterdodge(0.15, dodge.width = .8), alpha = .5) +  #xlim(c(0,100))+
  #ylim(c(0,100))+
  #coord_cartesian(ylim=c(-20,20))+
  #stat_cor(method = "pearson",aes(color=state))+
  #facet_wrap(~state+sample_id,nrow = 3)+#, scales = "free")+
  facet_wrap(~group)+#, scales = "free")+
  #geom_abline(intercept = 0,slope=1,linetype="dashed")+  #geom_boxplot() +
  #geom_bar(position="fill")+#, stat="identity")+
  #geom_text(stat = "prop", position = position_fill(.5),size=3.5)+
  #stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE)+
  xlab("")+#Cell type")+ 
  ylab("Age acceleration")+
  #scale_x_continuous(trans='log10')+
  #scale_y_continuous(trans='log10')+
  #ggtitle("cell-type deconvolution")+
  #facet_wrap(~state, nrow=1)+#scales = "free",nrow=2)+
  theme_classic()+
  #scale_fill_manual(values=c("black","firebrick2"))+
  #scale_color_manual(values=c("black","firebrick2"))+ 
  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")+
  scale_color_aaas()+
  scale_fill_aaas()+
  
  theme_bw() + theme(
    plot.title = element_text(size=16),
    #legend.position='in',
    legend.position='none', 
    #legend.position = c(0.80, 0.6),
    #legend.justification='left',
    #legend.direction='vertical',
    plot.background = element_blank(),
    axis.text.x = element_text(angle = 45, color = "black",size=12,hjust = 1),
    #axis.text.x=element_text(color = "black", size = 12),
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
ggboxplot3


