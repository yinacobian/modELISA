#stats
library(tidyverse)
library(ggpubr)
#require(scales)


setwd("C:/Users/Yina/Documents/Sync/Bushman Lab/IBD Lyndsey/IBD-serum-antigens/modELISAmerge/modELISA/")
DT <- read.csv("all_dilutions_at_threshold_dilution.csv")
#DT <- DT %>% filter(antigen=="HIV-p24")

#Plot all Baseline-Flare group

#pdf("BF-linear-all.pdf")
DT <- read.csv("all_dilutions_at_threshold_dilution.csv")
DT <- DT %>% filter(GROUP=="BF")
#DT <- DT %>% filter(antigen=="Norovirus")
ggpaired(DT, x = "VISIT_TYPE", y = "heat",
         color="VISIT_TYPE",line.color = "gray", line.size = 0.8, facet.by = "antigen", id="SUBJECT_ID",
         xlab="", ylab="Dilution to threshold", main="Baseline - Flare group",
         palette = "viridis", 1:10) +
  stat_compare_means(paired = TRUE, label.y.npc = 0.95) 
#stat_compare_means(paired = TRUE, label.y.npc = 0.95, method="t.test") 
#  yscale("log10", .format = TRUE)
#dev.off()

#Plot all Baseline-Quiescent group

DT <- read.csv("all_dilutions_at_threshold_dilution.csv")
DT <- DT %>% filter(GROUP=="BQ")
ggpaired(DT, x = "VISIT_TYPE", y = "heat",
         color="VISIT_TYPE",line.color = "gray", line.size = 0.8, facet.by = "antigen", id="SUBJECT_ID",
         xlab="", ylab="Dilution to threshold", main="Baseline - Quiescent group",
         palette = "viridis", 1:10) +
  stat_compare_means(paired = TRUE, label.y.npc = 0.95) 
#stat_compare_means(paired = TRUE, label.y.npc = 0.95, method="t.test") 
#  yscale("log10", .format = TRUE)

#  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#  annotation_logticks(sides="l")


##Norovirus looks insteresting
DT <- read.csv("all_dilutions_at_threshold_dilution.csv")
DT <- DT %>% filter(antigen=="Norovirus")
ggpaired(DT, x = "VISIT_TYPE", y = "heat",
         color="VISIT_TYPE",line.color = "gray", line.size = 0.8, facet.by = "GROUP", id="SUBJECT_ID",
         xlab="", ylab="Dilution to threshold", main="Norovirus",
         palette = "viridis") +
#yscale("log10", .format = TRUE) +
  stat_compare_means(paired = TRUE) 

##CMV looks insteresting
DT <- read.csv("all_dilutions_at_threshold_dilution.csv")
DT <- DT %>% filter(antigen=="CMV")
ggpaired(DT, x = "VISIT_TYPE", y = "heat",
         color="VISIT_TYPE",line.color = "gray", line.size = 0.8, facet.by = "GROUP", id="SUBJECT_ID",
         xlab="", ylab="Dilution to threshold", main="CMV",
         palette = "viridis") +
  #yscale("log10", .format = TRUE) +
  stat_compare_means(paired = TRUE) 


##Coronavirus
DT <- read.csv("all_dilutions_at_threshold_dilution.csv")
DT <- DT %>% filter(antigen=="Coronavirus HKU1")
ggpaired(DT, x = "VISIT_TYPE", y = "heat",
         color="VISIT_TYPE",line.color = "gray", line.size = 0.8, facet.by = "GROUP", id="SUBJECT_ID",
         xlab="", ylab="Dilution to threshold", main="Coronavirus HKU1",
         palette = "viridis") +
#  yscale("log10", .format = TRUE) +

  stat_compare_means(paired = TRUE, method="t.test") 
#  stat_compare_means(paired = TRUE) 

##HIV-p24
DT <- read.csv("all_dilutions_at_threshold_dilution.csv")
DT <- DT %>% filter(antigen=="HIV-p24")
ggpaired(DT, x = "VISIT_TYPE", y = "heat",
         color="VISIT_TYPE",line.color = "gray", line.size = 0.8, facet.by = "GROUP", id="SUBJECT_ID",
         xlab="", ylab="Dilution to threshold", main="HIV-p24",
         palette = "viridis") +
  #  yscale("log10", .format = TRUE) +
  
  stat_compare_means(paired = TRUE, method="t.test") 
#  stat_compare_means(paired = TRUE) 

##HIV-p24
DT <- read.csv("all_dilutions_at_threshold_dilution.csv")
DT <- DT %>% filter(antigen=="GFP")
ggpaired(DT, x = "VISIT_TYPE", y = "heat",
         color="VISIT_TYPE",line.color = "gray", line.size = 0.8, facet.by = "GROUP", id="SUBJECT_ID",
         xlab="", ylab="Dilution to threshold", main="GFP",
         palette = "viridis") +
  #  yscale("log10", .format = TRUE) +
  
  stat_compare_means(paired = TRUE, method="t.test") 
#  stat_compare_means(paired = TRUE) 





