library(tidyverse)
library(purrr)

source("modELISA_lib.R")
tableFile <- "file-to-info.csv"
if(! file.exists(tableFile)) stop('Error -- tableFile file not found.')
file_table <- read.csv(tableFile,stringsAsFactors=FALSE)




#c_row <- file_table[2,]
#tablepath=paste0("table_parts/",sprintf("%02d",c_row$Plate_in_file),"_",tools::file_path_sans_ext(c_row$File_name),".tab")
#raw_table <- read.table(tablepath)

#sample_cols <- c("row_A_antigen","row_B_antigen","row_C_antigen",
#                 "row_D_antigen","row_E_antigen","row_F_antigen",
#                 "row_G_antigen","row_H_antigen")
#standar_row <- min(which(c_row[,sample_cols]=="p24-standard"))
#control_row <- min(which(c_row[,sample_cols]=="p24-standard"))

#table <-reblank_table(raw_table,standar_row)
#get_control_threshold(r_table,standar_row)



#standar_table <- standarize_table(raw_table,standar_row)
#standar_table$antibody <- as.list(c_row[,sample_cols])
#standar_table$sample <- c_row$Sample_in_plate



#standar_table

kk <- apply(file_table,1, function(x){
  plate <- as.numeric(x["Plate_in_file"])
  filename <- x["File_name"] 
  tablepath <- paste0("table_parts/",sprintf("%02d",plate),"_",tools::file_path_sans_ext(filename),".tab")
  raw_table <- read.table(tablepath)
  sample_cols <- c(x["row_A_sample"],x["row_B_sample"],x["row_C_sample"],
                    x["row_D_sample"],x["row_E_sample"],x["row_F_sample"],
                    x["row_G_sample"],x["row_H_sample"])  
  antigen_cols <- c(x["row_A_antigen"],x["row_B_antigen"],x["row_C_antigen"],
                   x["row_D_antigen"],x["row_E_antigen"],x["row_F_antigen"],
                   x["row_G_antigen"],x["row_H_antigen"])
  standar_row <- min(which(antigen_cols=="p24-standard"))
  standar_table <- reblank_table(raw_table,standar_row)
  control_threshold <- get_control_threshold(standar_table,standar_row)
  standar_table$antigen <- unname(as.list(antigen_cols))
  standar_table$sample <- unname(as.list(sample_cols))
  standar_table$control_threshold <- rep(control_threshold,8)
  return(standar_table)
})
c_data <-reduce(kk,rbind)
c_data


kk2 <- apply(c_data,1, function(x){
  logdilutions <- log10(rep(c(1/300,1/900,1/2700,1/8100),3))
  measurments <- as.numeric(c(x["V1"],x["V2"],x["V3"],x["V4"],
                              x["V5"],x["V6"],x["V7"],x["V8"],
                              x["V9"],x["V10"],x["V11"],x["V12"]))
  model <- lm(measurments ~ logdilutions)
  threshold <- as.numeric(x['control_threshold'])
  return((threshold-model$coefficients[1])/model$coefficients[2])
})
kk2

c_data$heat <- kk2
kk3 <- c_data[c_data$sample!='BB-POS',] %>% 
  transmute(antigen=simplify(antigen),sample=simplify(sample),heat=10^as.numeric(heat)) %>%
  mutate(sorting=as.numeric(sprintf("%04d",as.numeric(str_extract(sample,'\\d+'))))) %>%
  arrange(sorting,antigen) %>%
  mutate(sorting=NULL,sample=paste0("L",sorting),
         heat=replace(heat,(heat<0.0001) | (heat>0.005),0)) 
  

#kk4 <- spread(kk3, antigen, heat, fill = 0, convert = TRUE)

png(height = 3, width = 8,units = 'in', res=300, file = 'chanchan.png')
#forcats::fct_rev(forcats::fct_inorder(sample)
ggplot(kk3, aes(forcats::fct_inorder(sample),antigen, fill= heat)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Sample")
#  labs() +
  coord_equal()
dev.off()




