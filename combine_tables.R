library(purrr)

source("modELISA_lib.R")
tableFile <- "file-to-info.csv"
if(! file.exists(tableFile)) stop('Error -- tableFile file not found.')
file_table <- read.csv(tableFile,stringsAsFactors=FALSE)




c_row <- file_table[2,]
tablepath=paste0("table_parts/",sprintf("%02d",c_row$Plate_in_file),"_",tools::file_path_sans_ext(c_row$File_name),".tab")
raw_table <- read.table(tablepath)

sample_cols <- c("row_A","row_B","row_C","row_D","row_E","row_F","row_G","row_H")
standar_row=min(which(c_row[,sample_cols]=="p24-standard"))

standar_table <- standarize_table(raw_table,standar_row)
standar_table$antibody <- as.list(c_row[,sample_cols])
standar_table$sample <- c_row$Sample_in_plate



standar_table

kk <- apply(file_table,1, function(x){
  plate <- as.numeric(x["Plate_in_file"])
  filename <- x["File_name"] 
  tablepath <- paste0("table_parts/",sprintf("%02d",plate),"_",tools::file_path_sans_ext(filename),".tab")
  raw_table <- read.table(tablepath)
  sample_cols <- c(x["row_A"],x["row_B"],x["row_C"],x["row_D"],x["row_E"],x["row_F"],x["row_G"],x["row_H"])
  standar_row <- min(which(sample_cols=="p24-standard"))
  standar_table <- standarize_table(raw_table,standar_row)
  standar_table$antibody <- as.list(sample_cols)
  standar_table$sample <- rep(x['Sample_in_plate'],8)
  return(standar_table)
})
c_data <-reduce(kk,rbind)
c_data
c_data %>% mutate(
  r1 <- mean(c(.V1,.V2,.V3,.V4))
#  r2 <- mean(V5,V6,V7,V8),
#  r3 <-mean(V9,V10,V11,V12)
)
sapply(t(c_data[,1:4]), mean)
