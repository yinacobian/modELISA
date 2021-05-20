library(dplyr)


#tableFile <- commandArgs(trailingOnly = TRUE)
tableFile <- "table_parts/03_2020-12-10_AC_IBD_4-plates_SET1-L149-L150-L224-L225.tab"
if(! file.exists(tableFile)) stop('Error -- tableFile file not found.')
raw_table <- read.table(tableFile)
blank_mean <- mean(as.numeric(raw_table[1,1:3]))

adj_table <- (raw_table-blank_mean) # %>% mutate_all(~replace(.,.x<0, 0))
#adj_table <- tmp_table %>% mutate_all(~replace(.,.x<0, 0))

control_reads <- as.numeric(adj_table[1,5:12])
control_consentration <- c(240,80,27,9,240,80,27,9)

linearMod <- lm(control_consentration ~ control_reads )
adj_table %>% mutate_all(~predict(linearMod,data.frame(control_reads = .)))

standarize_table