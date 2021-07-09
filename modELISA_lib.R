library(dplyr)

standarize_table <- function(table,control_row) {
  
  blank_mean <- mean(as.numeric(table[control_row,1:3]))
  
  adj_table <- (table-blank_mean) # %>% mutate_all(~replace(.,.x<0, 0))
  #adj_table <- tmp_table %>% mutate_all(~replace(.,.x<0, 0))
  
  control_reads <- as.numeric(adj_table[control_row,5:12])
  control_consentration <- c(240,80,27,9,240,80,27,9)
  linearMod <- lm(control_consentration ~ control_reads )
  adj_table %>% mutate_all(~predict(linearMod,data.frame(control_reads = .)))
}


reblank_table <- function(table,control_row) {
  blank_mean <- mean(as.numeric(table[control_row,1:3]))
  adj_table <- (table-blank_mean)
  return(adj_table)
}

get_control_threshold <- function(table,control_row,standar_file="neg-pos-concentrations.txt") {
  if (standar_file=='neg-pos-concentrations-2.txt') {
    control_reads <- as.numeric(table[control_row,5:12])
    control_concentration <- c(80,27,9,3,80,27,9,3)
  } else if (standar_file=='neg-pos-concentrations-3.txt') {
    control_reads <- as.numeric(table[control_row,5:12])
    control_concentration <- c(27,9,3,1,27,9,3,1)
  } else {
    control_reads <- as.numeric(table[control_row,5:12])
    control_concentration <- c(240,80,27,9,240,80,27,9)
  }
  linearMod <- lm(control_reads ~ control_concentration )
  return(predict(linearMod,data.frame(control_concentration = 10)))
}
