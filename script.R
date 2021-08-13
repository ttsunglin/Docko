#----library----
library(readr)
library(stringr)
library(dplyr)

#----variable----
import_exp <- "Your_project_name" #Name your project
file_location <- "C:\\Enter\\Your\\Docking\\Result_Folder" #The location of out.pdbqt and .log files

ligand_sdf_location <- "C:\\Enter\\Your\\separated\\sdf_location" #seperated sdf files containg ZINC ID
sort_threshold <- 5  # top _% of strongest binding interested

#----extract data from log file----
EXTRACTENERGY <- function(file_no){
  log_read <- readLines(paste(file_location, "\\", file_no, ".pdbqt_log.log", sep = ""))
  log_no <- str_extract(grep("Output will be", log_read, value = TRUE),"[0-9]+")
  line_of_energy <- log_read[which(is.na(str_locate(log_read,"affinity ")) == FALSE)[1]+3]
  log_energy <- str_trim(str_sub(str_trim(line_of_energy), start = 2, end = 20))
  
  ligand_read <- readLines(paste(ligand_sdf_location, "\\", file_no, ".sdf", sep = ""))
  ligand_no <- ligand_read[which(is.na(str_locate(ligand_read,"zinc_id")) == FALSE)[1]+1]
  out <- c(log_no, log_energy, ligand_no)
  out
}

energy_list <- t(apply(data.frame(c(1:1614)), 1, EXTRACTENERGY)) %>% as.data.frame()
colnames(energy_list) <- c("log_no", "log_energy", "ligand_no")

dir.create(import_exp)
write.table(energy_list, file = paste(".\\", import_exp, "\\", import_exp, 
                                      "_energy_list.csv", sep = ""), 
            sep = ",", col.names = TRUE, row.names = FALSE)

#----copying top matches out.pdbqt to a folder----
energy_list$log_energy <- sapply(energy_list$log_energy, as.numeric)
energy_sort <- arrange(energy_list, log_energy)

top_num  <- ceiling(nrow(energy_sort)/100*sort_threshold)
choose_limit <- energy_sort$log_energy[top_num]
energy_sort_chosen <- energy_sort[which(energy_sort$log_energy <= choose_limit),]

write.table(energy_sort_chosen, file = paste(".\\", import_exp, "\\", import_exp, "_top",
                                             sort_threshold,
                                             "p_energy_list.csv", sep = ""), 
            sep = ",", col.names = TRUE, row.names = FALSE)
file_list <- sapply(energy_sort_chosen$log_no, function(x) paste(x, "_out.pdbqt", sep = ""))

for(x in c(1:nrow(energy_sort_chosen))){
  
  file.copy(paste(file_location, "\\", file_list[x], sep = ""), 
            paste(".\\", import_exp, sep = ""))
}
