## this code times how long the Pan Cancer datasets load on Server 9

## load libraries
library("Seurat")

## change the path to wherever your RDS files are
file_path <- "/media/seq-srv-05/vrc/Project/Project_Theo/pan_cancer_shiny/"

## function which records how long datasets load
time_load <- function(object){
	start_time <- Sys.time()
	readRDS(file = paste0(file_path, object))
	end_time <- Sys.time()
	data_loading <- print(end_time - start_time)
}

## load RDS objects and print loading time
for (i in 1:length(list.files(file_path, pattern = "*.rds"))) {
print(list.files(file_path, pattern = "*.rds")[i])
time_load(list.files(file_path, pattern = "*.rds")[i])
}
