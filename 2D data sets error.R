# List all RData files starting with 'Datos_2D' in the specified directory

data_error=NULL

for (i in 1:4) {

  print(paste0("i = ",i))
  
  rdata_files <- list.files(path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/2do paper/codigo/Paper 2 Classification/2D/2D-Partially-observed-data", pattern = "^Datos_2D_c_15_miss.*\\.RData$", full.names = TRUE)
  
  file_name=basename(rdata_files[i])
  
  file_name_no_ext <- sub("\\.RData$", "", file_name)
  
  load(rdata_files[i])
  
  data_error=rbind(data_error,data.frame(case=file_name_no_ext, error=sum(miss_error)/N, AUC=AUC))
  
  rm(list = setdiff(ls(), "data_error"))

  gc()
  
}

data_error$case[order(data_error$error)]
