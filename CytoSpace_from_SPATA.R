# Run Cell type Deconvolution with Cytospace
#Linux
driver <- "/media/example..."

#Create Input files from SPATA and Seurat

#1. Single Cell ref Counts:
sc <- readRDS(paste0(driver, "/reference.RDS"))
library(Seurat)
library(tidyverse)
counts <- Seurat::GetAssayData(sc, "counts") %>% as.data.frame() %>% rownames_to_column("GENES")
write.table(counts, file=paste0(driver, "/SingleCellRef/Ref_mat.txt"), row.names=F, sep="\t",quote=F)


# 2. Cell Labels
lables <- data.frame(SpotID=sc@meta.data %>% rownames(), CellType=sc@meta.data$annotation_level_4)
write.table(lables, file=paste0(driver, "/SingleCellRef/lables.txt"), row.names=F, sep="\t",quote=F)


# Files will be placed in the Cytopace folder
CS_folder <- paste0(driver, "/CytoSpace")



#3. Spatial Transcriptomics
#Load table
Table_of_SPATA2_objects <- readRDS(paste0(driver, "/Tabellen/patients.RDS"))



run.list <- map(.x=1:nrow(Table_of_SPATA2_objects), .f=function(i){
  
  sample <- Table_of_SPATA2_objects$sample[i]
  st <- readRDS(paste0(driver, Table_of_SPATA2_objects$SPATA[i]))
  
  counts_st <- SPATA2::getCountMatrix(st) %>% as.data.frame() %>% rownames_to_column("V1")
  coords <- SPATA2::getCoordsDf(st) %>% dplyr::select(barcodes, row, col) %>% dplyr::rename("SpotID":=barcodes)
  
  #Create a folder for the data
  sample_dir <- paste0(CS_folder, "/", sample)
  if(dir.exists(sample_dir)==F){dir.create(sample_dir)}
  
  #output files
  write.table(counts_st, file=paste0(sample_dir, "/counts.txt"), row.names=F, sep="\t", quote=F)
  write.table(coords, file=paste0(sample_dir, "/coords.txt"), row.names=F, sep="\t", quote=F)
  
  
  
  #create the conda run script
  
  
  run <- paste0("conda run -n cytospace cytospace --scRNA-path '" , paste0(driver, "/SingleCellRef/Ref_mat.txt"), 
                "' --cell-type-path '", paste0(driver, "/SingleCellRef/lables.txt"),
                "' --st-path '", paste0(sample_dir, "/counts.txt"), 
                "' --coordinates-path '", paste0(sample_dir, "/coords.txt"), "' -o '", sample_dir, "'")
  
  return(run)
  
  
})
saveRDS(run.list, paste0(CS_folder, "/run.list.RDS"))


#create a bash run script

string.start <- " #!/usr/bin/env bash -l \n \n \n "
run.code <- map_chr(.x=run.list, .f=~ paste0(.x, " \n "))
string.out <- base::paste0(string.start,"\n", paste0(run.code, collapse = "\n"), "\n \n print(\"Done\") \n \n")
base::writeLines(string.out, paste0(CS_folder, "/cytoSpace.sh"))


























