#' @title  CytoSpace2SPATA2
#' @author Dieter Henrik Heiland
#' @description Import CytoSpace outputs to SPATA2
#' @inherit
#' @param object SPATA2 object
#' @param cytospace.out CytoSpace output folder
#' @param reference Seurat Objects of the sibgle cell reference
#' @return A list with a single-cell object (SPATA2) and the regular SPATA2 object
#' @examples
#' @export
#'

CytoSpace2SPATA2 <- function(object,
                             cytospace.out,
                             reference){
  
  st <- object
  
  message(paste0(Sys.time(), "--- Read in Cell Type fractions -----"))
  setwd(cytospace.out)
  sample <- st %>% SPATA2::getSampleName()
  
  #Add Cells to org SPATA2 object
  cell_types_spot <- read.csv(paste0(cytospace.out,"/cell_type_assignments_by_spot.csv"))
  names(cell_types_spot)[1] <- "barcodes"
  names(cell_types_spot)[2:ncol(cell_types_spot)] <- paste0(names(cell_types_spot)[2:ncol(cell_types_spot)], "_CS")
  
  cell_types_spot <- left_join(data.frame(barcodes=SPATA2::getBarcodes(st)[[1]]), cell_types_spot, by="barcodes")
  cell_types_spot[is.na(cell_types_spot)]=0
  st <- SPATA2::addFeatures(st, cell_types_spot,overwrite = T)
  
  message(paste0(Sys.time(), "--- Create single cell SPATA2 object -----"))
  sc_loc <- read.csv(paste0(cytospace.out,"/assigned_locations.csv"))
  
  # jitter position
  sc_loc <- sc_loc %>% left_join(., sc_loc %>% group_by(SpotID) %>% count())
  sc_loc <- map_dfr(.x=sc_loc$SpotID %>% unique(), .f=function(x){
    a <- sc_loc %>% filter(SpotID==x)
    if(a$n[1]==1){a$row.new=a$row; a$col.new=a$col}else{
      a$row.new <- runif(nrow(a), unique(a$row)-0.6, unique(a$row)+0.6) %>% sample()
      a$col.new <- runif(nrow(a), unique(a$col)-0.6, unique(a$col)+0.6) %>% sample()
    }
    return(a)
  }, .progress=T)
  
  
  #Get count mat
  mat <- Seurat::GetAssayData(reference, "counts")
  mat <- mat[,sc_loc$OriginalCID]
  mat <- mat %>% as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("OriginalCID")
  mat <- sc_loc %>% dplyr::select(UniqueCID,OriginalCID ) %>% left_join(.,mat)
  
  #UMAP
  umap <- reference@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column("OriginalCID")
  umap <- 
    sc_loc %>% 
    dplyr::select(UniqueCID,OriginalCID ) %>% 
    left_join(.,umap) %>% dplyr::select(-2) %>% 
    mutate(sample=sample) %>% 
    rename("umap1":=UMAP_1) %>% 
    rename("umap2":=UMAP_2) %>% 
    rename("barcodes":=UniqueCID) %>% 
    dplyr::select(barcodes, sample, umap1, umap2)
  
  meta_info <- left_join(sc_loc[,1:2],reference@meta.data %>% as.data.frame() %>% rownames_to_column("OriginalCID"))
  meta_info$barcodes <- meta_info$UniqueCID
  
  counts <- mat[,3:ncol(mat)]
  counts[is.na(counts)]=0
  counts[1:10, 1:10]
  
  counts <- counts %>% t() %>% Matrix::Matrix(sparse = T)
  colnames(counts) <- mat$UniqueCID
  
  coords <- data.frame(barcodes=sc_loc$UniqueCID, x=sc_loc$col.new, y=sc_loc$row.new)
  
  
  message(paste0(Sys.time(), "--- Reshape space of single cells -----"))
  coords_org <- getCoordsDf(st)
  coords$x <- scales::rescale(coords$x, c(min(coords_org$x), max(coords_org$x)))
  coords$y <- scales::rescale(coords$y, c(min(coords_org$y), max(coords_org$y)))
  coords_df <- coords
  
  
  
  message(paste0(Sys.time(), "--- Initiate object -----"))
  sc_obj <- SPATA2::initiateSpataObject_CountMtr(coords_df = coords_df, count_mtr = counts, sample_name = sample,
                                                 RunPCA=T, FindNeighbors=F, FindClusters=F,
                                                 RunTSNE=F, RunUMAP=F)
  
  sc_obj <- SPATA2::addFeatures(sc_obj, data.frame(barcodes=sc_loc$UniqueCID, Cell_lables=sc_loc$CellType, SpotID=sc_loc$SpotID))
  sc_obj@used_genesets <- st@used_genesets
  img <- st@images
  img[[sample]]@coordinates <- SPATA2::getCoordsDf(sc_obj)
  sc_obj@images <-img
  sc_obj <- SPATA2::flipCoordinates(sc_obj, axis="x")
  
  #Adapt the umap coords
  sc_obj@dim_red[[sample]]$umap <- umap
  names(meta_info) <- paste0(names(meta_info), "_meta")
  meta_info$barcodes <- meta_info$UniqueCID_meta
  
  sc_obj <- SPATA2::addFeatures(sc_obj, meta_info)
  
  #Add method info 
  sc_obj@information$pxl_scale_fct <- st@information$pxl_scale_fct
  sc_obj@information$method <- st@information$method
  
  return(list(st, sc_obj))
  
}