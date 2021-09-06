#!/usr/bin/R

######################################
# Seurat integration: cell selection #
######################################

# This scripts selects the cells that will be integrated for the AMICA paper

### Steps
# Select the cells: remove CD8+, maynard and puram only T cells (drop them)
# Doing QC: filtering genes
# Integration

# # Previous integration 1
# library(Seurat)
# setwd("/home/ciro/ad_hoc/tests")
# load("/home/ciro/large/simon/results/integration_tcells/seurat_9sets/PC15R0.3/integrated15PCs.RData")
# pancreas.integrated
# table(pancreas.integrated@meta.data$orig.celltype, useNA = 'ifany')
# table(pancreas.integrated@meta.data$orig.majorCluster, useNA = 'ifany')
# DimPlot(object = pancreas.integrated, reduction = "umap")
# dev.off()
# colnames(pancreas.integrated[[]])
# table(pancreas.integrated[['orig.treatment']])
# table(pancreas.integrated[['orig.BL_TRT']])

output_dir = "/home/ciro/large/amica/results/integration"
dir.create(output_dir); setwd(output_dir)

### Loading packages ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/file_reading.R")
source("/home/ciro/scripts/handy_functions/devel/filters.R")
source("/home/ciro/scripts/handy_functions/devel/utilities.R")

### Loading data ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# guo, zheng, zhang, sade, puram, arnon,
# discarding: Lambrechts, Savas, Li
metadata_list_f = c(
  "/home/ciro/large/simon/raw/guo_GSE99254_all_annotation.rdata",
  "/home/ciro/large/simon/raw/zheng_GSE98638_all_annotation.rdata",
  "/home/ciro/large/simon/raw/zhang_GSE108989_all_annotation.rdata",
  "/home/ciro/large/simon/raw/sadefeldman_GSE120575_all_annotation.rdata",
  "/home/ciro/large/simon/raw/puram_GSE103322_all_annotation.rdata",
  "/home/ciro/large/simon/raw/jerbyarnon_GSE115978_annotation.rdata",
  ## initial datasets ##
  "/home/ciro/large/simon/raw/zhang_GSE146771_10x_annotation.rdata",
  # "/home/ciro/large/simon/raw/zhang_GSE146771_ss2_annotation.rdata", # already in zhang_GSE108989
  "/home/ciro/large/simon/raw/luoma_GSE144469_CD4_annotation.rdata",
  "/home/ciro/large/simon/raw/maynard_GDrive_Tcells_annotation.rdata",
  "/home/ciro/large/simon/raw/oh_GSE149652_cd4_annotation.rdata",
  "/home/ciro/large/simon/raw/oh_GSE149652_cd8_annotation.rdata"
)
metadata_list <- lapply(X = metadata_list_f, FUN = readfile)
set_names <- gsub(
  pattern = "_all.*|_annot.*|_Tc.*|_cd.*|_10x.*|_tpm.*",
  replacement = "",
  basename(metadata_list_f),
  ignore.case = TRUE
)
metadata_list <- lapply(1:length(set_names), function(x){
  y <- metadata_list[[x]]
  y[["orig.set"]] <- set_names[x]
  y
})
names(metadata_list) <- set_names
t(sapply(metadata_list, function(x) table(x$orig.set, useNA = 'always') ))
metadata_celltypes <- t(sapply(metadata_list, function(x){
  sapply(c("CD4", "CD8", "Tcell"), function(y) any(grepl(y, c(x$orig.celltype, x$orig.majorCluster))) )
}))
metadata_celltypes <- data.frame(metadata_celltypes)
metadata_missed <- metadata_celltypes[which(rowSums(metadata_celltypes) < 2), ]
metadata_missed <- metadata_missed[!metadata_missed[, 1], ]
if(nrow(metadata_missed) > 0){
  print(lapply(metadata_list[rownames(metadata_missed)], function(x) table(x$orig.cluster) ))
  print(lapply(metadata_list[rownames(metadata_missed)], function(x) table(x$orig.celltype) ))
}
reshape2::melt(lapply(metadata_list, function(x) table(x$orig.celltype) ))
reshape2::melt(lapply(metadata_list, function(x) table(x$orig.library_method) ))
reshape2::melt(lapply(metadata_list, function(x) table(x$orig.treatment) ))
reshape2::melt(lapply(metadata_list, function(x) table(x$orig.tissue) ))

### Defining number of studies ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# patty <- paste0(rownames(metadata_missed), collapse = "|")
# metadata_list <- metadata_list[!grepl(patty, names(metadata_list))]
# table(metadata_list[['puram_GSE103322']]$orig.tissue)
# tvar <- reshape2::melt(table(metadata_list[['puram_GSE103322']][, c("orig.celltype", "orig.treatment", "orig.tissue")]))
# tvar[tvar[, 1] %in% "Tcell", ]

### Filtering ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter_columns = list(
  c("orig.celltype", "CD4", "CD8", "Tcell"),
  c("orig.treatment", "Pre"),
  c("orig.tissue", "T", "Tumor", "cancer", "Lung", "noncancer") # Lung is for maynard, noncancer for Puram
)
metadata_list_subset = lapply(
  X = names(metadata_list),
  FUN = function(mname){
    cat(crayon::red("=========================="), mname, "\n")
    filters_subset_df(
      x = filter_columns, df = metadata_list[[mname]], verbose = TRUE
    )
  }
)
names(metadata_list_subset) <- names(metadata_list)
tvar <- sapply(X = metadata_list_subset, FUN = length)
reshape2::melt(tvar)
metadata_list_subset <- metadata_list_subset[tvar > 100]
reshape2::melt(lapply(metadata_list[names(metadata_list_subset)], function(x) table(x$orig.library_method) ))

### Saving selected data ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sufix <- paste0(9, "studies") # length(metadata_list_subset)
saveRDS(object = metadata_list_subset, file = paste0("cells_selected_", sufix, ".rds"))
list.files()

### Checking database ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sc_database = readr::read_delim('http://www.nxn.se/single-cell-studies/data.tsv', delim = '\t')
sc_database_loc = readr::read_delim("/home/ciro/amica/info/data_table.csv", delim = ",")

colnames(sc_database)
colnames(sc_database_loc)
xcolumns = c("Title", "Authors", "Reported cells total", "Tissue")
find_by = c("Title", "Title")
xcolumns <- unique(c(xcolumns, find_by))
fetched_studies = lapply(
  X = sc_database_loc[, find_by[1], drop = TRUE],
  FUN = function(x_attr){
    found <- grepl(pattern = x_attr, x = sc_database[, find_by[2], drop = TRUE])
    y <- sc_database[found, xcolumns]
  }
)
fetched_studies <- fetched_studies[sapply(fetched_studies, nrow) > 0]
fetched_studies <- t(sapply(fetched_studies, function(x) x ))
fetched_studies
