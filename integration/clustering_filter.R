#!/usr/bin/R

cat("****** Alignment - filter\n") ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("/home/ciro/scripts/handy_functions/devel/filters.R") # filters_subset_df
source("/home/ciro/scripts/handy_functions/devel/utilities.R") # make_list

n_features = 2000; n_pcs = 30
output_dir = "/home/ciro/large/amica/results/integration"
cluster_dir = paste0(output_dir, "/hvg", n_features, "_pc", n_pcs)
setwd(cluster_dir)
cat("Working at:", getwd(), "\n")
system("ls -loh")

if(file.exists("object.rds") && !exists('tcells_combined')) tcells_combined = readRDS("object.rds")

scells_names = filters_subset_df(
  c("integrated_snn_res.0.4", "-6", "-9", "-10"),
  tcells_combined@meta.data,
  v = TRUE
)
table(tcells_combined@meta.data$orig.ident)
print(table(tcells_combined@meta.data[scells_names, "integrated_snn_res.0.4"]))

fname = head(list.files(
  path = "..",
  pattern = "cells_selected_[0-9]{1,}sets\\.",
  full.names = TRUE
), 1)
scells0 = readRDS(fname)
scells = make_list(
  tcells_combined@meta.data[scells_names, ],
  "orig.study", "cellname"
)
cat("Previous:"); str(scells0)
cat("New:"); str(scells)
for(i in 1:length(scells)){
  j <- which(grepl(gsub("_.*", "", names(scells)[i]), names(scells0)))
  cat(" *", names(scells)[i], "renamed to", names(scells0)[j], "\n")
}
names(scells) <- names(scells0)
names(scells)

saveRDS(scells, file = gsub("[0-9]{1,}sets", "filtered", fname))
