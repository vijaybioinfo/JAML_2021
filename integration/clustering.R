#!/usr/bin/R

################################################
# Seurat integration: alignment and clustering #
################################################

# This scripts performs the actual integration and calculates the clustering.

# screen -S sia
# qsub -I -q default -l nodes=1:ppn=1 -l mem=80gb -l  walltime=36:00:00 -N sia
cat("****** Alignment\n") ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_dir = "/home/ciro/large/amica/results/integration"
dir.create(output_dir); setwd(output_dir)
cat("Working at:", getwd(), "\n")
system("ls -loh")

source("/home/ciro/amica/scripts/integration_clustering_env.R")

cat("========================== Loading data\n") ### %%%%%%%%%%%%%%%%%%%%%%%%%%%
scells_f = tail(list.files(pattern = "cells_selected_"), 1)
sfeatures_f = list.files(path = "harmonise_features", pattern = "updated", full.names = TRUE)
edata_list_f = c(
  "/home/ciro/large/simon/raw/guo_GSE99254_all_tpm.rdata",
  "/home/ciro/large/simon/raw/zheng_GSE98638_all_tpm.rdata",
  "/home/ciro/large/simon/raw/zhang_GSE108989_all_tpm.rdata",
  "/home/ciro/large/simon/raw/sadefeldman_GSE120575_all_tpm.rdata",
  "/home/ciro/large/simon/raw/puram_GSE103322_all_tpm.rdata",
  "/home/ciro/large/simon/raw/jerbyarnon_GSE115978_tpm.rdata",
  ## initial datasets ##
  "/home/ciro/large/simon/raw/zhang_GSE146771_10x_cp10k.rdata",
  "/home/ciro/large/simon/raw/maynard_GDrive_Tcells_umi.rdata",
  "/home/ciro/large/simon/raw/oh_GSE149652_cd4_cp10k.rdata"
)
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
set_file_names = function(x){
  gsub(
    pattern = "_all.*|_annot.*|_Tc.*|_cd.*|_10x.*|_tpm.*",
    replacement = "", basename(x), ignore.case = TRUE
  )
}
names(edata_list_f) = set_file_names(basename(edata_list_f))
names(metadata_list_f) = set_file_names(basename(metadata_list_f))

cat("Using cells:", scells_f, "\n")
scells = readRDS(scells_f)
sfeatures = readRDS(sfeatures_f)
names(sfeatures) = set_file_names(names(sfeatures))
names(scells) = set_file_names(names(scells))

mdata_list = lapply(X = metadata_list_f, FUN = readfile)
edata_list = lapply(X = edata_list_f, FUN = readfile)

cat("========================== Creating objects\n") ### %%%%%%%%%%%%%%%%%%%%%%%
objects_list <- lapply(
  X = setNames(names(sfeatures), names(sfeatures)),
  FUN = function(x) {
    cat("-----------------------", x, "\n")
    expr_data = edata_list[[x]]
    cat("Renaming features\n")
    rownames(expr_data) = sfeatures[[x]]
    cat("Create object\n"); timestamp()
    sobject <- CreateSeuratObject(
      counts = expr_data[, scells[[x]]],
      project = paste0("seurat_aligned_", x),
      meta.data = mdata_list[[x]][scells[[x]], ]
    ); timestamp()
    cat("Normalise\n")
    sobject <- NormalizeData(sobject)
    cat("Highly variable genes\n")
    sobject <- FindVariableFeatures(sobject, selection.method = "vst")
    sobject
})

cat("========================== Integration\n") ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exists("edata_list")) rm(edata_list)
n_features = 2000; n_pcs = 30
sufix = gsub("cells_selected_|[0-9]{1,}sets|.rds", "", scells_f)
if(sufix != "") sufix <- paste0("_", sufix)
cluster_dir = paste0(output_dir, "/hvg", n_features, "_pc", n_pcs, sufix)
dir.create(cluster_dir); setwd(cluster_dir)
cat("Working at:", getwd(), "\n")

features <- SelectIntegrationFeatures(object.list = objects_list, nfeatures = n_features)
tcells_anchors <- FindIntegrationAnchors(object.list = objects_list, anchor.features = features)
# Create the integrated data assay
tcells_combined <- IntegrateData(anchorset = tcells_anchors)
str(tcells_combined@assays$integrated@data)

# Assay (corrected data) for downstream analysis. The original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(tcells_combined) <- "integrated"

cat("========================== Standard workflow for clustering\n") ### %%%%%%%
# no need to regress for quality: https://github.com/satijalab/seurat/issues/3579
tcells_combined <- ScaleData(tcells_combined)
str(tcells_combined@assays$integrated@scale.data)
tcells_combined <- RunPCA(tcells_combined, npcs = 40)

pdf("pc_sdev.pdf", width = 10)
ElbowPlot(object = tcells_combined, ndims = 40)
graphics.off()

resolutions = c(0.1, 0.2, 0.4, 0.5, 0.6, 0.8)
tcells_combined <- RunUMAP(tcells_combined, reduction = "pca", dims = 1:n_pcs)
tcells_combined <- FindNeighbors(tcells_combined, reduction = "pca", dims = 1:n_pcs)
tcells_combined <- FindClusters(tcells_combined, resolution = resolutions)

cat("----------------------- Writing object to disk\n")
saveRDS(tcells_combined, file = "object.rds")
DefaultAssay(tcells_combined) <- "RNA"

cat("========================== Visualisation\n") ### %%%%%%%%%%%%%%%%%%%%%%%%%%
str(tcells_combined@meta.data)
cat("----------------------- UMAP of known variables\n")
dir.create("umaps")
vars2viz = c(
  "orig.tissue", "orig.study", "orig.library_method", "orig.disease", "orig.celltype",
  grep(pattern = "integrated_snn_res", colnames(tcells_combined@meta.data), value = TRUE)
)
for(i in vars2viz){
  cat(i, "\n")
  p1 <- DimPlot(tcells_combined, reduction = "umap", group.by = i)
  pdf(paste0("umaps/", gsub("orig.", "", i), ".pdf"))
  print(p1)
  graphics.off()
}

cat("----------------------- Resolutions\n")
pdf("clustree.pdf", height = 10)
clustree::clustree(tcells_combined, prefix = "integrated_snn_res.0.")
graphics.off()

pp = lapply(X = paste0("integrated_snn_res.", resolutions), FUN = function(x) {
  DimPlot(tcells_combined, reduction = "umap", group.by = x, label = TRUE, repel = TRUE) +
    theme_void() + theme(legend.position = "none")
})
pdf("umaps/_resolutions.pdf", width = 10)
plot_grid(plotlist = pp, nrow = 2)
graphics.off()

cat("----------------------- Markers\n")
markers <- read.csv('/home/ciro/simon/info/markers.csv', stringsAsFactors = FALSE)
mymarkers <- c(
  'FOXP3', 'IL2RA', 'TNFRSF9', 'TNFRSF18', 'DUSP4', 'CCR8', 'IL1R2', 'IKZF2',
  'ENTPD1', 'LAG3', 'TIGIT', 'CTLA4', 'PDCD1','TOX', "CXCL13", "JAML",
  "ITGAE", "ITGA1", "CD69",
  "PDCD1", "ICOS", "HAVCR2", "TNFRSF18", "TNFRSF9", "TIGIT", "LAG3", "VSIR"
)
mymarkers = unique(c(mymarkers, markers[, 1]))
tvar = mymarkers %in% rownames(tcells_combined@assays$RNA)
cat("Not found:", mymarkers[!tvar], "\n")
mymarkers <- mymarkers[tvar]

dir.create("markers_explore")
mycols_list = list(
  grey_blue=c("lightgrey", "blue"),
  blue_maroon=c("#edf8fb", "#6e016b"),
  gray_pink=c("#f1eef6", "#91003f"),
  lemon_blue=c("#ffffcc", "#0c2c84"))
# lapply(unlist(mycols_list), plotrix::color.id)
for(pt_size in c(0.2, 0.4)[2]){ # original is 0.1
  for(j in names(mycols_list)[1]){
    for(i in 1:length(mymarkers)){
      cat(mymarkers[i], "\n")
      fname = paste0('markers_explore/umap_size', pt_size, '_', j, '_', mymarkers[i])
      p <- FeaturePlot(
        tcells_combined, features = mymarkers[i], min.cutoff = 0,
        cols = mycols_list[[j]], pt.size = pt_size)
      pdf(paste0(fname, '.pdf'))
      print(p)
      graphics.off()
      pdf(paste0(fname, '_blank.pdf'))
      print(plot_blank(p))
      graphics.off()
    }
  }
}
