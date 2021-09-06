#!/usr/bin/R

##########################
# SiEs05 TCR re-analysis #
##########################

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2021-06-08
# ---

# This script will create plots of the TCR data for our AMICA paper

output_dir = "/home/ciro/ad_hoc/amica/results/tcr/SiEs05/specific"
dir.create(output_dir); setwd(output_dir)
system("ls -loh")

cat("------------------ Packages\n")
loaded <- lapply(
  X = c("Seurat", "ggplot2", "cowplot")[-1],
  FUN = function(x){
    cat("*", x, "\n")
    suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
})
theme_set(theme_cowplot())
source("/home/ciro/scripts/handy_functions/devel/file_reading.R") # readfile
source("/home/ciro/scripts/handy_functions/devel/utilities.R") # joindf
source("/home/ciro/scripts/handy_functions/devel/plots.R") # couls_opt
source("/home/ciro/scripts/handy_functions/devel/filters.R") # features_add_tag
source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list

redu = list(umap = c("UMAP_1", "UMAP_2"))
trm_clusters = c("2", "3", "8", "9")
treg_clusters = c("1", "6")

clone_data_f = "/home/ciro/ad_hoc/amica/results/tcr/SiEs05/tcr_metadata_2021-06-08.rds"
clone_data = readfile(clone_data_f, row.names = 1, stringsAsFactors = FALSE)
clone_data$clone_size = clone_data$clon.size.tag
tvar <- 50 # quantile(clone_data$clone_size, probs = 0.97, na.rm = TRUE)
clone_data$clone_size[clone_data$clone_size>tvar] <- unname(tvar)
clone_data$cluster = as.character(clone_data$RNA_snn_res.0.6)

sc_sies05_f = "/home/ciro/ad_hoc/simon/results/clustering_seurat/SiEs05_notcr/clustering/zetInfo/clustCells17PCs_30Ks_0.06667JD.RData"
sc_sies05 = readfile(sc_sies05_f)
sc_sies05@meta.data$cluster = as.character(sc_sies05$RNA_snn_res.0.6)

### Re-run TCR analysis ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# REPORT=/home/ciro/ad_hoc/amica/results/tcr/SiEs05
# mkdir --parents ${REPORT}
#
# AGGR=${REPORT}/aggr.csv
# echo "library_id,clonotypes,annotations,sample.sffx" > ${AGGR}
# unset FNAMES; FNAMES=(`ls -rd /home/ciro/ad_hoc/hayley/raw/NV023/COUNTS/*TCR`)
# for I in `seq 0 $((${#FNAMES[@]}-1))`; do
#   SNAME=$(basename ${FNAMES[${I}]}); FNAME=${FNAMES[${I}]}
#   echo "${SNAME},${FNAME}/outs/clonotypes.csv,${FNAME}/outs/filtered_contig_annotations.csv,$((${I}+1))" >> ${AGGR}
# done
# cat ${AGGR}
# Rscript /home/vfajardo/scripts/TCR_data_analysis/10X_preliminary_analyses/aggr_vdj.2.2.R \
#   --ReportsPath ${REPORT} \
#   --GenInputPath ${REPORT} \
#   --AggrTable ${AGGR}\
#   --SampleCountThold 2 --FreqThold 2
# ls -loh ${REPORT}
# TCRCODE=/home/vfajardo/scripts/TCR_data_analysis/10X_preliminary_analyses/preliminary_TCR_data_analysis.1.5.4.R
# TCRCODE=/home/ciro/scripts/crtcr/analysis.R
# SOBJECT=/home/ciro/ad_hoc/simon/results/clustering_seurat/SiEs05_notcr/clustering/zetInfo/clustCells17PCs_30Ks_0.06667JD.RData
# Rscript ${TCRCODE} \
#   --ReportsPath=${REPORT} \
#   --TCRContigs=${REPORT}/filtered_contig_annotations_aggr.csv \
#   --TCRClonotypes=${REPORT}/clonotypes_aggr.csv \
#   --SeuratObj=${SOBJECT} \
#   --Tags "c('RNA_snn_res.0.6', 'orig.tissue')"

### Clone sharing ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "clone_sharing/"; dir.create(result_id)
clone_data$tmp <- clone_data$cluster
clone_data$tmp[clone_data$tmp %in% trm_clusters] <- paste0(clone_data$tmp[clone_data$tmp %in% trm_clusters], " TRM")
clone_data$tmp[clone_data$tmp %in% treg_clusters] <- paste0(clone_data$tmp[clone_data$tmp %in% treg_clusters], " TREG")
clones_list <- make_list(
  x = clone_data[!is.na(clone_data$clonotype.tag), ], colname = "tmp", col_objects = "clonotype.tag"
)
length(unique(unlist(clones_list)))
uddf <- as.data.frame(ComplexHeatmap::list_to_matrix(clones_list))
dim(uddf)
uddf <- as.data.frame(ComplexHeatmap::list_to_matrix(clones_list))
pdf(paste0(result_id, "upset.pdf"), onefile = FALSE, width = 10)
print(UpSetR::upset(data = uddf, sets = colnames(uddf)))
dev.off()

### Clonality UMAP plot ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "scatter/"; dir.create("scatter")

cols = gg_color_hue(0:9)
cols[c(2, 7)] = c("#8dd3c7", "#ffffb3")
# setNames(
#   c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd"),
#   0:9
# )
pairs2plot = list(
  list(x = "cluster", y = "clone_size"),
  list(x = "clone_size"),
  list(x = "CD8A"),
  list(x = "CD8B"),
  list(x = "CD4"),
  list(x = "FOXP3"),
  list(x = "JAML")
)
ddf = Seurat::FetchData(object = sc_sies05, vars = c(redu[[1]], unique(unlist(pairs2plot))))
ddf = joindf(ddf, clone_data)
set.seed(27); scells = sample(1:nrow(ddf), round(nrow(ddf) * 0.5))
for(k in seq(length(pairs2plot))[1]){
  i = pairs2plot[[k]]$x
  j = pairs2plot[[k]]$y
  aesy = aes_string(x = redu[[1]][1], y = redu[[1]][2], color = i, size = j)
  p1 = ggplot(data = ddf, mapping = aesy)
  if(!is.null(j)){
    p1 <- p1 + geom_point() +
      geom_point(shape = 1, color = "black", alpha = 0.1, stroke = 0.4) +
      scale_color_manual(values = cols) +
      scale_radius(breaks = scales::pretty_breaks(n = 5), range = c(0, 3)) +
      guides(colour = guide_legend(override.aes = list(size = 6)))
  }else{
    p1 <- p1 + geom_point(size = 0.3) + scale_color_gradient(low = "lightgrey", high = "blue")
      # scale_color_gradientn(colours = couls_opt$red_gradient$divisive, breaks = scales::pretty_breaks(n = 5))
  }
  p1 <- p1 + labs(color = NULL, size = NULL)

  fname <- paste0(result_id, paste0(c(names(redu)[1], i, j), collapse = "_"))
  cat(fname, "\n")
  pdf(paste0(fname, ".pdf"))
  print(p1)
  graphics.off()
  pdf(paste0(fname, "_blank.pdf"))
  print(plot_blank(p1))
  graphics.off()
}

### Bar plot comparing JAML+ versus JAML- TRM ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id = "barplot_"
result_id = "barplot_trm/"; dir.create("barplot_trm")
clones_list_ov = overlap_list(clones_list)
clones_list_ov <- clones_list_ov[sapply(clones_list_ov, length)>0]
clones_list_ov_trm = clones_list_ov[grep("TRM", names(clones_list_ov))]
tmp <- which(clone_data$cluster %in% trm_clusters); length(tmp)
tmp <- which(clone_data$cluster %in% trm_clusters | clone_data$clonotype.tag %in% unname(unlist(clones_list_ov_trm)))
length(tmp)
trm_ddf = clone_data[tmp, ]
table(trm_ddf$cluster)
tmp = features_add_tag(
  lgenes = "JAML", annot = trm_ddf,
  mat = expm1(sc_sies05@assays$RNA@data) * 100,
  thresh = 1, v = TRUE
)
# table(tmp[, 1], tmp_1[, 1])
trm_ddf$TRM_JAML = tmp[, 1]
cols = c("Non-expanded" = "gray", "Expanded" = "red")
trm_ddf$Expansion = factor(trm_ddf$TCR.tag, names(cols))
saveRDS(trm_ddf, file = paste0(gsub("\\/$", "_", result_id), "metadata.rds"))

pct_res <- plot_pct(trm_ddf, groups = c("Expansion", "TRM_JAML"), v = TRUE, return_table = TRUE)
p2 <- pct_res$plot + scale_fill_manual(values = cols) + labs(fill = NULL) + theme_cowplot()
pct_res$table

fname <- paste0(result_id, "TRM_JAML")
cat(fname, "\n")
pdf(paste0(fname, ".pdf"), width = 6)
print(p2)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 5)
print(plot_blank(p2))
graphics.off()

tmp = reshape2::melt(table(trm_ddf$Expansion, trm_ddf$TRM_JAML))
tmp <- tmp[nrow(tmp):1, ]
tmp <- tmp %>% group_by(Var2) %>%
  mutate(
    Percent = scales::percent(value / sum(value)),
    label_ypos = cumsum(value)
  )
tmp$PN = paste(tmp$Percent, "/", tmp$value)

p2 = ggplot(data = tmp, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat='identity') +
  geom_text(aes(y = label_ypos, label = PN), vjust = 1.6, color = "white", size = 3.5)+
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_manual(values = cols) + theme_cowplot()

fname <- paste0(result_id, "TRM_JAML_stacked")
cat(fname, "\n")
pdf(paste0(fname, ".pdf"), width = 6)
print(p2)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 5)
print(plot_blank(p2))
graphics.off()

# Clonotype table
sum(trm_ddf$clonotype.tag == "clonotype1", na.rm=T)
str(clone_data[which(clone_data$clonotype.tag == "clonotype1"), ]) # I don't know what it's 150 but I only find 141
table(clone_data[which(clone_data$clonotype.tag == "clonotype1"), "cluster"])
clone_summary_list = lapply(
  X = c("TRM_JAML", "cluster", "orig.tissue"),
  FUN = function(x){
    reshape2::dcast(trm_ddf, paste0("clonotype.tag ~", x), value.var = "clon.size.tag")
})
colnames(clone_summary_list[[2]])[-1] <- paste0("Cluster_", colnames(clone_summary_list[[2]]))[-1]
clone_summary <- clone_summary_list %>% Reduce(function(x, y) left_join(x, y, by = "clonotype.tag"), .)
head(clone_summary)
write.csv(clone_summary, file = "clonotypes_table.csv")

dir.create("violin")
fconfigs = list(list(result_id = "violin/", ddf2plot = trm_ddf[which(trm_ddf$clone_size>1), c("TRM_JAML", "clone_size")]))
tmp = reshape2::melt(clone_summary[, 2:3])
tmp = tmp[which(tmp$value>1), ]
colnames(tmp) <- c("TRM_JAML", "clone_size")
fconfigs[[2]] <- list()
fconfigs[[2]]$result_id = "violin/recalculated"
fconfigs[[2]]$ddf2plot = tmp
tvar <- quantile(tmp$clone_size, probs = 0.97, na.rm = TRUE)
tmp$clone_size[tmp$clone_size>tvar] <- unname(tvar)
fconfigs[[3]] <- list()
fconfigs[[3]]$result_id = paste0("violin/recalculated_capped", tvar)
fconfigs[[3]]$ddf2plot = tmp

for(fconfig in fconfigs){
  for(i in c("t.test", "wilcox.test")){
    p2 = violin(fconfig$ddf2plot, "TRM_JAML", y = "clone_size") +
      ggpubr::stat_compare_means(method = i) +
      labs(x = NULL, y = "Clone size") + theme(legend.position = "none")
    fname <- paste0(c(fconfig$result_id, "TRM_JAML", i), collapse = "_")
    cat(fname, "\n")
    pdf(paste0(fname, ".pdf"), width = 6)
    print(p2)
    graphics.off()
    pdf(paste0(fname, "_blank.pdf"), width = 5)
    print(plot_blank(p2))
    graphics.off()
  }
}

### Volcano plot JAML+ versus JAML- TRM cells ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_id0 = "NV033_murine_CD8-Tcell_"
fcthr_i = 1
dgea_res_f = "/home/ciro/ad_hoc/amica/results/dgea/NV033_murine/comprs/treatment_CD8-Tcell/antiPD1vsantiAMICA1/results_antiPD1vsantiAMICA1_deseq2.csv"
groups = NULL; genes = NULL; result_id = paste0(result_id0, basename(dirname(dgea_res_f)), "_")
dgea_res_f = "/home/ciro/ad_hoc/amica/results/dgea/NV033_murine/comprs/treatment_CD8-Tcell/isotypevsantiPD1/results_isotypevsantiPD1_deseq2.csv"
groups = NULL; genes = NULL; result_id = paste0(result_id0, basename(dirname(dgea_res_f)), "_")
# Rscript /mnt/BioHome/ciro/scripts/dgea/R/dgea.R \
#   --method=mastlog2cpm \
#   --metadata=/home/ciro/ad_hoc/amica/results/tcr/SiEs05/specific/barplot_metadata.rds \
#   --expression_data=/home/ciro/ad_hoc/simon/results/clustering_seurat/SiEs05_notcr/clustering/zetInfo/clustCells17PCs_30Ks_0.06667JD.RData \
#   --group1='JAML+' --group2='JAML-' \
#   --hname='TRM_JAML' \
#   --output_dir=/home/ciro/ad_hoc/amica/results/tcr/SiEs05/specific \
#   --ctrans=log2cpm
dgea_res_f = "/home/ciro/ad_hoc/amica/results/tcr/SiEs05/specific/JAML+vsJAML-/mastlog2cpm_results.csv"
groups = c("JAML-" = "JAML-", "JAML\\+" = "JAML+")
genes = c("GZMB", "CXCL13", "PRF1", "CD3D", "IFNG", "HLA-DMA", "HLA-DRB1", "LAG3", "GPR25")
result_id = ""; fcthr_i = 0.25

dgea_res = readfile(dgea_res_f, check.names = FALSE, stringsAsFactors = FALSE)
dgea_res <- dgea_res[!is.na(dgea_res$log2FoldChange), ]
colnames(dgea_res) <- gsub("/HOME.*", "", colnames(dgea_res))
dgea_res <- dgea_res[, !grepl("pct_diff", colnames(dgea_res))]
mygenes <- show_found(genes, dgea_res$gene, verbose = TRUE)
source("/home/ciro/scripts/dgea/R/report.R")
if(is.null(groups)) groups = strsplit(basename(dirname(dgea_res_f)), "vs")[[1]]
preport = dgea_report(
  results = dgea_res,
  group1 = groups[1],
  group2 = groups[2],
  fcthr = fcthr_i,
  cols_filt = NULL,
  dtype = "CPM",
  return_report = "volcano",
  showgenes = genes
)

pvolc = preport$volcano
if(result_id == "") pvolc <- pvolc + coord_cartesian(ylim = c(0, max(preport$volcano$data$FDR)))

fname = paste0(result_id, "volcano")
pdf(paste0(fname, ".pdf"), width = 10, height = 10)
print(pvolc)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 10, height = 10)
print(plot_blank(pvolc))
graphics.off()

### dotplots ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir.create("dotplots")
sc_sies05$celltype = sc_sies05@meta.data$cluster
sc_sies05$celltype[sc_sies05$celltype %in% trm_clusters] <- "TRM"
sc_sies05$celltype[sc_sies05$celltype %in% treg_clusters] <- "TREG"
table(sc_sies05$celltype)
scells = filters_subset_df(c("celltype", "TRM", "TREG"), sc_sies05@meta.data, v = TRUE)

result_id = "dotplots/costim_"
genes = c("ICOS", "TNFRSF4", "TNFRSF9", "TNFRSF18", "TIGIT", "JAML")
result_id = "dotplots/coinhibit_"
genes = c("CTLA4", "LAG3", "PDCD1", "HAVCR2")

mygenes = show_found(genes, rownames(sc_sies05), v = TRUE)

k = "celltype"
p2 <- Seurat::DotPlot(
  object = sc_sies05[, scells], features = mygenes,
  cols = c('#fff4ba', '#ff0000'), group.by = k
) + geom_point(mapping = aes(size = pct.exp), shape = 1, color = "black") +
  coord_flip() + labs(x = NULL, y = NULL) + Seurat::RotatedAxis()
p2$guides$colour$title <- "Z-scored\nAverage\nExpression"
p2$guides$size$title <- "Percent\nExpressed"
pdf(paste0(result_id, k, ".pdf"))
print(p2)
graphics.off()
pdf(paste0(result_id, k, "_blank.pdf"), width = 6)
print(plot_blank(p2))
graphics.off()

### SiEs12 ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sc_sies12_clust = "RNA_snn_res.0.1"
sc_sies12_f = "/home/ciro/ad_hoc/amica/results/clustering/SiEs12_Mo_1_CD45_3H_xdoublets/.object_stem_seurat_mean0.01_pct20.rds"
sc_sies12 = readfile(sc_sies12_f)
sc_sies12_mdata_f = "/home/ciro/ad_hoc/amica/results/clustering/SiEs12_Mo_1_CD45_3H_xdoublets/.object_meta.data_seurat_mean0.01_pct20_pc15.rds"
sc_sies12_mdata = readfile(sc_sies12_mdata_f)
sc_sies12_dr_f = "/home/ciro/ad_hoc/amica/results/clustering/SiEs12_Mo_1_CD45_3H_xdoublets/.object_reductions_seurat_mean0.01_pct20_pc15.rds"
sc_sies12_dr = readfile(sc_sies12_dr_f)
names(sc_sies12_dr) <- gsub("_.*", "", names(sc_sies12_dr))
sc_sies12@meta.data = sc_sies12_mdata[rownames(sc_sies12@meta.data), ]
sc_sies12@reductions = sc_sies12_dr

## Dim red ## ------------------------------------------------------------------
result_id = "SiEs12_umap"
p1 <- Seurat::DimPlot(
  object = sc_sies12, reduction = "umap",
  group.by = sc_sies12_clust, pt.size = 0.8
) + labs(x = NULL, y = NULL)
pdf(paste0(result_id, ".pdf"))
print(p1)
graphics.off()
pdf(paste0(result_id, "_blank.pdf"))
print(plot_blank(p1))
graphics.off()

## Heatmap ## ------------------------------------------------------------------
ntop = Inf
result_id = "SiEs12_heatmap_markers"
marknames = paste0(
  gsub(".rds", "", gsub(".object_.*seurat_me", "seurat_me", sc_sies12_mdata_f)),
  gsub(".*res.", "_res", sc_sies12_clust), "/dgea_MAST_fc0.25_padj0.05_summary_stats.csv"
)
ident_info_i = list()
myobject = "sc_sies12"
master_column = sc_sies12_clust

mygenes <- readfile(marknames, stringsAsFactors = FALSE, row.names = 1)

if(!grepl("uniq", result_id)){
  mygenesl <- make_list(mygenes, colname = "cluster", col_objects = "gene")
  str(mygenesl)
  uddf <- as.data.frame.matrix(table(mygenes[, c("gene", "cluster")]))
  if(length(ident_info_i) > 0){
    uddf <- uddf[, names(ident_info_i[[1]])]
    colnames(uddf) <- paste0("C", colnames(uddf), " ", ident_info_i[[1]][colnames(uddf)])
    str(uddf)
  }
  pdf(paste0(result_id, "_upset.pdf"), onefile = FALSE)
  print(UpSetR::upset(data = uddf, sets = colnames(uddf)))
  dev.off()
}

mygenes$gene <- gsub("'", "", mygenes$gene_name)
tvar <- mygenes$Dpct > .0; table(tvar)
mygenes <- mygenes[tvar, ]
tvar <- mygenes$avg_logFC > .25; table(tvar)
mygenes <- mygenes[tvar, ]
if(grepl("uniq", result_id)){
  tvar <- !grepl("&", mygenes$sCluster); table(tvar)
  mygenes <- mygenes[tvar, ]
}
tvar <- if(length(ident_info_i) > 0) names(ident_info_i[[1]]) else gtools::mixedsort(unique(mygenes$cluster))
mygenes$cluster <- factor(mygenes$cluster, tvar)
mygenes <- mygenes[order(mygenes$cluster), ]
source("/home/ciro/scripts/clustering/R/utilities.R")
topgenes <- get_top_n(x = mygenes, n = ntop)
fname <- paste0(result_id, ifelse(nrow(topgenes) != nrow(mygenes), paste0("_top", ntop), ""))
topgenes <- topgenes[!duplicated(topgenes$gene), ]
genes <- gsub("'", "", topgenes$gene_name)
genes <- show_found(genes, rownames(eval(parse(text = myobject))), v = TRUE)
genesl <- make_list(x = topgenes, colname = "cluster", col_objects = "gene_name")

annor = reshape2::melt(lapply(genesl, gsub, pattern = "'", replacement = ""))
annor <- data.frame(RowGroup = as.character(annor$L1), row.names = as.character(annor$value))
write.csv(annor, file = paste0(fname, ".csv"))

source("/home/ciro/scripts/handy_functions/devel/plots.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
# png(paste0(fname, ".png"), width = 1500, height = 1700, res = 300)
pdf(paste0(fname, ".pdf"), width = 7, height = 9)
custom_heatmap(
  object = eval(parse(text = myobject)),
  rnames = genes,
  orderby = master_column,
  use_mean = master_column,
  # sample_it = c(cname = master_column, maxln = "-1000"),
  scale_row = TRUE,
  categorical_col = c(master_column),
  feature_order = TRUE,
  couls = ident_info_i$colours,
  hcouls = c('yellow', 'black', 'blue'),
  regress = c('nCount_RNA', 'percent.mt'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annor, annotation_names_row = FALSE, annotation_names_col = FALSE
)
graphics.off()

## Violins ## ------------------------------------------------------------------
dir.create("SiEs12_violins")
genes = c(
  "Tnfrsf9", "Lck", "Gzmb", "Gzmk", "Prf1", "Ifng", "Mki67", "Top2a", "Pdcd1",
  "Lag3", "Havcr2", "Tox"
)
mymarkers <- show_found(genes, rownames(sc_sies12), v = TRUE)
ddf = Seurat::FetchData(
  object = sc_sies12,
  vars = c(sc_sies12_clust, mymarkers),
  cells = filters_subset_df(c(sc_sies12_clust, "0", "2"), sc_sies12[[]], v = TRUE)
)
for(i in 1:length(mymarkers)){
  cat(" -", mymarkers[i], "\n")
  fname <- paste0('SiEs12_violins/', sc_sies12_clust, "_", mymarkers[i])
  p <- violin(ddf, sc_sies12_clust, mymarkers[i], colour_by = "pct") +
    labs(x = NULL, y = "Seurat Normalized")
  pdf(paste0(fname, '.pdf'))
  print(p)
  graphics.off()
  pdf(paste0(fname, '_blank.pdf'))
  print(plot_blank(p))
  graphics.off()
}
