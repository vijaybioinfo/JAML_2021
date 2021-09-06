#!/usr/bin/R

cat("****** Alignment - extras\n") ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("/home/ciro/amica/scripts/integration_clustering_env.R")

n_features = 2000; n_pcs = 30
output_dir = "/home/ciro/large/amica/results/integration"
cluster_dir = paste0(output_dir, "/hvg", n_features, "_pc", n_pcs, "_filtered")
setwd(cluster_dir)
cat("Working at:", getwd(), "\n")
system("ls -loh")

if(file.exists("object.rds") && !exists('tcells_combined')) tcells_combined = readRDS("object.rds")
DefaultAssay(tcells_combined) <- "RNA"
resolution_i = "integrated_snn_res.0.4"
tcells_combined$cluster.study = ident_combine(tcells_combined[[]], c("orig.study", resolution_i))

if(!file.exists("clustree_compared.pdf")){
  tcells_init_mdata = readRDS("../hvg2000_pc30/object.rds")@meta.data
  tvar <- tcells_init_mdata[, grep("res.0.4|study", colnames(tcells_init_mdata))]
  colnames(tvar) <- paste0(colnames(tvar), "0")
  mdata_combn = joindf(tvar, tcells_combined@meta.data)
  mdata_combn <- mdata_combn[, grep("res.0.4|study", colnames(mdata_combn))]
  for(i in 1:ncol(mdata_combn)) mdata_combn[, i] <- as.character(mdata_combn[, i])
  mdata_combn[is.na(mdata_combn)] <- "Filtered"
  pdf("clustree_compared.pdf")
  print(clustree::clustree(mdata_combn, prefix = "integrated_snn_res.0."))
  graphics.off()
}

cat("Creating 4 celltype groups\n")
celltype4groups_cols = setNames(
  RColorBrewer::brewer.pal(n = 4, name = 'Dark2'),
  c("non-TREG", "TREG", "TRM", "non-TRM")
)
tvar <- as.character(tcells_combined@meta.data[, resolution_i])
tmp <- tvar %in% c("2", "3", "4", "7")
tcells_combined$celltype4groups = ifelse(tmp, "non-TRM", "non-TREG")
tcells_combined$celltype4groups[tvar %in% c("0", "8")] <- "TREG"
tcells_combined$celltype4groups[tvar %in% c("2", "4")] <- "TRM"
tcells_combined$celltype4groups <- factor(
  tcells_combined$celltype4groups, names(celltype4groups_cols)
)
print(table(tcells_combined$celltype4groups))

p1 <- DimPlot(
  object = tcells_combined, reduction = "umap",
  group.by = "celltype4groups",
  cols = celltype4groups_cols
) + labs(x = NULL, y = NULL)
str(p1$data)
pdf("umaps/_celltype4groups.pdf", width = 10, height = 10)
print(p1)
graphics.off()

p1 <- DimPlot(
  object = tcells_combined, reduction = "umap",
  group.by = resolution_i, split.by = resolution_i,
  pt.size = 0.01, ncol = 4
) + theme(legend.position = "none") + labs(x = NULL, y = NULL)
str(p1$data)
pdf("umaps/selected_resolution.pdf", width = 10, height = 10)
print(p1)
graphics.off()

### Markers ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

vars2viz <- c(resolution_i, "celltype4groups")
ddf = FetchData(tcells_combined, vars = c(vars2viz, mymarkers))
for(k in vars2viz){
  cat(k, "\n")
  for(i in 1:length(mymarkers)){
    cat(" -", mymarkers[i], "\n")
    fname <- paste0('markers/violin_', k, "_", mymarkers[i], '.pdf')
    # if(file.exists(fname)) next
    p <- violin(ddf, k, mymarkers[i], colour_by = "pct") +
      labs(x = NULL, y = "Seurat Normalized")
    pdf(fname)
    print(p)
    graphics.off()
    fname <- paste0('markers/violin_', k, "_", mymarkers[i], '_blank.pdf')
    pdf(fname)
    print(plot_blank(p))
    graphics.off()
  }
}

### Markers dotplots ### --------
dir.create("dotplots_explore")
# JAML (AMICA1), PDCD1 (PD-1), ICOS, HAVCR2 (TIM3), TNFRSF18 (GITR),
# TNFRSF9 (4-1BB), TIGIT, LAG3, VSIR (VISTA)
fconfigs = list(
  list(
    genes = c("ICOS", "TNFRSF18", "TNFRSF9", "TIGIT", "VSIR", "HAVCR2", "LAG3", "JAML", "PDCD1"),
    result_id = "dotplots_explore/"
  ),
  list(
    genes = c("ICOS", "TNFRSF4", "TNFRSF9", "TNFRSF18", "JAML"),
    result_id = "dotplots_explore/costim_"
  ),
  list(
    genes = c("CTLA4", "LAG3", "PDCD1", "HAVCR2", "TIGIT"),
    result_id = "dotplots_explore/coinhibit_"
  )
)

for(dot_scale in c(12, 15)){
  for(fconfig in fconfigs){
    for(k in c("celltype4groups")){
      cat(k, "\n")
      p2 <- Seurat::DotPlot(
        object = tcells_combined, features = fconfig$genes,
        cols = c('#fff4ba', '#ff0000'), group.by = k, dot.scale = dot_scale
      ) + geom_point(mapping = aes(size = pct.exp), shape = 1, color = "black") +
        coord_flip() + labs(x = NULL, y = NULL)
      p2 <- p2 + RotatedAxis()
      p2$guides$colour$title <- "Z-scored\nAverage\nExpression"
      p2$guides$size$title <- "Percent\nExpressed"
      pdf(paste0(fconfig$result_id, "seurat_size", dot_scale, "_", k, ".pdf"))
      print(p2)
      graphics.off()
      pdf(paste0(fconfig$result_id, "seurat_size", dot_scale, "_", k, "_blank.pdf"), width = 6)
      print(plot_blank(p2))
      graphics.off()
    }
  }
}

### GSEA ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/R/gsea_signature.R")
source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_title, custom_heatmap

result_id = "signatures_"
signatures_files = c(
  "/home/ciro/scripts/handy_functions/data/signatures_vijay_db.rdata"
)

globalsign = list()
signatures_files <- signatures_files[file.exists(signatures_files)]
tvar <- lapply(signatures_files, readfile, stringsAsFactors = FALSE)
tvar <- unlist(lapply(tvar, as.list), recursive = FALSE)
tvar <- sapply(tvar, function(x){ y <- x[!is.na(x)]; y[which(y != "")] })
tvar <- tvar[sapply(tvar, length) > 0]; #str(tvar)
signatures_list <- c(globalsign, tvar[!names(tvar) %in% names(globalsign)])
signatures_list = signatures_list[grep("clarke", names(signatures_list), ignore.case = TRUE)]
signatures_list = signatures_list[!grepl("ter5|ter3|ter4", names(signatures_list), ignore.case = TRUE)]
str(signatures_list)
signatures_subset <- clean_feature_list(
  mat = tcells_combined@assays$RNA@data, features = signatures_list, filterby = "p~2", v = TRUE
)
dim(tcells_combined@assays$RNA@data)
dim(tcells_combined@assays$integrated@data)

ddf <- vlist2df_diff(
  x = signatures_subset,
  y = signatures_list,
  delim = "Not found or < 2 pct"
)
str(headtail(ddf))
tvar <- which(ddf == "Not found or < 2 pct", arr.ind = TRUE)
ddf[unique(tvar[,1]), unique(tvar[,2])]
write.csv(ddf, file = paste0(result_id, "table.csv"))

table(tcells_combined$cluster.study)
tcells_combined <- signature_scoring(
  object = tcells_combined,
  prefix = paste0(result_id, "modulescore_integrated/"),
  lsignatures = signatures_subset,
  confounders = c(resolution_i, "orig.study", "cluster.study"),
  violins_color = "mean",
  verbose = TRUE
)

dname = paste0(result_id, "violins")
dir.create(dname)
ddf = read.csv(paste0(result_id, "modulescore/signatures.csv"), row.names = 1)
ddf = joindf(tcells_combined@meta.data, ddf)
ddf$Study = stringr::str_to_title(make_title(ddf$orig.study))
colby = "mean"
colby = "pct"
p <- violin(ddf, "cluster.study", "TRM_CLARKE.Score", colour_by = colby, couls = 2) +
  facet_wrap(~Study, scale = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(
    labels = setNames(gsub(".*\\-", "", levels(ddf$cluster.study)), levels(ddf$cluster.study))
  )
summary(p$data$pct)
pdf(paste0("violins/cluster.study_", colby, ".pdf"), width = 10)
print(p)
graphics.off()

gsea_results <- gsea_matrix(
  mat = expm1(tcells_combined@assays$RNA@data),
  groups = resolution_i,
  metadata = tcells_combined@meta.data,
  gsea_list = signatures_list,
  path = paste0(result_id, "gsea_integrated/"),
  verbose = TRUE
)
for(i in c(0.05, 0.2)){
  for(j in c(0.5, 1, 1.5)){
    x <- gsea_plot_summary(
      tests_list = gsea_results,
      path = paste0(result_id, "gsea_integrated/", resolution_i, "_"),
      padjthr = i,
      nesthr = j
    )
  }
}

### Markers ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefix = paste0('markers_', resolution_i, "_pairwise")
prefix = paste0('markers_', resolution_i)
dir.create(prefix)
idents_each <- matrix(levels(tcells_combined@meta.data[, resolution_i]))
idents <- if(grepl("pairwise", prefix)){
  gtools::combinations(nrow(idents_each), r = 2, v = idents_each[, 1], set = TRUE, repeats.allowed = FALSE)
}else{ idents_each }
Idents(tcells_combined) <- resolution_i; table(Idents(tcells_combined))
for(i in 1:nrow(idents)){
  identy <- idents[i, 1]
  if(ncol(idents) > 1) identy2 <- idents[i, 2] else identy2 <- NULL
  cat('================= Group(s)', show_commas(idents[i, ]), '\n')
  fname <- paste0(c(
    paste0(prefix, '/dgea'), identy, "vs",
    ifelse(is.null(identy2), "REST", identy2),
    '.csv'), collapse = "_")
  if(!file.exists(fname)){
    cmarkers <- try(FindConservedMarkers(
      object = tcells_combined, ident.1 = identy, ident.2 = identy2,
      grouping.var = "orig.study", logfc.threshold = 0.1))
    if(class(cmarkers) == "try-error") next
    cmarkers$min_avg_logFC <- apply(
      X = cmarkers[, grep('avg_logFC', colnames(cmarkers))],
      MARGIN = 1,
      FUN = function(x){
        ifelse(all(min(x) * x > 0), min(x), 0)
    })
    write.csv(cmarkers, file = fname)
  }else cmarkers <- read.csv(
    fname, stringsAsFactors = F, check.names = F, row.names = 1)
  head(cmarkers); dim(cmarkers)
  degs <- getDEGenes(
    cmarkers, pv = 0.05, upreg = T,
    pvtype = 'minimump_p_val', lfc.type = 'min_avg_logFC')
  degs <- rownames(cmarkers)[rownames(cmarkers) %in% degs]
  # str(cmarkers[degs, ])
  cat('Done!\n')
}

# Merge results
mymarkersf <- list.files(prefix, pattern = '^dgea.*csv', full.names = T)
cmarkers <- lapply(mymarkersf, read.csv, stringsAsFactors = F, check.names = F, row.names = 1)
