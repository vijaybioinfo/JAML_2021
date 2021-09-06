#!/usr/bin/R

###############
# Global data #
###############

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2021-08-18
# ---

# This script creates plots for our CD45+ JAML+ mice cells


## Locking object
source("/home/ciro/scripts/clustering/R/utilities.R")
configs = list(
  list(name = "SiEs12_Mo_1_CD45_3H_xdoublets", sufix = "_seurat_mean0.01_pct20", pc = 15)
)
for (config in configs) {
  mycells_f = paste0(
    "/mnt/beegfs/", Sys.info()[["user"]], "/", basename(getwd()), "_",
    gsub(".*seurat", "object_lock", config$sufix), "_pc", config$pc)
  if(file.exists(mycells_f)) next
  cat("Locking clustering\n")
  setwd(paste0("/home/ciro/ad_hoc/amica/results/clustering/", config$name))
  mycells = complex_object_fetch(
    object_or_dir = paste0(".object_stem", config$sufix, ".rds"),
    id = paste0(config$sufix, "_pc", config$pc, ".rds"), verbose = TRUE
  )
  tvar <- grep("dim_", colnames(mycells@meta.data), inver = TRUE)
  mycells@meta.data <- mycells@meta.data[, tvar]
  names(mycells@reductions) <- gsub("_.*", "", names(mycells@reductions))

  cat("Saving to:", mycells_f, "\n")
  saveRDS(mycells, file = paste0(mycells_f, ".rds"))
}

if(requireNamespace("crayon", quietly = TRUE)){
  cyan = crayon::cyan; red_bold = crayon::red$bold
}else{ cyan = red_bold = c }

{ cat(red_bold("------------------ Setting global parameters -----------------\n"))

  ### Objects names ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cat(cyan("------------------ File names and initial objects\n"))
  fig_dir = "/home/ciro/ad_hoc/amica/results/SiEs12_figures"
  redu = list(umap = c("UMAP_1", "UMAP_2"))
  sc_cd45jaml_clust = "RNA_snn_res.0.1"
  global_objects_f = c(
    gcolours = "/home/ciro/amica/info/global_colours.csv",
    sc_cd45jaml = "/mnt/beegfs/ciro/SiEs12_Mo_1_CD45_3H_xdoublets_object_lock_mean0.01_pct20_pc15.rds"
  )
  if(!exists("include")) include = names(global_objects_f)

  ### Environment/packages ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dir.create(fig_dir, showWarnings = FALSE); setwd(fig_dir)
  cat("Working at (fig_dir):", fig_dir, "\n")

  cat(cyan("------------------ Packages and functions\n"))
  packages_funcs = c(
    "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
    "/home/ciro/scripts/handy_functions/devel/filters.R",
    "/home/ciro/scripts/handy_functions/devel/utilities.R",
    "/home/ciro/scripts/handy_functions/devel/plots.R",
    "/home/ciro/scripts/figease/source.R", # fig_global_objects
    "ggplot2", "cowplot", "patchwork", "Seurat"
  )
  loaded <- lapply(X = packages_funcs, FUN = function(x){
    cat("*", x, "\n")
    if(!file.exists(x)){
      suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
    }else{ source(x) }
  }); theme_set(theme_cowplot())

  ### Global objects ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cat(cyan("------------------ Objects\n"))
  object_names = fig_global_objects(
    global_objects_files = global_objects_f,
    object_names = ls(pattern = "_f$|^redu|_clust$|_ident$"),
    reading_fun = readfile, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1
  )

  source("/home/ciro/scripts/clustering/R/utilities.R")
  celltype_subset = c(
    "0" = "CD45[Cd3g]",
    "1" = "CD45[Siglech]",
    "2" = "CD45[Pdcd1]",
    "3" = "CD45[Mafb]", # myeloid
    "4" = "CD45[Ifitm1]",
    "5" = "CD45[Xcr1]"
  )
  sc_cd45jaml_ident = list(
    celltype = gsub("\\[.*", "", celltype_subset),
    celltype_subset = celltype_subset,
    order = names(celltype_subset)
  )
  sc_cd45jaml_ident$colours = v2cols(names(sc_cd45jaml_ident[[2]]), NULL)
  sc_cd45jaml = ident_set_names(
    object = sc_cd45jaml,
    ident = sc_cd45jaml_ident,
    cluster_column = sc_cd45jaml_clust
  )
  sc_cd45jaml_ident <- ident_colours(
    sc_cd45jaml_ident, mdata = sc_cd45jaml@meta.data)

  print(format(
    object_names,
    justify = "centre", quote = FALSE)
  ); rm(include)
}

{ cat(red_bold("### Secondary global variables ### %%%%%%%%%%%%%%%%%%%%%%%%\n"))
  # As you progress they appear in many tasks and so become more or less global
}

{ cat(red_bold("### Markers: violins/dim. red. ### %%%%%%%%%%%%%%%%%%%%%%%%\n"))
  dir.create("dim_reduction_markers"); dir.create("violins")
  mdata = joindf(sc_cd45jaml@meta.data,
   as.data.frame(sc_cd45jaml@reductions$umap@cell.embeddings))
  fconfigs = list(
    list(result_id = "dim_reduction_markers/sc_cd45jaml/",
      edata = "sc_cd45jaml@assays$RNA@data", metadata = "mdata",
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("Pdcd1", "Tcf7"),
      colors = c("lightgrey", "blue")
    )
  )
  pp_markers = fig_plot_base(
    fconfigs, return_plot = FALSE, verbose = TRUE,
    theme_extra = function(x) x + geom_point(size = 0.4) + labs(x = "Dim 1", y = "Dim 2")
  )
  # pp_markers = fig_plot_scatters(fconfigs, return_plot = FALSE, verbose = TRUE) # Same as
  for (i in 1:length(fconfigs)) {
    fconfigs[[i]]$result_id = "violins/sc_cd45jaml/"
  }
  fconfigs[[1]]$axis_x = "cluster"
  fconfigs[[1]]$sample_filter = c("cluster", "0", "2")
  fconfigs[[1]]$features = c("Pik3cd", "Nkg7")
  source("/home/ciro/scripts/handy_functions/devel/plots.R")
  pp_violins = fig_plot_violins(
    fconfigs, return_plot = TRUE,
    verbose = TRUE,
    theme_extra = labs(x = NULL, y = "Seurat Normalized"),
    colour_by = "pct",
    couls = couls_opt$red_gradient$white
  );
}
