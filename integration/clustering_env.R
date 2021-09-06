#!/usr/bin/R

cat("========================== Loading packages\n") ### %%%%%%%%%%%%%%%%%%%%%%%
packages_funcs = c(
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  "/home/ciro/scripts/handy_functions/devel/filters.R",
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  "/home/ciro/scripts/handy_functions/devel/plots.R",
  "/home/ciro/scripts/handy_functions/devel/plots_dotplot.R",
  "/home/ciro/scripts/clustering/R/utilities.R",
  "ggplot2", "cowplot", "ggraph", "Seurat"
)
loaded <- lapply(X = packages_funcs, FUN = function(x){
  cat("*", x, "\n")
  if(!file.exists(x)){
    suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
  }else{ source(x) }
}); theme_set(theme_cowplot())
