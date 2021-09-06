#!/usr/bin/R

#############################################
# Seurat integration: feature harmonisation #
#############################################

# This scripts harmonises the feature names across datasets

# if(!require(geneSynonym)) devtools::install_github('oganm/geneSynonym')
output_dir = "/home/ciro/large/amica/results/integration"
dir.create(output_dir); setwd(output_dir)
dir.create("harmonise_features/")

### Loading packages ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("/home/ciro/scripts/handy_functions/devel/file_reading.R")
source("/home/ciro/scripts/handy_functions/devel/overlap.R")

### Loading data ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

features_names_f = "harmonise_features/features_names.rds"
if(file.exists(features_names_f)){
  edata_list <- lapply(X = edata_list_f, FUN = readfile)
  set_names <- gsub(
    pattern = "_all.*|_annot.*|_Tc.*|_cd.*|_10x.*|_tpm.*",
    replacement = "",
    basename(edata_list_f),
    ignore.case = TRUE
  )
  names(edata_list) <- set_names
  features_names = lapply(X = edata_list, FUN = function(x) unname(rownames(x)) )
  saveRDS(features_names, file = features_names_f)
}else{
  features_names <- readRDS(file = features_names_f)
}
str(features_names)
features_all = unique(unlist(features_names, use.names = FALSE))

# Finding synonyms
features_a_syn_f = "harmonise_features/features_a_syn.rds"
features_a_syn_f = "harmonise_features/features_a_ths.rds"
if(!file.exists(features_a_syn_f)){
  features_a_syn = if(grepl("ths", features_a_syn_f)){
    Seurat::GeneSymbolThesarus(symbols = features_all, timeout = 15, several.ok = TRUE)
  }else{ Seurat::GeneSymbolThesarus(symbols = features_all, timeout = 15) }
  saveRDS(features_a_syn, file = features_a_syn_f)
}else{
  features_a_syn <- readRDS(file = features_a_syn_f)
}
cat("synonyms", length(features_a_syn), "\n"); str(features_a_syn[sapply(features_a_syn, length) > 1])
cat("More than one:\n"); str(features_a_syn[sapply(features_a_syn, length) > 1])
cat("Not found:"); str(features_all[!features_all %in% names(features_a_syn)])
sapply(features_names, function(x) "GAGE12F" %in% x )
features_a_syn[["GAGE12F"]] # not found
# Seurat::GeneSymbolThesarus("GAGE-12F")
# Seurat::GeneSymbolThesarus("GAGE12F")
# Seurat::GeneSymbolThesarus("GAGE7B") # from GeneCards
# Seurat::GeneSymbolThesarus("GAGE12I")

# finding overlapping synonyms across datasets
if(!exists("features_u_all")) features_u_all = overlap_list(features_names)
features_u <- features_u_all[sapply(features_u_all, length) > 0]
str(features_u, vec.len = 2) # now take the ones not in every dataset
features_u <- features_u[stringr::str_count(names(features_u), "~") + 1 < length(features_names)]
features_symbols <- unique(unlist(features_u, use.names = FALSE))
cat("All\n"); str(features_all); cat("Not in all datasets\n"); str(features_symbols)

# Seurat::GeneSymbolThesarus(symbols = "CXXC11", timeout = 15, several.ok = TRUE)
features_u_syn_f = "harmonise_features/features_u_ths.rds";
features_u_syn_f = "harmonise_features/features_u_syn.rds";
if(!file.exists(features_u_syn_f)){
  features_u_syn = if(grepl("ths", features_u_syn_f)){
    Seurat::GeneSymbolThesarus(symbols = features_symbols, timeout = 15, several.ok = TRUE)
  }else{ Seurat::GeneSymbolThesarus(symbols = features_symbols, timeout = 15) }
  saveRDS(features_u_syn, file = features_u_syn_f)
}else{
  features_u_syn <- readRDS(file = features_u_syn_f)
}
cat("Not in all datasets\n"); length(x = features_symbols); cat("Non-updateable\n"); length(x = features_u_syn)
mean(features_u_syn %in% features_symbols)

check_gene = function(gpattern, features_names_list, features_thes){
  cat("Symbols", grep(gpattern, features_symbols, value = TRUE), "\n")
  cat("Updates", grep(gpattern, features_u_syn, value = TRUE), "\n")
  synl = features_thes[grep(gpattern, names(features_thes))]
  cat("Thesaurus "); str(synl)
  find_them = unique(c(unlist(synl), names(synl), unlist(strsplit(gpattern, "\\|"))))
  synl_insets = lapply(features_names_list, function(x) sort(x[x %in% find_them]) )
  cat("In eac dataset: "); str(synl_insets)
  # checking if there's recursivity
  gpattern_ext = paste0(find_them, collapse = "|")
  tvar <- grep(gpattern_ext, names(features_thes))
  if(length(tvar) > 1){ cat("Recursive? "); str(features_thes[tvar]) }
  invisible(x = NULL)
}
check_genes = c(
  "AMICA|JAML",
  "CXXC11", # looks like the first one could be right
  "DUSP27|DUSP29|STYXL2", # synonyms not found across datasets... wait!
  "DUPD1$|DUSP27", # shared synonym??!1
  "CCBP2", # looks like the first one could be right
  "AGPAT9", # conflictive... I find all combinations
  "GABPB2", # GABPB2 > GABPB1 > GABPAP
  "LUST", # LUST > RBM5-AS1 > SEMA3F-AS1
  "POLQ" # POLQ > POLK > TENT4A / POLQ ^ POLK in the same dataset
)
check_gene(gpattern = check_genes[4], features_names_list = features_names, features_thes = features_a_syn)
grep("DUSP27|DUSP29|STYXL2", unlist(features_a_syn), value = TRUE) # shared synonyms!!

# checking if there's recursivity
recursive_names = c()
for(i in names(features_a_syn)){ # for each 'original' gene name [from datasets]
  synl = features_a_syn[grep(paste0("^", i, "$"), names(features_a_syn))] # exact name
  find_them = unique(c( # take all possible names
    unlist(lapply(synl, head, n = Inf)), # only the first one in each thesaurus (app., it doesn't matter)
    names(synl), # also carry the name of the found thesaurus (the same as just i if no regex was used)
    unlist(strsplit(i, "\\|")) # separate the regex'd names
  ))
  gpattern_ext = paste0(paste0("^", find_them, "$"), collapse = "|") # exact matches
  # look for ALL collected names in the 'original' gene names (thesaurus list for all genes)
  j <- features_a_syn[grep(gpattern_ext, names(features_a_syn))]
  if(length(j) > 1) recursive_names <- c(recursive_names, i)
}; tvar = length(unique(recursive_names))
cat("Recursive names: ", length(recursive_names), " (unique, ", tvar, ")\n", sep = "")
head(recursive_names); tail(recursive_names)
cat(stringr::str_wrap(paste0(recursive_names, collapse = ", "), width = 70), "\n")
try2save = grep("FABP5|MS4", recursive_names, value = TRUE)
features_a_syn[try2save]
sapply(try2save, function(x) sapply(features_names, function(y) x %in% y ) )

# Now don't take the conflictive ones for homogenisation...
features_a_syn_clean0 <- features_a_syn[!names(features_a_syn) %in% recursive_names]
length(features_a_syn_clean0)
features_a_syn_clean = lapply( # eliminate shared synonyms
  X = setNames(seq(length(features_a_syn_clean0)), names(features_a_syn_clean0)),
  FUN = function(i){
    y <- features_a_syn_clean0[[i]]
    if(length(y) == 1) return(y) # cool, it's just one
    z <- y %in% unlist(features_a_syn_clean0[-i]) # discarding the ones found as other's synonym
    return(y[!z])
})
affected = sapply(features_a_syn_clean, length) != sapply(features_a_syn_clean0, length)
features_a_syn_clean0[affected]
features_a_syn_clean[affected]
table(sapply(features_a_syn_clean, length)) # CCBP2 has 0 now
features_a_syn_clean[sapply(features_a_syn_clean, length) == 0]
features_a_syn_clean[sapply(features_a_syn_clean, length) == 2]
features_a_syn_clean <- features_a_syn_clean[sapply(features_a_syn_clean, length)>0]
features_a_syn_clean = lapply(features_a_syn_clean, head, 1) # only take one synonym

features_names_update = lapply(
  X = setNames(seq(length(features_names)), names(features_names)),
  FUN = function(i){
    cat(names(features_names)[i], "\n"); genes_i <- features_names[[i]]
    synl = features_a_syn_clean[names(features_a_syn_clean) %in% genes_i] # take only genes in set
    cat(length(synl), "/", length(features_a_syn_clean), "thesaurus'd genes in set\n")
    cat("Replacing\n"); updated_q = 0
    for(j in names(synl)){
      # eliminate already present synonyms so it doesn't duplicate
      y <- synl[[j]][!synl[[j]] %in% genes_i]
      if(length(y) > 0){ # gene can be updated, do it!
        genes_i[which(genes_i %in% j)] <- y[[1]]; updated_q <- updated_q + 1
      }
    }; cat(updated_q, "updated\n\n")
    if(updated_q > length(synl)) warning("More replacements than found genes?")
    if(sum(duplicated(genes_i))) warning("duplicates created in", i, "\n")
    return(genes_i)
})
str(features_names)
str(features_names_update)
check_gene(gpattern = check_genes[4], features_names_list = features_names, features_thes = features_a_syn_clean)
check_gene(gpattern = check_genes[4], features_names_list = features_names_update, features_thes = features_a_syn_clean)

# sapply(features_names, function(y) sum(grepl("ORF", y)) )
# sapply(features_names, function(y) sum(grepl("orf", y)) )
features_clean_shared0 = sapply(features_a_syn_clean, function(x) sapply(features_names, function(y) x %in% y ) )
features_clean_shared = sapply(features_a_syn_clean, function(x) sapply(features_names_update, function(y) x %in% y ) )
if(1){ # 354
  cat("How many genes are actually across all datasets?\n")
  cat("Before: "); str(colnames(features_clean_shared0)[apply(X = features_clean_shared0, MARGIN = 2, FUN = all)])
  cat("After: "); str(colnames(features_clean_shared)[apply(X = features_clean_shared, MARGIN = 2, FUN = all)])
}
tvar = "JAML|AMICA"
features_clean_shared0[, grep(tvar, colnames(features_clean_shared0))]
features_clean_shared[, grep(tvar, colnames(features_clean_shared))]
features_a_syn_clean[grep(tvar, names(features_a_syn_clean))]

cat("Same length?")
reshape2::melt(sapply(features_names, length) == sapply(features_names_update, length))
str(features_names_update)
saveRDS(features_names_update, file = "harmonise_features/features_updated.rds")


mat0 = sapply(X = features_names, FUN = function(x){
  sapply(X = features_names, FUN = function(y) mean(x %in% y) )
}); mat0 <- round(mat0 * 100, 2);
mat = sapply(X = features_names_update, FUN = function(x){
  sapply(X = features_names_update, FUN = function(y) mean(x %in% y) )
}); mat <- round(mat * 100, 2);
matdiff = mat - mat0
hcr = hclust(dist(mat))
hcc = hclust(dist(t(mat)))
annor = data.frame(Size = sapply(X = features_names, FUN = length))
annor$Size = format(annor$Size, nsmall = 1, big.mark = ",")
anncols = list(Size = setNames(RColorBrewer::brewer.pal(n = nrow(annor), name = "Reds"), sort(annor$Size)))

pdf(paste0("harmonise_features/overlap_gene_names.pdf"))
hm <- pheatmap::pheatmap(
  mat = mat0[hcr$labels[hcr$order], hcc$labels[hcc$order]],
  color = colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "YlGnBu"))(100),
  angle_col = 90, annotation_row = annor, annotation_colors = anncols,
  cluster_rows = FALSE, cluster_cols = FALSE
); #hm$tree_row$labels[hm$tree_row$order]
pheatmap::pheatmap(
  mat = mat[hcr$labels[hcr$order], hcc$labels[hcc$order]],
  color = colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "YlGnBu"))(100),
  angle_col = 90, annotation_row = annor, annotation_colors = anncols,
  cluster_rows = FALSE, cluster_cols = FALSE
)
pheatmap::pheatmap(
  mat = matdiff[hcr$labels[hcr$order], hcc$labels[hcc$order]],
  color = colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Greys"))(100),
  angle_col = 90, annotation_row = annor, annotation_colors = anncols,
  cluster_rows = FALSE, cluster_cols = FALSE
)
graphics.off()

fnames = list.files(pattern="feature.*rds|gene_")
fnames
file.rename(fnames, paste0("harmonise_features/", fnames))
