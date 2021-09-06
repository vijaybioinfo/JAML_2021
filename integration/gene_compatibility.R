#!/usr/bin/R

################################
# Seurat integration: cleaning #
################################

# This scripts selects the genes that will be integrated for the AMICA paper

output_dir = "/home/ciro/large/amica/results/integration/"
# dir.create(output_dir); setwd(output_dir)
dir.create("harmonise_features/")
fname = list.files(path = output_dir, pattern = "cells_selected", full.names = TRUE)
metadata_list_subset = readRDS(file = fname)
sufix = paste0(length(metadata_list_subset), "sets")

### Loading packages ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(tidyverse)
})
theme_set(theme_cowplot())
source("/home/ciro/scripts/handy_functions/devel/file_reading.R")
source("/home/ciro/scripts/handy_functions/devel/filters.R")
source("/home/ciro/scripts/handy_functions/devel/utilities.R")

### Loading expression data ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patty <- paste0(names(metadata_list_subset), collapse = "|")
edata_list_f <- list.files(path = "/home/ciro/large/simon/raw", pattern = patty, full.names = TRUE)
edata_list_f <- grep(pattern = "\\.", x = edata_list_f, value = TRUE)
edata_list_f <- grep(pattern = "RData", x = edata_list_f, value = TRUE, invert = TRUE)
edata_list_f <- grep(pattern = "annot", x = edata_list_f, value = TRUE, invert = TRUE)

tvar <- "_all.*|_Tc.*|_cd.*|_10x.*|_tpm.*"
set_names <- gsub(tvar, "", basename(edata_list_f), ignore.case = TRUE)
names(edata_list_f) <- set_names
edata_list <- lapply(
  X = names(edata_list_f),
  FUN = function(i){
    y <- readfile(myfile = edata_list_f[[i]])
    y[, metadata_list_subset[[i]]]
  }
)
names(edata_list) <- names(edata_list_f)
all_genes <- lapply(X = edata_list, FUN = rownames)
all_genes_total <- unique(x = unlist(x = all_genes))
str(all_genes); length(x = all_genes_total)
# saveRDS(object = all_genes, file = paste0("total_genes_", sufix, ".rds"))
gene_names_df = reshape2::melt(lapply(all_genes, function(x) grep(pattern = "JAML|AMICA", x, value = TRUE) ))

### Checking annotations ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## checking if I can 'translate' ## --------------------------------------------
setwd(tempdir())
# system("grep -E 'JAML|AMICA' /mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCh38P13_annotation.csv")
# system("grep -E 'JAML|AMICA' /mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCh37_annotation.csv")
fnames = c(
  hg19_ucsc = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
  hg38_ucsc = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",
  GRCh38v103_ensembl = "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz",
  GRCh37v87_ensembl = "ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz",
  v27_gencode = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz",
  v32_gencode = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz",
  v19_gencode = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
  GRCh37_ncbi = "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20201022/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz",
  GRCh38_ncbi = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/current/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz",
  hg19_cellranger = "/mnt/BioAdHoc/Groups/vd-vijay/references/refdata-cellranger-hg19-3.0.0/genes/genes.gtf", # GRCh37v87_ensembl
  GRCh37_ndu = "/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCh37_annotation.csv", # GENCODE Release 19
  GRCh38_ndu = "/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCh38P13_annotation.csv"
)
gene_annotated = lapply(
  X = names(fnames),
  FUN = function(i){
    new_name = paste0(i, gsub(".*(txt.*)|.*(gtf.*)", ".\\1\\2", fnames[i]))
    cat(new_name)
    fname <- if(file.exists(new_name)){
      new_name
    }else if(!file.exists(fnames[[i]])){
      system(paste0("wget -c ", fnames[[i]], " -O ", new_name)) # options: continue Output
      new_name
    }else{
      fnames[[i]]
    }
    genes <- if(grepl("gtf", fname)){
      cat(" - rtracklayer::import\n")
      annot <- rtracklayer::import(con = fname)
      if(is.null(annot$gene)) annot$gene_name else annot$gene
    }else{
      cat(" - data.table::fread\n")
      annot <- data.table::fread(input = fname)
      if(is.null(annot$V13)) annot$gene_name else annot$V13
    }
  }
); names(gene_annotated) <- names(fnames)
gene_annotated <- lapply(X = gene_annotated, FUN = unique)
str(gene_annotated)
sapply(X = gene_annotated, FUN = length) > length(x = all_genes_total)
gene_names_pct <- sapply(
  X = gene_annotated,
  FUN = function(gene_annot){
    round(sapply(X = all_genes, FUN = function(x){
      # mean(x %in% gene_annot)
      mean(gsub("\\-1$", "", x) %in% gene_annot)
    }) * 100, 1)
  }
)
which_name <-  sapply(X = gene_annotated, FUN = function(x){
  if(any(grepl("JAML", x))){ "JAML" }else if(any(grepl("AMICA", x))){ "AMICA" }else{ "NONE" }
})
colnames(gene_names_pct) <- paste0(colnames(gene_names_pct), ";", which_name)
gene_names_pct_df <- cbind(gene_names_df, gene_names_pct)
fname <- paste0("harmonise_features/gene_annotation_", sufix, ".csv")
write.csv(gene_names_pct_df, file = fname, row.names = FALSE)
gene_total_annot = reshape2::melt(lapply(
  X = gene_annotated,
  FUN = function(gene_annot){
    y <- c(mean(all_genes_total %in% gene_annot), mean(gsub("\\-1$", "", all_genes_total) %in% gene_annot))
    list(as_is = round(y[1] * 100, 1), rm1 = round(y[2] * 100, 1))
  }
))
gene_total_annot$gene_name = which_name[gene_total_annot$L1]
fname <- paste0("harmonise_features/gene_total_", sufix, ".csv")
write.csv(gene_total_annot, file = fname, row.names = FALSE)

# str(all_genes$sadefeldman_GSE120575[!all_genes$sadefeldman_GSE120575 %in% gene_annotated$GRCh37_ndu])
# str(all_genes$oh_GSE149652[!all_genes$oh_GSE149652 %in% gene_annotated$hg19_cellranger])
# str(all_genes$oh_GSE149652[!gsub("\\-1$", "", all_genes$oh_GSE149652) %in% gene_annotated$hg19_cellranger])
# str(all_genes$oh_GSE149652[!gsub("\\-[0-9]{1,}$", "", all_genes$oh_GSE149652) %in% gene_annotated$hg19_cellranger])
