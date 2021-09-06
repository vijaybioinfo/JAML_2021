#!/usr/bin/bash

################################
# SiEs12 mouse data processing #
################################

# This script was used to process the CD45+JAML+ murine cells

### Metadata ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp /Users/ciro/Documents/liai/cancer/amica/SiEs12_metadata_library.csv /Volumes/ciro/amica/info/
cp /Users/ciro/Documents/liai/cancer/amica/SiEs12_fbarcodes.csv /Volumes/ciro/amica/info/
Rscript /home/ciro/scripts/functions/csvCorrect.R /home/ciro/amica/info/SiEs12_metadata_library.csv
Rscript /home/ciro/scripts/functions/csvCorrect.R /home/ciro/amica/info/SiEs12_fbarcodes.csv

### Mapping ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cp /home/ciro/preethi/scripts/tfh_cellranger.yaml /home/ciro/amica/scripts/SiEs12_cellranger.yaml
# atom --add /Volumes/ciro/amica/scripts/SiEs12_cellranger.yaml
sh /home/ciro/scripts/cellranger_wrappeR/run.sh -y /home/ciro/amica/scripts/SiEs12_cellranger.yaml -v

### Doublets/Hashtags ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUTDIR=/home/ciro/ad_hoc/amica
SNAME=SiEs12_Mo_1_CD45_3H
mkdir -p ${OUTDIR}/results/ab_capture
Rscript /home/ciro/scripts/ab_capture/demux_seurat.R \
  --edata ${OUTDIR}/raw/cellranger_SiEs12/count/${SNAME}_Gex/outs/filtered_feature_bc_matrix \
  --capture ${OUTDIR}/raw/cellranger_SiEs12/count/${SNAME}_CITE/outs/raw_feature_bc_matrix \
  --outdir ${OUTDIR}/results/ab_capture \
  --prefix ${SNAME}

### Quality Control ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cp /home/ciro/scripts/quality_control/config.yaml /home/ciro/amica/scripts/SiEs12_qc.yaml
# atom --add /Volumes/ciro/amica/scripts/SiEs12_qc.yaml
Rscript /home/ciro/scripts/quality_control/single_cell.R -y /home/ciro/amica/scripts/SiEs12_qc.yaml

### Doublets/Scrublet ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conda activate doublets
echo "---
project_name: doubdet_test
input_matrix: ${OUTDIR}/raw/cellranger_SiEs12/count/${SNAME}_Gex/outs/filtered_feature_bc_matrix
metadata: ${OUTDIR}/results/quality_control/${SNAME}/metadata_prefilter.rdata
output_dir: ${OUTDIR}/results/doublets/${SNAME}
vars:
  library: origlib
  class: orig.HT_ID.global
  dims: umap
filters:
  score: 0.4
..." > /home/ciro/amica/scripts/SiEs12_scrublet.yaml
Rscript /home/ciro/scripts/doubdet/scrublet.R -y /home/ciro/amica/scripts/SiEs12_scrublet.yaml
conda deactivate

### Clustering ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cp /home/ciro/scripts/clustering/config.yaml /home/ciro/amica/scripts/SiEs12_clustering.yaml
# atom --add /Volumes/ciro/amica/scripts/SiEs12_clustering.yaml
sh /home/ciro/scripts/clustering/run.sh -y /home/ciro/amica/scripts/SiEs12_clustering.yaml
# cp /home/ciro/amica/scripts/SiEs12_clustering.yaml /home/ciro/amica/scripts/SiEs12_clustering_filt.yaml
# atom --add /Volumes/ciro/amica/scripts/SiEs12_clustering_filt.yaml
sh /home/ciro/scripts/clustering/run.sh -y /home/ciro/amica/scripts/SiEs12_clustering_filt.yaml
