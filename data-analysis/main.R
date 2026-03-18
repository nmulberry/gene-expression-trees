library(tidyverse)
library(RColorBrewer)
require(data.table)
library(Seurat)
library(Matrix)
library(ape)
library(ggtree)
library(TreeDist)
library(phangorn)
library(phytools)
library(ggplot2)
theme_set(theme_bw())
################################
################################
# set-up
tree_path <- "gastruloid-data/BEAST2_runs/"
data_path <- "gastruloid-data/Data/filtered_feature_bc_matrix/"

## Load 
raw_counts <- readMM(paste0(data_path, "matrix.mtx"))
barcodes <- read.table(paste0(data_path, "barcodes.tsv"), stringsAsFactors=F)[,1]
features <- read.csv(paste0(data_path, "features.tsv"), stringsAsFactors=F, sep="\t", header=F)
counts <- raw_counts
##combine count data with corresponding features/barcode annotations
rownames(counts) <- features[,1]
colnames(counts) <- barcodes
# remove tape lead
counts <- counts[rownames(counts) != "TAPE-lead",]
#----------------------#

## generate trees
clade <- 0 # all 


## get parsimony scores

#cell state annotations
states <- read.csv("gastruloid-data/mGASv2_Lane2_Group1_cell_annotation.csv")



#------ options --------#
# run over all options
df_runs <- crossing(preprocessing = c("sanity", "sctransform", "raw"),
    clade = c(0, 340),
    genes = c("all". "gemli")
)

#---- run ----#
# generate trees

mcc_tree <- ape::read.nexus(paste0(tree_path,"4-mGASv2-skyline-ou-40K.1000Kresampled.MCC.CommonAncestorHeights.tree"))
ccd_tree <- ape::read.nexus(paste0(tree_path,"4-mGASv2-skyline-ou-40K.1000Kresampled.CCD0.CommonAncestorHeights.tree"))



# get parsimony scores


