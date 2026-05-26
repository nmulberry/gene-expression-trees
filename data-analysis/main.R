library(tidyverse)
library(xml2)
library(RColorBrewer)
require(data.table)
library(Seurat)
library(Matrix)
library(ape)
library(scran)
library(ggtree)
library(TreeDist)
library(phangorn)
library(phytools)
library(ggplot2)
library(SingleCellExperiment)
library(scran)
library(scater)
library(stats)
library(patchwork)
library(scTreeSim)
theme_set(theme_bw())

source("funcs.R")
################################
################################
# set-up
tree_path <- "gastruloid-data/BEAST2_runs/"
data_path <- "gastruloid-data/Data/filtered_feature_bc_matrix/"
parsimony_dir <- "parsimony-xml/"
xml_path <- "4-mGASv2-skyline-ou-40K-parsimony.xml" #original file 
beast_path <- "/Applications/BEAST 2.7.7/bin/beast"


## Load 
raw_counts <- readMM(paste0(data_path, "matrix.mtx"))
barcodes <- read.table(paste0(data_path, "barcodes.tsv"), stringsAsFactors=F)[,1]
features <- read.csv(paste0(data_path, "features.tsv"), stringsAsFactors=F, sep="\t", header=F)
counts <- raw_counts

##combine count data with corresponding features/barcode annotations
rownames(counts) <- features[,1]
colnames(counts) <- barcodes
# remove tape lead
counts <- counts[rownames(counts) != "TAPE_lead",]
colnames(counts) <-  sub("-1$", "", colnames(counts))

colnames(raw_counts) <- barcodes
colnames(raw_counts) <-  sub("-1$", "", colnames(raw_counts))



ccd_tree <- ape::read.nexus(paste0(tree_path,"4-mGASv2-skyline-ou-40K.1000Kresampled.CCD0.CommonAncestorHeights.tree"))
states <- read.csv("gastruloid-data/mGASv2_Lane2_Group1_cell_annotation.csv")%>% dplyr::select(-X)
#--subsample

label_ids <- read.csv(paste0(tree_path, "mGASv2_cellbc_tiplabel.csv"))
ntip <- length(ccd_tree$tip.label)
label_ids$Cell <- sub("-1$", "", label_ids$cell_bc)

map_labels <- setNames(label_ids$Cell, label_ids$tip_label)
ccd_tree$tip.label <- map_labels[ccd_tree$tip.label]

#--keep only annotated cells
duplicated <- states$Cell[duplicated(states$Cell)]
to_keep <- intersect(intersect(states$Cell, label_ids$Cell), colnames(counts))
to_keep <- setdiff(to_keep, duplicated)

sm_ccd_tree <- ape::keep.tip(ccd_tree, to_keep)
states <- filter(states, Cell %in% to_keep)
counts <- counts[, colnames(counts) %in% to_keep]


#----------------------#
# normalization
#----------------------#

# scTransform
seurat_obj <- CreateSeuratObject(counts)
seurat_obj <- SCTransform(seurat_obj)
seurat_counts <- seurat_obj[["SCT"]]$scale.data

# load in sanity results
sanity_counts <- data.table::fread("gastruloid-data/sanity-analysis/log_transcription_quotients.txt", sep="\t",
    header=T, data.table=F)
colnames(sanity_counts) <- colnames(raw_counts)
sanity_counts <- sanity_counts[, colnames(sanity_counts) %in% to_keep]


# deconvolution method from scran
sce <- SingleCellExperiment(assays = list(counts = counts))
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce)
scran_counts <- logcounts(sce)
#----------------------#
# TREE VISUALIZATION
#----------------------#
## generate trees
clade <- 0 # all 

types_to_keep <- c("ExE endoderm","ExE mesoderm","Forebrain/Midbrain/Hindbrain", "Pharyngeal mesoderm",
	"Endothelium", "Intermediate mesoderm", "Surface ectoderm", "Parietal endoderm", "Rostral neurectoderm")
states <- states %>% mutate(
	types = case_when(pijuan_celltype %in% c(types_to_keep) ~ pijuan_celltype,
						TRUE ~ "Other"))
#cells_for_plotting <- unique(filter(states, types %in% types_to_keep)$Cell) 
cells_for_plotting <- unique(states$Cell)
rownames(states) <- states$Cell

# visualize
gg_tree <- ggtree(keep.tip(sm_ccd_tree, cells_for_plotting), layout="circ", size=0.1)
gg_tree <- gheatmap(gg_tree, states["types"], offset=0, width=0.05, colnames=F)+
	scale_fill_brewer(palette="Spectral")+labs(fill="")+
	ggtitle("Sciphy CCD Tree")


raw_tree <- upgma(dist(t(counts)))

gg_raw_tree <- ggtree(keep.tip(raw_tree, cells_for_plotting), layout="circ", size=0.1)
gg_raw_tree <- gheatmap(gg_raw_tree, states["types"], offset=0, width=0.05, colnames=F)+
	scale_fill_brewer(palette="Spectral")+labs(fill="")+
	ggtitle("Raw")

sct_tree <- upgma(dist(t(seurat_counts)))

gg_sct_tree <- ggtree(keep.tip(sct_tree, cells_for_plotting), layout="circ", size=0.1)
gg_sct_tree <- gheatmap(gg_sct_tree, states["types"], offset=0, width=0.05, colnames=F)+
	scale_fill_brewer(palette="Spectral")+labs(fill="")+
	ggtitle("scTransform")


scran_tree <- upgma(dist(t(scran_counts)))

gg_scran_tree <- ggtree(keep.tip(scran_tree, cells_for_plotting), layout="circ", size=0.1)
gg_scran_tree <- gheatmap(gg_scran_tree, states["types"], offset=0, width=0.05, colnames=F)+
	scale_fill_brewer(palette="Spectral")+labs(fill="")+
	ggtitle("Deconvolution")


sanity_tree <- upgma(dist(t(sanity_counts)))

gg_sanity_tree <- ggtree(keep.tip(sanity_tree, cells_for_plotting), layout="circ", size=0.1)
gg_sanity_tree <- gheatmap(gg_sanity_tree, states["types"], offset=0, width=0.05, colnames=F)+
	scale_fill_brewer(palette="Spectral")+labs(fill="")+
	ggtitle("Sanity")


(gg_tree | gg_raw_tree)+plot_layout(guides="collect")&theme(legend.position="bottom")
my_layout <- "
AB
CG
"
gg_sanity_tree + gg_sct_tree + gg_scran_tree + guide_area()+
	plot_layout(design=my_layout,guides="collect") & theme(legend.position="right")



#------------------------------------------------------------------#
# Within-type trees
#------------------------------------------------------------------#

# get cells by types w. > 15 cells
type <- "Pharyngeal mesoderm"
to_keep <- intersect(filter(states, pijuan_celltype == type)$Cell, label_ids$Cell)
Pm_tree <- ape::keep.tip(sm_ccd_tree, to_keep)

type <- "ExE ectoderm"
to_keep <- intersect(filter(states, pijuan_celltype == type)$Cell, label_ids$Cell)
Exeecto_tree <- ape::keep.tip(sm_ccd_tree, to_keep)

type <- "ExE endoderm"
to_keep <- intersect(filter(states, pijuan_celltype == type)$Cell, label_ids$Cell)
Exeendo_tree <- ape::keep.tip(sm_ccd_tree, to_keep)

type <- "Forebrain/Midbrain/Hindbrain"
to_keep <- intersect(filter(states, pijuan_celltype == type)$Cell, label_ids$Cell)
fbm_tree <- ape::keep.tip(sm_ccd_tree, to_keep)



#------------------------------------------------------------------#
# GET TREES
#------------------------------------------------------------------#
ccd_clade <- ccd_clade <- ape::extract.clade(sm_ccd_tree, 340) #"Mixed" clade w. 14 tips 

reference_trees <- c(sm_ccd_tree, ccd_clade, Pm_tree, Exeecto_tree, Exeendo_tree, fbm_tree)
tree_ids <- c("all", "clade-340", "Pm", "ExEecto", "ExEendo", "FBM") 

# initialize
scT_counts0 <- seurat_counts
sanity_counts0 <- sanity_counts
scran_counts0 <- scran_counts
counts0 <- counts

 

run_over_cells <- function(tree, tree_id, with_K){

    to_keep <- tree$tip.label
	ntip <- length(to_keep)
    # filter counts too
    sanity_counts <- sanity_counts0[, colnames(sanity_counts0) %in% to_keep]
    scT_counts <- scT_counts0[, colnames(scT_counts0) %in% to_keep]
    counts <- counts0[, colnames(counts0) %in% to_keep]
    scran_counts <- scran_counts0[, colnames(scran_counts0) %in% to_keep]

    #---- Get K threshold genes
    if (with_K){
    print("Running K thresholds on seurat data")
    K_genes <- get_K_cutoff(scT_counts, tree, cutoff=0.95)
    K_scT_counts <- scT_counts[rownames(scT_counts) %in% K_genes,]
    print("Running K thresholds on sanity data")
    K_genes <- get_K_cutoff(sanity_counts, tree, cutoff=0.95)
    K_sanity_counts <- sanity_counts[rownames(sanity_counts) %in% K_genes,]
    print("Running K thresholds on raw counts")
    K_genes <- get_K_cutoff(counts, tree, cutoff=0.95)
    K_raw_counts <- counts[rownames(counts) %in% K_genes,]
    print("Running K thresholds on scran counts")
    K_genes <- get_K_cutoff(scran_counts, tree, cutoff=0.95)
    K_scran_counts <- scran_counts[rownames(scran_counts) %in% K_genes,]
    }
    #------------------------------------------------------------------#
    # GENERATE TREES FROM DIFFERENT DATA SOURCES 
    #------------------------------------------------------------------#
    print(paste("Running parsimony score for", tree_id))
    # using all genes
	score1 <- score_over_params("all", "scT", scT_counts, tree_id, tree, parsimony_dir, xml_path, beast_path)  
    score2 <- score_over_params("all", "sanity", sanity_counts, tree_id, tree, parsimony_dir, xml_path, beast_path)  
    score3 <- score_over_params("all", "raw", counts, tree_id, tree, parsimony_dir, xml_path, beast_path)  
    score4 <- score_over_params("all", "deconvolve", scran_counts, tree_id, tree, parsimony_dir, xml_path, beast_path)  
    res <- rbind(score1, rbind(score2, rbind(score3, score4))) 
    if (with_K){
    # using K thresh
    score1 <- score_over_params("K", "scT", K_scT_counts, tree_id, tree, parsimony_dir, xml_path, beast_path)  
    score2 <- score_over_params("K", "sanity", K_sanity_counts, tree_id, tree, parsimony_dir, xml_path, beast_path)  
    score3 <- score_over_params("K", "raw", K_raw_counts, tree_id, tree, parsimony_dir, xml_path, beast_path)  
    score4 <- score_over_params("K", "deconvolve", K_scran_counts, tree_id, tree, parsimony_dir, xml_path, beast_path)  
    res <- rbind(res, rbind(score1, rbind(score2, rbind(score3, score4))))
    }
    #------------------------------------------------------------------#
	
	# generate tree & randomly assign tip labels
	# note: only topology matters
	n <- 10
	random_score <- rep(0,n)
	random_dist <- rep(0,n)
	samp <- ntip/2^15
	rtree <- sim_adb_ntaxa_samp(ntip, 1, 20, 0.0, samp)@phylo
	for (i in 1:n){
		rtree$tip.label <- sample(unname(tree$tip.label), replace=F)	
		write.tree(rtree, paste0(parsimony_dir, "tree_barcodes_random_", tree_id, ".nwk"))
		random_score[i] <- get_parsimony_score(paste0("tree_barcodes_random_",tree_id), parsimony_dir, xml_path, beast_path)
		random_dist[i] <- TreeDistance(rtree, tree)
	}
	random_res <- data.frame(norm="random", tree_id=tree_id, genes="none", parsimony_score =
		median(random_score), num_genes_used = 0, tree_dist = median(random_dist))
	res <- rbind(res, random_res)
	
	# sciphy reference score
    write.tree(tree, paste0(parsimony_dir, "tree_barcodes_reference_", tree_id, ".nwk"))
    reference_score <- get_parsimony_score(paste0("tree_barcodes_reference_",tree_id), parsimony_dir, xml_path, beast_path)
	
	res$relative_parsimony_score <- (res$parsimony_score - reference_score)/reference_score

    res$tree_size <- length(to_keep)

    return(res)
}


res <- mapply(run_over_cells, reference_trees, tree_ids, with_K=T, SIMPLIFY=F)
res <- do.call("rbind", res)


res <- res %>%
    mutate(label=case_when(tree_id == "all" ~ "Whole Tree (285 cells)",
        tree_id == "clade-340" ~ "Mixed Clade (14 cells)",
        tree_id == "Pm" ~ "Pharyngeal Mesoderm (45 cells)",
        tree_id == "ExEecto" ~ "ExE Ectoderm (13 cells)",
        tree_id == "ExEendo" ~ "ExE Endoderm (69 cells)",
        tree_id == "FBM" ~ "Forebrain/Midbrain/Hindbrain (15 cells)"))

#---plot
res$norm <- factor(res$norm, levels=c("random", "raw", "scT", "deconvolve", "sanity"))

gg <- ggplot(res, aes(y=relative_parsimony_score, x=norm, fill=genes))+
    geom_bar(stat="identity", position=position_dodge2(width=.5, preserve="single"), col="black", linewidth=0.75, width=.5)+
    facet_wrap(~label)+
    scale_fill_brewer(palette="Blues")+
    scale_y_continuous(expand=c(0,0,0.05, 0))+
    labs(x="Method", y="Relative Parsimony Score", fill="Data")
#ggsave("parsimony-scores.pdf", height=6, width=12)    

gg2 <- ggplot(res, aes(y=tree_dist, x=norm, fill=genes))+
    geom_bar(stat="identity", position=position_dodge2(width=.5, preserve="single"), col="black", linewidth=0.75, width=0.5)+
    facet_wrap(~label)+
    scale_fill_brewer(palette="Greens")+
    scale_y_continuous(expand=c(0,0,0.05, 0))+
    labs(x="Method", y="CI Distance", fill="Data")



## compare # of genes used too
res2 <- filter(res, norm != "random") %>% group_by(norm, tree_id, label) %>%
	summarize(prop_K = num_genes_used[genes=="K"]/num_genes_used[genes=="all"]) %>% 
	ungroup()

gg3 <- ggplot(res2, aes(y=prop_K, x=norm))+
    geom_bar(stat="identity", position=position_dodge2(width=.5, preserve="single"), col="black", linewidth=0.75, width=0.5)+
    facet_wrap(~label)+
    scale_y_continuous(expand=c(0,0,0.05, 0))+
    labs(x="Method", y="Proportion genes with K>0.95")











