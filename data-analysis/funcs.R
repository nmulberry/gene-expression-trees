###################################
#--Helper funcs for analysing data
###################################

write_tree_from_data <- function(dat, dir, fname){
    dist_mat <- dist(t(dat))
    tree <- upgma(dist_mat)
    write.tree(tree, paste0(dir,fname, ".nwk")) 
}


get_K <- function(g_vec, tree){
	# calc Bloomberg's K for each trait on tree
    res <- phylosig(tree, g_vec, method="K")
	return(res[1])
}


get_K_cutoff <- function(tr_dat, tree, cutoff=0.8){
    res <- apply(tr_dat, 1, get_K, tree=tree)
    res <- data.frame(K=res, gene=rownames(tr_dat))
	#apply cutoff 
    return(filter(res, K > cutoff)$gene)
}

#helper func to call parsimony score

score_over_params <- function(genes, norm, dat, tree_id, reference_tree, parsimony_dir, xml_path, beast_path){
    #--write tree
	tree_file <- paste0("tree_", genes, "_", norm, "_", tree_id)
	write_tree_from_data(dat, parsimony_dir, tree_file)
	# also get PI distance to reference tree (before running parsimony)
	new_tree <- read.tree(paste0(parsimony_dir, tree_file, ".nwk"))
	tree_dist <- TreeDistance(new_tree, reference_tree)
	#--get score
	score <- get_parsimony_score(tree_file, parsimony_dir, xml_path, beast_path)
	#--return df
    return(data.frame(norm=norm, genes=genes, tree_id=tree_id, parsimony_score=score,
        num_genes_used = nrow(dat), tree_dist = tree_dist))
}

# main parsimony func
get_parsimony_score <- function(id, parsimony_dir, xml_path, beast_path){
	# read in dat from xml
	doc <- read_xml(paste0(parsimony_dir, xml_path))
	tree_file <- paste0(parsimony_dir, id, ".nwk")
    tree <- read.tree(tree_file)

	tree_sequences <- tree$tip.label 

    all_sequences <- xml_find_all(doc, ".//sequence")

	seq_info <- data.frame(
  		id = xml_attr(all_sequences, "id"),
  		taxon = xml_attr(all_sequences, "taxon"),
  		stringsAsFactors = FALSE
	)

	# get cell ids
	seq_info <- seq_info %>%
		mutate(cell_id = str_extract(id, "(?<=cell_)[A-Z]+(?=-1)"))


	# match w/ tree
	not_in_tree <- setdiff(seq_info$cell_id, tree_sequences)
	# also match taxa in tree
	new_labels <- seq_info$taxon[match(tree$tip.label, seq_info$cell_id)]

	#now remove extra cells
	for (seq_node in all_sequences) {
 	 	current_taxon <- xml_attr(seq_node, "taxon")
  
  		if (!(current_taxon %in% new_labels)) {
    		xml_remove(seq_node)
  		}
	}
	# change log file too
	logger_node <- xml_find_first(doc, ".//logger[@id='tracelog']")	
	log_file <- paste0(parsimony_dir, "parsimony_score_", id, ".log")	
	xml_set_attr(logger_node, "fileName", log_file)

	# save
	tree$tip.label <- new_labels

	xml_file <- paste0(parsimony_dir, "filtered_parsimony_", id, ".xml")

	write.tree(tree, file=tree_file)
	write_xml(doc, xml_file)

	# RUN BEAST
	tree_str <- trimws(readChar(tree_file, file.info(tree_file)$size))
	newick_tree <- paste0("newick-tree=",shQuote(tree_str))
	system2(beast_path, args=c("-overwrite -D", newick_tree, xml_file))		

	# read log file
	res <- read.table(log_file, sep="\t", header=T)
	return(res$parsimony)
}

get_tree_dist<- function(tree_id, reference_tree_id, dir){
    # read in trees
    tree <- read.tree(paste0(dir, "tree_", tree_id, ".nwk"))
    ref_tree <- read.tree(paste0(dir, "tree_", reference_tree_id, ".nwk"))
    # get PI dist
    return(TreeDistance(tree, ref_tree))
}





