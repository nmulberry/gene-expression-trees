library(tidyverse)
library(ape)
library(scTreeSim)
library(phangorn)
# VARY: ALPHA, SIGMA, NUM GENERATIONS, K (SHAPE PARAM) 
args = commandArgs(trailingOnly=TRUE)

set.seed(as.numeric(args[1]))
message(paste("running with seed", args[1]))

generate_sync_tree <- function(ng, k, rho){
    res <- sim_adb_origin_samp(1.0, a=1/(k*ng), b=k, d=0.0, rho=rho)
    tree <- res@phylo
    return(tree)
}

get_dist <- function(iter, num_sites, tree, sd, alpha){
	
    # ----CHOOSE MODEL---
    if (alpha == 0.0){
        res <- sapply(rep(0, num_sites), function(x) {
	    	r <- 1000
		    sim_trait <- rTraitCont(tree, model="BM",
		    sigma=sd, root.value=r)
		return(abs(sim_trait))})
	} else {

        res <- sapply(rep(0, num_sites), function(x) {
		    r <- 1000
	            #draw random optimum for each trait
        	    mu <- 1000+runif(1,10,1000)
		    sim_trait <- rTraitCont(tree, model="OU",
		    sigma=sd, alpha=alpha, theta=mu,root.value=r)
		return(abs(sim_trait))})
    }
    #--- create distance matrix 
	dist <- as.matrix(dist(res, method="manhattan"))
	recreated_tree <- upgma(dist)
	dist_topo <- dist.topo(tree, recreated_tree)
	return(dist_topo)
}


get_accuracy <- function(niter,ng, alpha, sigma, kappa, rho=1) {
	# Get tree
	tree <- generate_sync_tree(ng, kappa, rho)
	# run sims over tree
	accuracy <- 0	
    	M_vec <- ceiling(10^seq(2,3.5, by=0.15))
	accuracy <- 0*M_vec
    	it <- 1
    	for (m in M_vec){
		res <- sapply(1:niter, get_dist, m, tree, sigma, alpha)
		accuracy[it] <- sum(res==0)/niter
        it <- it+1 
	}
	message(paste("done:",ng, ceiling(m), alpha, sigma, kappa, "\n"))
	return(data.frame(accuracy=accuracy,M=M_vec, ng=ng, alpha=alpha,sigma=sigma, kappa=kappa, rho=rho))
}

## test
alpha <- c(1.0, 1.5, 2.0)
sigma <- c(0.3, 1.0)
kappa <- c(10,50,500)
niter <- 50
ng <- seq(3,9,by=2)
pars <- crossing(niter=niter, ng=ng, alpha=alpha, sigma=sigma, kappa=kappa)

df <- pars %>% pmap_dfr(., get_accuracy)
saveRDS(df, paste0("sim_res_OU_", args[1], ".RDS"))
