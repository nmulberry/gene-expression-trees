library(tidyverse)
library(ape)
library(scTreeSim)
library(phangorn)
# VARY: ALPHA, SIGMA, NUM GENERATIONS, K (SHAPE PARAM) 

args = commandArgs(trailingOnly=TRUE)

set.seed(as.numeric(args[1]))
message(paste("running with seed", args[1]))

generate_sync_tree <- function(ng, k, rho){
#    res <- sim_adb_origin_samp(1.0, a=1/(k*ng), b=k, d=0.0, rho=rho)
    res <- sim_adb_ntaxa_samp(ntaxa=rho*(2^ng), a=1/(k*ng), b=k, d=0, rho=rho)
    tree <- res@phylo
    return(tree)
}

get_dist <- function(iter, num_sites, tree, sd, alpha, pc, shape_gamma, mu){
    # ----CHOOSE MODEL---
    if (alpha == 0.0){
        res <- sapply(rep(0, num_sites), function(x) {
	        r <- 10
		    sim_trait <- rTraitCont(tree, model="BM",
		        sigma=sd, root.value=r)
            val <- abs(sim_trait)
           # gamma-distributed error
            measured_val <- sapply(val, function (x) rgamma(1, shape=shape_gamma,scale=pc*x/shape_gamma)) 
            return(measured_val)})
	} else {
        res <- sapply(rep(0, num_sites), function(x) {
		    r <- 10
            #draw random optimum for each trait
      
		    sim_trait <- rTraitCont(tree, model="OU",
		    sigma=sd, alpha=alpha, theta=mu,root.value=r)
            val <- abs(sim_trait)
            measured_val <- sapply(val, function (x) rgamma(1, shape=shape_gamma,scale=pc*x/shape_gamma)) 
        return(measured_val)})
    }
    #--- create distance matrix 
	dist <- as.matrix(dist(res))
	recreated_tree <- upgma(dist)
	dist_topo <- RF.dist(tree, recreated_tree, normalize=T)
	return(dist_topo)
}


get_accuracy <- function(niter,ng, alpha, sigma, kappa, pc, shape_gamma,mu, rho=1) {
	# Get tree
	tree <- generate_sync_tree(ng, kappa, rho)
	# run sims over tree
    M_vec <- ceiling(10^seq(2,4, by=0.1))
	
    accuracy <- 0*M_vec
    med_RF <- 0*M_vec
    up_RF <- 0*M_vec
    low_RF <- 0*M_vec
    it <- 1
    
    for (m in M_vec){
	res <- sapply(1:niter, get_dist, m, tree, sigma, alpha, pc, shape_gamma,mu)
	accuracy[it] <- sum(res==0)/niter
        med_RF[it] <- quantile(res, prob=0.5)
        up_RF[it] <- quantile(res, prob=0.95)
        low_RF[it] <- quantile(res, prob=0.05)
        it <- it+1 
	}
	message(paste("done:",ng, alpha, sigma, kappa, pc, max(accuracy), "\n"))
	return(data.frame(mu=mu, accuracy=accuracy,median=med_RF, upper=up_RF, lower=low_RF, ng=ng, M=M_vec, alpha=alpha,sigma=sigma, kappa=kappa, rho=rho, pc=pc, shape=shape_gamma))
}
#--------------------#
# add drop-outs and 
## test
alpha <- c(0.0)
sigma <- c(0.1, 0.3, 1.0)
kappa <- c(500)
niter <- 50
pc <- c(0.1, 0.2)
shape_gamma <- c(5,10,20,50,100)
ng <- c(3,5,7,9)
mu <- c(100)
pars <- crossing(niter=niter, ng=ng, alpha=alpha, sigma=sigma, kappa=kappa, pc=pc, shape_gamma=shape_gamma, mu=mu)

df <- pars %>% pmap_dfr(., get_accuracy)

saveRDS(df, paste0("sim_res_sampling_", args[1], ".RDS"))
