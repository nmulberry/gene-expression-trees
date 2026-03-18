library(tidyverse)
library(geomtextpath) 
library(ggplot2)
library(viridis)
library(ape)
library(phangorn)
library(wesanderson)
library(sciscales)
#theme_set(theme_minimal(base_size=12)+
#    theme(
#      panel.border = element_rect(
#        colour = "black",
#        fill = NA,
#        linewidth = 1.2
#      )))

theme_set(theme_minimal(base_size = 10)+theme(
	  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
	  strip.placement = "outside",
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "bold"))
)
######################################
p_trip_OU <- function(M,dv,du, alpha, theta){
    gu <- 1-exp(-2*alpha*(1-du))
    gv <- 1-exp(-2*alpha*(1-dv))
    num <- sqrt(gu/alpha+2*theta^2)-sqrt(gv/alpha+2*theta^2)
    denom <- sqrt((1-2/pi)*(1/alpha*(gu+gv)+4*theta^2))
    C <- sqrt(2*M/pi)
    p0 <- pnorm(C*num/denom)
    
    return(p0)
}

p_trip_BM <- function(M,dv,du, theta){
    gu <- 2*(1-du)
    gv <- 2*(1-dv)
    num <- sqrt(gu+2*theta^2)-sqrt(gv+2*theta^2)
    denom <- sqrt((1-2/pi)*(1*(gu+gv)+4*theta^2))
	C <- sqrt(2*M/pi)
    p0 <- pnorm(C*num/denom)
    
    return(p0)
}

p_trip <- function(M,dv,du,alpha,theta){
	if (alpha > 0){
		prob <- p_trip_OU(M,dv,du,alpha,theta)
	} else {
		prob <- p_trip_BM(M,dv,du,theta)
	}
	return(prob)
}

get_tree_accuracy <- function(M,ng, alpha, theta, N){
    ell <- 1/(ng+1)
	min_du <- ell
    min_dv <- 2*ell
    p0 <- p_trip(M, min_dv, min_du, alpha, theta)
    bound <- max(0, 1-choose(N,3)*(1-p0))
    return(bound)
}


get_tree_accuracy_sync <- function(M,ng,alpha,theta){
    n <- 2^ng #num tips (fully samped)
    ell <- 1/(ng+1) #conservative 
    p0 <- 0
    for (k in (1:(ng-1))){
        du <- k*ell
		dv <- du+ell
       	pk <- 1-p_trip(M, dv,du,alpha,theta)#prob not resolving split at level k
		coeff <- n*choose(n/2^k,2) 
        p0 <- p0+(coeff*pk)
    }
    return(max(0, 1-p0))
}

get_M_pred_sync <- function(ng, eps, alpha,theta, limit=1e5){
    p_eps <- function(M){
        get_tree_accuracy_sync(M,ng,alpha,theta) - eps
    }
    # use optim
    res <- uniroot(p_eps, c(1,limit))
    return(ceiling(res$root))
}




######################################
pars <- crossing(M=10^seq(3, 4.5, by=0.1),
    ng = seq(4,14,by=1),
    alpha = c(0.0,0.2,1.0),
    theta = c(0, 0.4, 0.9))

df <- pars
pars$N <- (2^pars$ng)*0.1
df$res <- pars %>% pmap_dbl(get_tree_accuracy)

clevels <- c(0.0, 5, 50, 98, 100, 200)
library(scales)
gg_contour <- ggplot(df, aes(x=ng, y=log10(M), z=res*100))+
    facet_grid(alpha~theta, labeller=label_bquote(cols=theta:.(theta), rows=alpha:.(alpha)))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    geom_contour_filled(breaks=clevels)+
    geom_textcontour(breaks=clevels)+
    scale_fill_viridis_d(option="plasma")+
    theme(legend.position="none")+
    labs(y="# Traits", x="# generations")

######################################
pars <- crossing(M=10^seq(2, 4.5, by=0.1),
    ng = seq(2,16,by=1),
    alpha = c(0.0,0.1,0.5),
    theta = c(0, 0.4, 0.9))
df_sync <- pars
df_sync$res <- pars %>% pmap_dbl(get_tree_accuracy_sync)

gg_contour_sync <- ggplot(df_sync, aes(x=ng, y=log10(M), z=100*res))+
    facet_grid(alpha~theta, labeller=label_bquote(cols=eta:.(theta), rows=alpha:.(alpha)))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    geom_contour_filled(breaks=clevels)+
#    geom_textcontour(breaks=clevels)+
    scale_fill_viridis_d(option="plasma")+
 #   theme(legend.position="none")+
    labs(y="# traits (log10)", x="# generations")


###############################################
# min m req to reach eps = 0.95

get_M_pred <- function(ng, eps, alpha, theta,method="sync",limit=1e5){
    p_eps <- function(M){
		if (method=="sync"){
        	get_tree_accuracy_sync(M,ng,alpha,theta) - eps
    	} else {
			get_tree_accuracy(M, ng, alpha, theta, 2^ng) - eps 
		}
	}
    # use optim
	root <- tryCatch(uniroot( function(x) p_eps(x), lower=1, upper=limit)$root, 
		error=function(e){warning(conditionMessage(e)); NA})
    return(ceiling(root))
}


pars <- crossing(ng = seq(2,9),
	eps = c(0.95),
    alpha = c(0,0.1,0.5),
    theta = c(0, 0.2, 0.5),
	method=c("async"))

df_m <- pars
df_m$min_m <- pars %>% pmap_dbl(get_M_pred)
#df_m$generations <- as.factor(df_m$ng)
df_m$theta <- as.factor(df_m$theta)


gg_m <- ggplot(filter(df_m,alpha==0), aes(x=ng, y=log10(min_m), col=theta, group=theta))+
	geom_point()+
#	facet_wrap(~alpha, labeller=label_bquote("Brownian Motion"), scales="free")+
	scale_color_manual(values=wes_palette("Darjeeling1"))+
	labs(y="Min # traits (log10)", x="# generations", col=expression(eta))+
#	scale_y_continuous(labels = label_scientific())+
	theme(legend.position="bottom")
#ggsave("min_m_BM.pdf", gg_m, height=5, width=8)

pars <- crossing(ng = c(2,4,6,8),
	eps = c(0.95),
    alpha = seq(0,1.5, by=0.1),
    theta = c(0, 0.2, 0.5),
	method=c("async"))

df_m <- pars
df_m$min_m <- pars %>% pmap_dbl(get_M_pred)
#df_m$generations <- as.factor(df_m$ng)
df_m$theta <- as.factor(df_m$theta)


gg_m2 <- ggplot(df_m, aes(x=alpha, y=log10(min_m), col=theta, group=theta))+
	geom_line(linewidth=1.0)+
	facet_wrap(~ng, labeller=label_bquote("generations":.(ng)))+
	scale_color_manual(values=wes_palette("Darjeeling1"))+
	labs(y="Min # traits (log10)", x=expression("Strength of selection " * alpha), col=expression(eta))+
#	scale_y_continuous(labels = label_scientific())+
	theme(legend.position="bottom")

cowplot::plot_grid(gg_m,gg_m2, labels=c("A","B"),rel_widths=c(0.45,0.55))
ggsave("min_m_comb.pdf", height=5, width=10)

#########################
# SIMULATION RESULTS
#########################
#sim_res <- readRDS("sim_study/sim_res_OU_6.RDS")
#sim_res <- readRDS("sim_res_OU.RDS")

# READ IN ALL DATA
files <- list.files("sim_study", pattern="\\.RDS", full.names=T)
res_df <- do.call(rbind,
	lapply(files, function (f) {
		df <- readRDS(f)
		tree_id <- as.numeric(sub(".*_(\\d+)\\.RDS$", "\\1", basename(f)))
		df$tree_id <- tree_id
		return(df)
	}))

# over all trees, ng=9
df_9 <- filter(res_df, ng==9)
df_9$sigma <- factor(df_9$sigma)
df_9 <- df_9 %>% group_by(M,sigma,alpha,kappa, rho)%>%
	summarize(min_acc = min(accuracy), max_acc = max(accuracy), med_acc = median(accuracy))


gg_9 <- ggplot(df_9, aes(x=log10(M), y=med_acc,fill=sigma))+
	geom_line(aes(col=sigma), linewidth=1)+
	geom_ribbon(alpha=0.5, aes(ymin=min_acc, ymax=max_acc, group=sigma))+
	facet_grid(alpha~kappa, labeller=label_bquote(cols=kappa:.(kappa), rows=alpha:.(alpha)))+
	scale_color_manual(values=c("darkgoldenrod1", "deepskyblue1"))+
	scale_fill_manual(values=c("darkgoldenrod1", "deepskyblue1"))+
	labs(fill=expression(sigma), col=expression(sigma), y="accuracy", x="# traits (log10)")+
	theme(legend.position="bottom")

######################
sim_res2 <- res_df %>%
	group_by(M,ng,alpha,sigma,kappa, rho) %>%
	summarize(accuracy=min(accuracy)) %>%
	ungroup()

gg_res1 <- ggplot(filter(sim_res2,sigma==0.3), aes(x=ng, y=log10(M), fill=100*accuracy, col=100*accuracy))+
	geom_tile()+
	facet_grid(alpha~kappa, labeller=label_bquote(cols=kappa:.(kappa), rows=alpha:.(alpha)))+
	scale_x_continuous(expand=c(0,0), breaks=c(3,5,7,9))+
	scale_y_continuous(expand=c(0,0))+
	scale_fill_viridis(option="magma")+
	scale_colour_viridis(option="magma")+
	labs(x="# generations", y="# traits (log10)", fill="accuracy", col="accuracy")

gg_res1 <- ggplot(filter(sim_res2,sigma==0.3), aes(x=ng, y=log10(M), fill=100*accuracy, col=100*accuracy))+
	geom_tile()+
	facet_grid(kappa~alpha, labeller=label_bquote(rows=kappa:.(kappa), cols=alpha:.(alpha)))+
	scale_x_continuous(expand=c(0,0), breaks=c(3,5,7,9))+
	scale_y_continuous(expand=c(0,0))+
	scale_fill_viridis(option="magma")+
	scale_colour_viridis(option="magma")+
	labs(x="# generations", y="# traits (log10)", fill="Accuracy", col="Accuracy")+
	theme(legend.position="bottom", legend.ticks = element_blank())+
	guides(fill=guide_colourbar(
		frame.colour="black",
		barheight=unit(0.4, "cm")
	))

