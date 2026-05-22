library(tidyverse)
library(geomtextpath) 
library(ggplot2)
library(viridis)
library(ape)
library(phangorn)
library(wesanderson)

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

bound_p_trip <- function(M,ell,alpha,theta){
    #ell <- 1/(ng+1)
	du <- ell
    dv <- 2*ell
	if (alpha > 0){
		P <- p_trip_OU(M,dv,du,alpha,theta)
	} else {
		P <- p_trip_BM(M,dv,du,theta)
	}
	return(max(0,1-2*(1-P)))

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
    bound <- max(0, 1-2*choose(N,3)*(1-p0))
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
		coeff <- 2*n*choose(n/2^k,2) 
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


gg_m <- ggplot(dplyr::filter(df_m,alpha==0), aes(x=ng, y=log10(min_m), col=theta, group=theta))+
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


####################################
# PLOT MIN PTRIP
###################################
ng <- seq(2, 15, by=1)
ell <- seq(1/15, 1/3, by=0.001)
alpha <- c(0, 1,2)
theta <- c(0,0.5)
M <- 10^seq(0,4,by=0.1)

pars <- crossing(M=M,ell=ell, alpha=alpha,theta=theta)
res <- pars %>% pmap_dbl(bound_p_trip)

min_ptrip <- pars
min_ptrip$p <- res


gg_ptrip <- ggplot(min_ptrip, aes(x=ell, y=log10(M), fill=p, col=p))+
	facet_grid(theta~alpha, labeller=label_bquote(cols=alpha:.(alpha), rows=eta:.(theta)))+
	#facet_wrap(~alpha, nrow=1, labeller=label_bquote(alpha:.(alpha)))+
	geom_tile()+
	scale_fill_viridis(option="magma")+
	scale_color_viridis(option="magma")+
	theme(legend.position="bottom")+
	scale_x_continuous(expand=c(0,0))+
	scale_y_continuous(expand=c(0,0))+
	labs(y="# traits (log10)", x="Minimal branch length", fill="Triplet reconstruction guarantee",
		col="Triplet reconstruction guarantee")+
		guides(fill=guide_colourbar(
				frame.colour="black",
				barheight=unit(0.4, "cm")
				))



cowplot::plot_grid(
	gg_ptrip,
	gg_m,
	labels=c("A","B"), rel_widths=c(1,0.7))
