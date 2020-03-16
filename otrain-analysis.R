rm(list=ls())
libraries_check <- c("data.table", "Rcpp", "RcppArmadillo","igraph",
                     "ggplot2","directlabels","gbm","xtable")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs)
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

sourceCpp("RandomizationInference.cpp")

otrain <- read.csv("data_processed_202001.csv")
otrain <- cbind("id"=1:nrow(otrain),otrain)
head(otrain)
sort(names(otrain))

# convert to data.table format
onedat <- data.table(otrain)
setkey(onedat)
# remove any missing data
onedat <- onedat[rowSums(is.na(onedat[, list(age,gender,priorrel)]))==0]
## check any observations with missing outcomes
onedat[,all(!is.na(control))]
(n <- nrow(onedat)) # sample size

# summaries of baseline covariates ============================================
by(onedat,INDICES=list(onedat$con,onedat$role),function(x) nrow(x))
by(onedat$age,INDICES=list(onedat$con,onedat$role),function(x) 
  c(mean(x),sd(x)))
by(onedat$gender,INDICES=list(onedat$con,onedat$role),function(x) 
  c(mean(x=="Female"),sd(x=="Male")))
by(onedat$priorrel,INDICES=list(onedat$con,onedat$role),function(x) 
  c(mean(x),sd(x)))

# construct adjacency matrix using classroom, group, and roles (source/target)
a_mtx <- matrix(NA,nrow=n,ncol=n)
for (i in 1:n) {
  alpha_i <- onedat[i,list(sample,group,role)]
  for (j in 1:n) {
    if (j != i) {
      alpha_j <- onedat[j,list(sample,group,role)]
      # check if same classroom and same group
      same <- (alpha_i$sample==alpha_j$sample) & (alpha_i$group==alpha_j$group)
      # only targets can be connected to sources within the same group
      ## one direction of possible interference: target <-- source
      same <- same & (alpha_i$role=="T" & alpha_j$role=="S")
      a_mtx[i,j] <- as.integer(same)
    }
  }
}
diag(a_mtx) <- 0L # set diagonal entries to zero
a_mtx_rowsums <- rowSums(a_mtx) # number of neighbors
a_mtx_rowsums.nonzero <- sapply(rowSums(a_mtx), max, 1) # for use as denominator
a_mtx_edgelist <- as.matrix(get.data.frame(graph.adjacency(a_mtx)))
table(rowSums(a_mtx))
table(rowSums(a_mtx),onedat$role) # check that sources have no neighbors
onedat[, neighs := a_mtx_rowsums]
setkey(onedat)

# generate set of all possible randomized treatment assignments ===============
## only sources can be assigned to treatment (exclusion)
onedat[, "X" := con*(role=="S")]
(m <- sum(onedat$X)) ## number of sources who will carry out exclusion
table(onedat[,sample],onedat[,X]) ## exclusion and inclusion conditions
(sample_id <- by(onedat,onedat[,sample],function(x) 
  length(unique(x[X==1,group])))) ## groups assigned to exclusion in each class
onedat[, X := NULL]
setkey(onedat)
mc_m <- 2e3 # generate approximation of sample space Omega
omega_mc <- sapply(1:mc_m, function(i) {
  x <- rep(0L, n)
  for (id in 1:length(sample_id)) {
    # unique groups within each classroom
    n.sample_id <- onedat[,sample]==names(sample_id)[id]
    group_id <- unique(onedat[n.sample_id,group])
    # randomly assign groups to exclusion condition within each class
    m.group_id <- sort(sample(group_id,sample_id[[id]],replace=FALSE))
    # indices of sources in the exclusion groups within each class
    m_id <- which( n.sample_id & (onedat[,group] %in% m.group_id) & 
                     (onedat[,role] == "S") )
    x[m_id] <- 1L # assign sources to exclusion
  }
  return(x)
})
dim(omega_mc)
## check that same number of groups assigned to treatment within each sample
unique(t(apply(omega_mc, 2, function(x) onedat[x==1,list(sample,group)][
  ,lapply(.SD, function(x1) length(unique(x1))), by=sample][,group])))
as.matrix(sample_id)[,1]

# observed outcomes ===========================================================
onedat[, Y := control]
setkey(onedat)

# carry out procedure for testing 
res.c <- OneCausalModel_edgelist_MC(
  obs_Ys = onedat$Y, obs_Xs = onedat$X,
  age = onedat$age, gender = as.integer(onedat$gender=="Male"),
  priorrel = onedat$priorrel,
  mc_Xs = omega_mc, 
  interfere_edges = a_mtx_edgelist,
  d_H0 = seq(from=0,to=1.5,by=.01), 
  t_H0 = seq(from=-1.5,to=0,by=.01),
  model = 0)

res <- data.table(do.call(rbind, lapply(res.c, unlist)))
resnames <- c("model", "delta", "tau", "pv.ssr_y0lm")
setnames(res, resnames)

save(res,file="otrain-all-pvalues-control.Rdata")

# results =====================================================================
plotfile <- "otrain-all-pvalues-control"
load(file=paste0(plotfile, ".Rdata"))
res[delta==0 & tau==0]
res[which.max(pv.ssr_y0lm)]
## CIs for causal effects
res[pv.ssr_y0lm >= 0.05, range(delta)] # 95% individual effect
res[pv.ssr_y0lm >= 0.05, range(tau)] # 95% social influence effect
res[pv.ssr_y0lm >= 0.01, range(delta)] # 99% individual effect
res[pv.ssr_y0lm >= 0.01, range(tau)] # 99% social influence effect
## assuming no inteference
res[pv.ssr_y0lm>=0.05 & tau==0, range(delta)] # 95% CI

results_sims <- res

pvs <- results_sims[,list(delta,tau,pv.ssr_y0lm)]
setnames(pvs,3,"pv"); setkey(pvs)

plot_halfwidth <- max(diff(range(pvs[,delta])),diff(range(pvs[,tau])))/2
onep <- ggplot(pvs, aes(x=delta,y=tau))+
  geom_vline(xintercept=0,colour="grey50",size=.2)+
  geom_hline(yintercept=0,colour="grey50",size=.2)+
  geom_point(size=.25,stroke=0, shape=20, color="grey50")+
  scale_x_continuous(name=expression(paste("Individual (",delta,")")),
                     limits=sum(range(pvs[,delta]))/2+c(-1,1)*plot_halfwidth)+
  scale_y_continuous(name=expression(paste("Social influence (",tau,")")),
                     limits=sum(range(pvs[,tau]))/2+c(-1,1)*plot_halfwidth)+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor=element_blank())

if(pvs[,any(pv>=.05)]) {
  onepl <- onep+
    geom_contour(aes(z=pv, color = ..level..),
                 breaks=c(.05,c(1,2,4,6,8)/10),size=.2)+
    scale_colour_gradient(na.value="white", low="blue", high="black", 
                          guide=F)
  onepl <- direct.label(onepl,list("top.points", cex=.5, fontfamily="sans"))
} else {
  onepl <- onep+geom_point(colour="white",size=0)
}
allplots <- list()
allplots[[1]] <- onepl +
  ggtitle(paste("95% confidence set")) +
  geom_point(data=pvs[which.max(pv)],aes(x=delta,y=tau),
             colour="black",size=.2,shape=4)

pdf(file=paste0(plotfile,".pdf"),width=3,height=3)
grid.arrange(grobs=allplots,ncol=1)
dev.off()
