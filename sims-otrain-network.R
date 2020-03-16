rm(list=ls())
libraries_check <- c("data.table", "Rcpp", "RcppArmadillo","igraph","mvtnorm")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs)
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

sourceCpp("RandomizationInference.cpp")

# all possible simulation settings
sim_settings <- expand.grid(d=0,t=0,m=1/2,mean_neighs=(1:3)/4)
nreps <- 10
# initialize for parallel MC jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
  sim_settings <- sim_settings[rep(1:nrow(sim_settings),each=nreps),]
}
(seed <- as.integer(args[1]))

(delta <- sim_settings[seed,"d"]) # individual effect
(tau <- sim_settings[seed,"t"]) # social influence effect
(m <- sim_settings[seed,"m"]) ## number of treated
(mean_neighbors <- sim_settings[seed,"mean_neighs"]) # av. number of neighbors

# use observed covariate data from otrain study ===============================
otrain <- read.csv("data_processed_202001.csv")
otrain <- cbind("i"=1:nrow(otrain),otrain)
head(otrain)
sort(names(otrain))

# convert to data.table format
onedat <- data.table(otrain)
setkey(onedat)
# remove any missing data
onedat <- onedat[rowSums(is.na(onedat[, list(age,gender,priorrel)]))==0]
setkey(onedat)
(n <- nrow(onedat)) # sample size
onedat[, con := NULL] # remove observed treatment

# simulate adjacency matrix for a social network ==============================
## generate single realization of a linear preferential attachment network
(mean_neighbors <- round(mean_neighbors*n))
g <- sample_pa(n, power = 1, m = mean_neighbors/2, directed = FALSE)
a_mtx <- as.matrix(as_adjacency_matrix(g)); rm(g)
i_perm <- sample(n,n,replace=FALSE) ## randomly permute indices
a_mtx <- a_mtx[i_perm,i_perm]
all(diag(a_mtx) == 0L) # check diagonal entries are all zero
isSymmetric.matrix(a_mtx) # check if adjacency matrix is symmetric
a_mtx_rowsums <- rowSums(a_mtx) # number of neighbors
a_mtx_rowsums.nonzero <- sapply(rowSums(a_mtx), max, 1) # for use as denominator
a_mtx_edgelist <- as.matrix(get.data.frame(graph.adjacency(a_mtx)))
summary(rowSums(a_mtx))
onedat[, neighs := a_mtx_rowsums]
setkey(onedat)

# generate set of all possible randomized treatment assignments ===============
(m <- round(m*n))
mc_m <- 1e4 # generate approximation of sample space Omega
omega_mc <- sapply(1:mc_m, function(i) {
  x <- rep(0L, n)
  ## assume complete randomization with fixed number assigned to treatment
  x[sample.int(n, m, replace = FALSE)] <- 1L
  return(x)
})
all(colSums(omega_mc)==m)
dim(omega_mc)

# generate uniformity trial outcomes ==========================================
## generate new covariate for prior relationship
(p.priorrel <- mean(onedat$priorrel==0))
onedat[, priorrel := sapply(onedat$neighs, function(A) 
  rbinom(1,A,prob = p.priorrel))]
setkey(onedat)
all(onedat[, priorrel <= neighs])

## parameter estimates from fitted model with 'sense of control' as outcome
onedat[, y0 := -0.16*neighs + 0.01*age + 0.15*as.integer(gender=="Male") +
         0.12*priorrel]
## correlated uniformity outcomes between individuals
y0_corr <- rmvnorm(n=1,mean=onedat$y0,sigma=crossprod(a_mtx)/(n*n))[1,]
onedat[, y0 := y0_corr]; rm(y0_corr)
setkey(onedat)

# simulate observed data ======================================================
nsims <- 1e3 # number of simulated datasets
sims_list <- NULL
ptm=proc.time()[3]
for (id in 1:nsims) {
  sim_s <- sample(mc_m,1) # randomly sample a treatment assignment
  sim_X <- omega_mc[,sim_s]
  ## treated neighbors for each treatment assignment
  sim_Z <- (a_mtx %*% sim_X)[,1]/a_mtx_rowsums.nonzero
  ## check proportion of treated neighbors <= 1
  max(sim_Z)<=1
  
  simdat <- cbind(onedat[,list(i,age,gender,priorrel,neighs,y0)],
                  "X"=sim_X,"Z"=sim_Z)
  setkey(simdat)
  # generate observed outcomes under true causal model
  simdat[, Y := y0 + delta*X + tau*Z]
  setkey(simdat)
  
  # carry out procedure for testing 
  res.c <- OneCausalModel_edgelist_MC(
    obs_Ys = simdat$Y, obs_Xs = simdat$X,
    age = simdat$age, gender = as.integer(simdat$gender=="Male"),
    priorrel = simdat$priorrel,
    mc_Xs = omega_mc[,sample(x = mc_m,size = nsims, replace = FALSE)], 
    interfere_edges = a_mtx_edgelist,
    d_H0 = seq(from=-.02,to=.02,by=.001),
    t_H0 = seq(from=-.10,to=.10,by=.005),
    model = 0)
  
  res <- data.table(do.call(rbind, lapply(res.c, unlist)))
  resnames <- c("model", "delta", "tau", "pv.ssr_y0lm")
  setnames(res, resnames)
  
  # linear regression assuming independent observations
  fit_lm <- lm(Y~X+Z+neighs+age+gender+priorrel, data=simdat)
  fit_null <- lm(Y~neighs+age+gender+priorrel, data=simdat)
  
  sim_res <- list()
  sim_res[[1]] <- c("true_d"=delta,"true_t"=tau,
                    "m"=m,"mean_neighs"=mean_neighbors)
  sim_res[[2]] <- res
  sim_res[[3]] <- fit_lm
  sim_res[[4]] <- anova(fit_lm,fit_null) # F-test of both effects being zero
  
  sims_list[[id]] <- sim_res
  
  cat(id, "sims; time taken (mins)", round((proc.time()[3]-ptm)/60), "\n")
}
save(sims_list,file=paste0("sims-network-res-", seed,".Rdata"))
q()
