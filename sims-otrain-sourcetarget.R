rm(list=ls())
libraries_check <- c("data.table", "Rcpp", "RcppArmadillo","igraph",
                     "lme4","lmerTest")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs)
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

sourceCpp("RandomizationInference.cpp")

# all possible simulation settings
sim_settings <- expand.grid(d=c(0,6,8,10)/10,t=c(0,-4,-8,-12)/10)
nreps <- 10
# initialize for parallel MC jobs
args <- 4
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  sim_settings <- sim_settings[rep(1:nrow(sim_settings),each=nreps),]
}
(seed <- as.integer(args[1]))

(delta <- sim_settings[seed,"d"]) # individual effect
(tau <- sim_settings[seed,"t"]) # social influence effect

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
onedat[, con := NULL] # remove observed treatment
(sample_id <- by(onedat,onedat[,sample],function(x) 
  length(unique(x[X==1,group])))) ## groups assigned to exclusion in each class
onedat[, X := NULL]
setkey(onedat)
mc_m <- 1e4 # generate approximation of sample space Omega
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

# generate uniformity trial outcomes ==========================================
## parameter estimates from fitted model with 'sense of control' as outcome
onedat[, y0 := -0.16*neighs + 0.01*age + 0.15*as.integer(gender=="Male") + 
         0.12*priorrel + rnorm(n,sd = 0.85)]
setkey(onedat)

# simulate observed data ======================================================
nsims <- 1e3 # number of simulated datasets; 2 mins per sim
sims_list <- NULL
ptm=proc.time()[3]
for (id in 1:nsims) {
  sim_s <- sample(mc_m,1) # randomly sample a treatment assignment
  sim_X <- omega_mc[,sim_s]
  ## neighbors assigned to exclusion for each randomized treatment
  sim_Z <- (a_mtx %*% sim_X)[,1]/a_mtx_rowsums.nonzero
  ## check proportion of treated neighbors <= 1
  max(sim_Z)<=1
  ## only targets have neighbors assigned to exclusion (Z>0);
  ## i.e., all sources must have Z=0
  all(sim_Z[onedat$role=="S"]==0)
  
  simdat <- cbind(onedat[,list(i,sample,group,age,gender,priorrel,neighs,y0)],
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
    d_H0 = seq(delta-1,delta+1,by = .05), 
    t_H0 = seq(tau-1.2,tau+1.2,by = .05),
    model = 0)
  
  res <- data.table(do.call(rbind, lapply(res.c, unlist)))
  resnames <- c("model", "delta", "tau", "pv.ssr_y0lm")
  setnames(res, resnames)
  
  # multilevel model with random intercept for each group
  simdat[, group.label := paste(sample,group,sep = ".")]
  fit_lmer <- lmerTest::lmer(
    Y~X+Z+neighs+age+gender+priorrel+(1|group.label)+(1|sample), 
    data=simdat)
  
  sim_res <- list()
  sim_res[[1]] <- c("delta"=delta,"tau"=tau)
  sim_res[[2]] <- res
  sim_res[[3]] <- fit_lmer
  
  sims_list[[id]] <- sim_res
  
  cat(id, "sims; time taken (mins)", round((proc.time()[3]-ptm)/60), "\n")
}
save(sims_list,file=paste0("sims-sourcetarget-res-", seed,".Rdata"))
q()
