rm(list=ls())
libraries_check <- c("data.table", "Rcpp", "RcppArmadillo","igraph")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs)
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# analyze treatment effects ===================================================
load(file="friends-data.Rdata")

# convert covariate factors to integers
## AGE
class(Data[,AGEC_NEW])
## GENDER
Data[, GENC.f := as.character(levels(GENC))[GENC]]
Data[sample(nrow(Data),10),list(GENC.f,GENC)] # random sample to check
Data[, GENC := sapply(gsub(gsub(
  pattern="[(]",Data[, GENC.f],replacement=""),
  pattern="[)]",Data[, GENC.f],replacement=""),function(x) 
    as.integer(strsplit(x, " ")[[1]][1]))]
sapply(0:1, function(t) Data[GENC==t,unique(GENC.f)])
Data[, GENC.f := NULL]
# composite outcome ===========================================================
out.subset <- paste0("PN",3:6)
out.composite <- c(out.subset,paste0(out.subset,"W2"))
all(out.composite %in% names(Data)) # check all outcome variables in data
out.numeric <- NULL
for (yy in 1:length(out.composite)) {
  yy.factor <- unlist(Data[, out.composite[yy], with=FALSE])
  yy.factor <- as.character(levels(yy.factor))[yy.factor]
  yy.numeric <- sapply(gsub(gsub(
    pattern="[(]",yy.factor,replacement=""),
    pattern="[)]",yy.factor,replacement=""),function(x) 
      as.numeric(strsplit(x, " ")[[1]][1]))
  # random sample to check
  print(cbind(Data[, out.composite[yy], with=FALSE],
              yy.numeric)[sample(nrow(Data),10),])
  names(yy.numeric) <- NULL
  # reverse-code based on signs of factor loadings in Table S3 of supp. materials
  yy.numeric <- 5L - yy.numeric
  
  out.numeric <- cbind(out.numeric,yy.numeric)
  rm(yy.factor,yy.numeric)
}
colnames(out.numeric) <- out.composite

Y <- rowMeans(out.numeric[,paste0(out.subset,"W2")],na.rm=TRUE)
Y0 <- rowMeans(out.numeric[,out.subset],na.rm=TRUE)
Data[, "Y" := Y]
Data[, "Y0" := Y0]
rm(Y, Y0)

# drop missing IDs and missing data
i_tokeep <- Data[!is.na(AGEC_NEW) & !is.na(GENC) & !is.na(Y0) & !is.na(Y),
                 id]
i_tokeep <- sort(i_tokeep)
Data <- Data[id %in% i_tokeep,]
nrow(Data)
Y <- Data[,Y]
summary(Y)
a_mtx_el <- all_a_mtx[rowSums(apply(all_a_mtx, 2, function(x) 
  x %in% i_tokeep))==ncol(all_a_mtx),]
## relabel IDs in adjacency matrix
a_mtx_edgelist <- matrix(NA,nrow=nrow(a_mtx_el),ncol=ncol(a_mtx_el))
new_ids <- cbind("id"=Data[,id],"id.new"=1:nrow(Data))
for (si in 1:nrow(a_mtx_el)) {
  a_mtx_edgelist[si,] <- c(
    "i"=as.integer(new_ids[new_ids[,"id"]==a_mtx_el[si,"i"],"id.new"]),
    "j"=as.integer(new_ids[new_ids[,"id"]==a_mtx_el[si,"j"],"id.new"])
  )
  if ((si %% 10000)==0) cat(si, "out of", nrow(a_mtx_el), "\n")
}
max(a_mtx_edgelist) <= length(Y)

# generate set of all possible randomized treatment assignments ===============
X <- Data[,TREAT]
table(X)
(n <- length(X))
## treated schools in each block
block.sch <- Data[,list("treated"=sum(TREAT)>0),by=list(SCHRB,SCHID)] 
setkey(block.sch)
(m.sch <- unique(block.sch[,list("treated"=sum(treated)),by=SCHRB][,treated]))
mc_m <- 1e2 # generate approximation of sample space Omega
omega_mc <- sapply(1:mc_m, function(i) {
  x <- rep(0L, n)
  ## schools selected to receive treatment
  x.sch <- sort(as.integer(unlist(
    by(data=block.sch$SCHID, INDICES=block.sch$SCHRB, sample,
       size=m.sch, replace=FALSE))))
  for (s in 1:length(x.sch)) {
    si <- Data[,SCHID]==x.sch[s] # students in treated school
    seed_eligible <- Data[si,ELIGIBLE>0] # seed-eligible students
    m_s <- Data[si,sum(TREAT)] # number treated in observed sample
    if (m_s==0) {
      m_s <- round(sum(seed_eligible)/2) # assign half of eligibles to treatment
    }
    x_s <- rep(0,sum(si))
    x_s[seed_eligible][sample(sum(seed_eligible),m_s,replace=FALSE)] <- 1L
    x[si] <- x_s
  }
  return(x)
})
dim(omega_mc)
## check that same number of schools assigned to treatment within each sample
m.sch==unique(apply(omega_mc, 2, function(x) {
  x.mc <- cbind(Data[,list(SCHRB,SCHID)],"TREAT"=as.integer(x))
  x.mc.block.sch <- x.mc[,list("treated"=sum(TREAT)>0),by=list(SCHRB,SCHID)]
  setkey(x.mc.block.sch)
  unique(x.mc.block.sch[,list("treated"=sum(treated)),by=SCHRB][,treated])
}))
## check that (almost) half of eligible students received treatment
summary(apply(apply(omega_mc, 2, function(x) {
  x.mc <- cbind(Data[,list(SCHID,ELIGIBLE)],"TREAT"=as.integer(x))
  x.mc[,sum(TREAT==1 & ELIGIBLE>0)/sum(ELIGIBLE>0),by=SCHID][,V1]
}),1,function(x) sort(unique(x))))
## check that ineligible students were not assigned to treatment
all(apply(omega_mc, 2, function(x) sum(Data[,ELIGIBLE==0]*x))==0)

summary(colSums(omega_mc))

d_H0.vals <- seq(from=-0.20, to=0.15, by=0.005); length(d_H0.vals)
t_H0.vals <- seq(from=-0.02, to=0.65, by=0.001); length(t_H0.vals)
# d_H0.vals <- t_H0.vals <- 0 # for testing sharp null

# initialize for parallel cluster jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(args[1]))
rm(args)

sourceCpp("RandomizationInference.cpp")

# carry out procedure for testing 
ptm=proc.time()[3]
res.c <- OneCausalModel_edgelist_MC(
  obs_Ys = Y, obs_Xs = X,
  age = Data$AGEC_NEW, gender = Data$GENC, priorrel = Data$Y0,
  mc_Xs = omega_mc, 
  interfere_edges = a_mtx_edgelist,
  d_H0 = d_H0.vals, t_H0 = t_H0.vals,
  model = 0)
proc.time()[3]-ptm

res <- data.table(do.call(rbind, lapply(res.c, unlist)))
resnames <- c("model", "delta", "tau", "pv.ssr_y0lm")
setnames(res, resnames)

save(res,file=paste0("friends-pvalues-",seed,".Rdata"))
q()
