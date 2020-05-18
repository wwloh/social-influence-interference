rm(list=ls())
libraries_check <- c("data.table","igraph",
                     "ggplot2","directlabels","gbm","xtable")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs)
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# results =====================================================================
subfolder <- "res/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
res_list <- list()
for (ll in myfiles) {
  # results
  load(paste0(subfolder,ll))
  res_list <- c(res_list, list(res))
  rm(res)
  cat(ll,"\n")
}
res <- rbindlist(res_list)
setkey(res)

# average over all MC p-values
res <- res[, colMeans(.SD), by=list(model,delta,tau)]
setnames(res,"V1","pv.ssr_y0lm")
setkey(res)

plotfile <- "friends-pvalues"
res[delta==0 & tau==0]
res[which.max(pv.ssr_y0lm)]
## CIs for causal effects
res[pv.ssr_y0lm >= 0.05, range(delta)] # 95% individual effect
res[pv.ssr_y0lm >= 0.05, range(tau)] # 95% social influence effect

results_sims <- res

pvs <- results_sims[,list(delta,tau,pv.ssr_y0lm)]
setnames(pvs,3,"pv"); setkey(pvs)

plot_halfwidth <- max(diff(range(pvs[,delta])),diff(range(pvs[,tau])))/2
onep <- ggplot(pvs, aes(x=delta,y=tau))+
  geom_vline(xintercept=0,colour="grey50",size=.2)+
  geom_hline(yintercept=0,colour="grey50",size=.2)+
  geom_point(size=.1,stroke=0, shape=20, color="grey50")+
  scale_x_continuous(name=expression(paste("Individual (",delta,")")),
                     limits=sum(range(pvs[,delta]))/2+c(-1,1)*plot_halfwidth)+
  scale_y_continuous(name=expression(paste("Social influence (",tau,")")),
                     limits=sum(range(pvs[,tau]))/2+c(-1,1)*plot_halfwidth)+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor=element_blank())

if(pvs[,any(pv>=.05)]) {
  onepl <- onep+
    geom_contour(aes(z=pv, color = ..level..),
                 breaks=c(.05,c(15,25,45,65,85)/100),size=.2)+
    scale_colour_gradient(na.value="white", low="blue", high="black", 
                          guide=F)
  onepl <- direct.label(onepl,list("top.points", cex=.4, fontfamily="sans"))
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
