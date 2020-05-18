rm(list=ls())
libraries_check <- c("data.table","igraph")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs)
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

load(file="ICPSR_37070-R/DS0001/37070-0001-Data.rda") # raw data
dim(da37070.0001)

# convert to data.table format
Data <- data.table(da37070.0001)
setkey(Data)

# convert ID factors to integers
Data[, ID.f := ID]
Data[, ID := NULL]
Data[, ID := as.numeric(levels(ID.f))[ID.f]]
Data[sample(nrow(Data),10),list(ID.f,ID)] # random sample to check
Data[is.na(ID),list(ID.f,ID)]
Data[, ID.f := NULL]
## remove 999 or NA IDs
Data <- Data[!is.na(ID) & ID < 999]

# convert treatment factors to integers
Data[, TREAT.f := as.character(levels(TREAT))[TREAT]]
Data[sample(nrow(Data),10),list(TREAT.f,TREAT)] # random sample to check
Data[, TREAT := sapply(gsub(gsub(
  pattern="[(]",Data[, TREAT.f],replacement=""),
  pattern="[)]",Data[, TREAT.f],replacement=""),function(x) 
    as.integer(strsplit(x, " ")[[1]][1]))]
sapply(0:2, function(t) Data[TREAT==t,unique(TREAT.f)])
Data[, TREAT.f := NULL]
Data <- Data[!is.na(TREAT)]

Data[, SCHTREAT.f := as.character(levels(SCHTREAT))[SCHTREAT]]
Data[sample(nrow(Data),10),list(SCHTREAT.f,SCHTREAT)] # random sample to check
Data[, SCHTREAT := sapply(gsub(gsub(
  pattern="[(]",Data[, SCHTREAT.f],replacement=""),
  pattern="[)]",Data[, SCHTREAT.f],replacement=""),function(x) 
    as.integer(strsplit(x, " ")[[1]][1]))]
sapply(0:1, function(t) Data[SCHTREAT==t,unique(SCHTREAT.f)])
Data[, SCHTREAT.f := NULL]

# convert grade factors to integers
Data[, GRC.f := as.character(levels(GRC))[GRC]]
Data[sample(nrow(Data),10),list(GRC.f,GRC)] # random sample to check
Data[, GRC := sapply(strsplit(x=Data[, GRC.f],split=" "),function(x) 
  as.integer(strsplit(x[2],split="th")[[1]][1]))]
Data[sample(nrow(Data),10),list(GRC.f,GRC)] # random sample to check
Data[,any(is.na(GRC))]
Data[, GRC.f := NULL]

# summaries of the schools ====================================================
length(unique(Data$SCHID)) # number of unique schools
sum(Data[, all(SCHTREAT==1,na.rm=TRUE), by=SCHID][,V1])
# schools assigned to treatment and treated/control students within each school
Data[, list(unique(SCHRB),.N,sum(SCHTREAT),sum(TREAT==1),sum(TREAT==2)),
     by=SCHID]
summary(Data[,list(sum(TREAT>0),.N), by=SCHID])
Data[, as.list(5:8 %in% unique(GRC)), by=SCHID] # unique grade levels

# redefine treatment ==========================================================
setnames(Data,"TREAT","ELIGIBLE") # eligible to be in seed group
setkey(Data)
Data[, TREAT := (ELIGIBLE==1)*(SCHTREAT==1)] # individual treatment
Data[, list(.N,sum(SCHTREAT),sum(ELIGIBLE==1),sum(ELIGIBLE==2),sum(TREAT==1)),
     by=SCHID]

# treated schools only ========================================================
Data <- Data[SCHTREAT==1]
setkey(Data)


# construct adjacency matrix ==================================================
## adjancency matrix by school ================================================
sch_dat <- by(Data,Data$SCHID, function(x) data.table(x))
sch_a_mtx <- lapply(sch_dat, function(schdat) {
  # construct adjacency matrix using students' self-reported friends
  friends <- schdat[,c("ID",paste0("ST",1:10)),with=FALSE]
  friends <- data.table(apply(friends, 2, as.integer))
  a_mtx <- NULL
  for (j in 1:max(unlist(friends),na.rm=TRUE)) {
    i <- as.integer(friends[j, -1, with=FALSE])
    ## unique values after removing 999 or NA or 0 or own IDs
    i <- unique(i[!is.na(i) & (i > 0) & (i < 999) & (i != friends[j,ID])])
    if (length(i)>0) {
      # note direction of effect: alter (being reported) --> ego (reporting)
      a_mtx <- rbind(a_mtx,cbind("i"=friends[j,ID],"j"=i))
    }
  }
  a_mtx <- unique(a_mtx) # remove duplicates
  return(a_mtx)
})
any(unlist(lapply(sch_a_mtx, function(a_mtx) {
  # check if adjacencies are symmetric
  isSymmetric(as.matrix(as_adj(graph_from_edgelist(a_mtx))))
})))

# create unique IDs by combining school and within-school IDs
max(unlist(lapply(sch_dat, nrow))) # max number of students in a school
Data[, id := SCHID*1000L+ID]
# Adjacency matrix for all schools
all_a_mtx <- do.call(rbind,lapply(1:length(sch_a_mtx), function(s) {
  as.integer(names(sch_a_mtx)[s])*1000L+sch_a_mtx[[s]]
}))

# school with the densest adjacency matrix
sch_id <- names(sch_a_mtx)[
  which.max(unlist(lapply(1:length(sch_a_mtx), function(s) {
    nrow(sch_a_mtx[[s]])/nrow(sch_dat[[s]])
  })))]
onesch <- sch_dat[[sch_id]]
setkey(onesch)
onesch_a_mtx.el <- sch_a_mtx[[sch_id]]
onesch_a_mtx.GRC <- matrix(NA,nrow=nrow(onesch_a_mtx.el),ncol=2)
a_mtx.n <- max(onesch[,max(ID)],max(onesch_a_mtx.el))
onesch_a_mtx <- matrix(0,nrow=a_mtx.n,ncol=a_mtx.n)
for (si in 1:nrow(onesch_a_mtx.el)) {
  onesch_a_mtx[onesch_a_mtx.el[si,"i"],
               onesch_a_mtx.el[si,"j"]] <- 1L
  # across grade levels?
  onesch_a_mtx.GRC[si,] <- c(
    onesch[ID==onesch_a_mtx.el[si,"i"],GRC],
    onesch[ID==onesch_a_mtx.el[si,"j"],GRC])
}
all.equal(sum(onesch_a_mtx),nrow(sch_a_mtx[[sch_id]]))
sum(onesch_a_mtx.GRC[,1]<onesch_a_mtx.GRC[,2]) # friends with younger student
sum(onesch_a_mtx.GRC[,1]>onesch_a_mtx.GRC[,2]) # friends with older student
sum(onesch_a_mtx.GRC[,1]==onesch_a_mtx.GRC[,2]) # friends within same grade

pdf(paste0("friends-indivadj-SCHID_", sch_id, ".pdf"), width=9, height=9)
par(pty="s",mar=rep(6,4), bg=NA)
i_s <- 1:nrow(onesch_a_mtx)
a_mtx_dt <- as.matrix(get.data.frame(graph.adjacency(onesch_a_mtx[i_s,i_s])))
plot(range(i_s),rev(range(i_s)),
     type="n",xlab="j",ylab="i",cex.lab=3,
     axes = FALSE, ylim=rev(range(i_s)))
title(main=paste0("School ID ", sch_id),cex.main=3.5)
axis(1,range(i_s),range(i_s),col="grey80",lwd=2,cex.axis=2.5)
axis(2,range(i_s),range(i_s),col="grey80",lwd=2,cex.axis=2.5)
box("plot",col = "grey80")
points(y=a_mtx_dt[,"from"],x=a_mtx_dt[,"to"],pch=15, cex=.5)
rm(a_mtx_dt)
dev.off()

save(Data,all_a_mtx,file="friends-data.Rdata")
q()
