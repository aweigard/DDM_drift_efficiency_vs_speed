rm(list=ls())

### load data
dat<-read.csv("numerosity_trimmed_combined.csv")

### remove irrelevant rows
dat<-dat[dat$type!="",]

# make stim, response, RT, acc columns
dat$S<-dat$type

dat$R<-NA
dat[dat$key_press==".",]$R<-"MANY"
dat[dat$key_press==",",]$R<-"FEW"


dat$C<-NA
dat[dat$difficulty=="easy","C"]<-"easy"
dat[dat$difficulty=="hard","C"]<-"hard"


dat$RT<-round(as.numeric(dat$rt)/1000,3)
dat[dat$key_press==-1,c("R","RT")]<-NA

dat$s<-as.factor(dat$subject_id)
dat$S<-as.factor(dat$S)
dat$C<-as.factor(dat$C)
dat$R<-as.factor(dat$R)


dat$acc<-dat$S==tolower(dat$R)

# remove fast guesses
dat.c<-dat[dat$RT>=0.200 | is.na(dat$RT),]
# excludes <1% of the data, very high actually
length(dat.c$s)/length(dat$s)
#[1]  0.9958669


# summary data: 
sumdat<-data.frame(s=unique(dat.c$s))

sumdat$easy.crt<-NA
sumdat$hard.crt<-NA
sumdat$easy.ert<-NA
sumdat$hard.ert<-NA
sumdat$easy.acc<-NA
sumdat$hard.acc<-NA

sumdat$easy.rtvar<-NA
sumdat$hard.rtvar<-NA
sumdat$easy.n<-NA
sumdat$hard.n<-NA

sumdat$n<-NA
sumdat$ov.acc<-NA
sumdat$p.MANY<-NA
sumdat$omit.p<-NA
sumdat$fast.guess.p<-NA

for (n in unique(dat.c$s)){
  
  tmp<-dat.c[dat.c$s==n,]
  tmp.raw<-dat[dat$s==n,]
  
  sumdat[sumdat$s==n,]$easy.crt<-mean(tmp[tmp$C=="easy" & tmp$acc,]$RT,na.rm=TRUE)
  sumdat[sumdat$s==n,]$hard.crt<-mean(tmp[tmp$C=="hard" & tmp$acc,]$RT,na.rm=TRUE)
  sumdat[sumdat$s==n,]$easy.ert<-mean(tmp[tmp$C=="easy" & !tmp$acc,]$RT,na.rm=TRUE)
  sumdat[sumdat$s==n,]$hard.ert<-mean(tmp[tmp$C=="hard" & !tmp$acc,]$RT,na.rm=TRUE)
  sumdat[sumdat$s==n,]$easy.acc<-mean(tmp[tmp$C=="easy",]$acc,na.rm=TRUE)
  sumdat[sumdat$s==n,]$hard.acc<-mean(tmp[tmp$C=="hard",]$acc,na.rm=TRUE)
  
  sumdat[sumdat$s==n,]$easy.rtvar<-var(tmp[tmp$C=="easy" & tmp$acc,]$RT,na.rm=TRUE)
  sumdat[sumdat$s==n,]$hard.rtvar<-var(tmp[tmp$C=="hard" & tmp$acc,]$RT,na.rm=TRUE)
  sumdat[sumdat$s==n,]$easy.n<-length(tmp[tmp$C=="easy" & !is.na(tmp$acc),]$RT)
  sumdat[sumdat$s==n,]$hard.n<-length(tmp[tmp$C=="hard" & !is.na(tmp$acc),]$RT)
  
  sumdat[sumdat$s==n,]$n<-length(tmp$acc)
  sumdat[sumdat$s==n,]$ov.acc<-mean(tmp$acc,na.rm=TRUE)
  sumdat[sumdat$s==n,]$omit.p<-mean(is.na(tmp$acc),na.rm=TRUE)
  sumdat[sumdat$s==n,]$p.MANY<-mean(tmp[!is.na(tmp$R),]$R=="MANY")
  sumdat[sumdat$s==n,]$fast.guess.p<-(1-length(tmp$s)/length(tmp.raw$acc))
  
}

colnames(sumdat)<-c("s",paste0("num.",colnames(sumdat[,colnames(sumdat)!="s"])))
sumdat.numerosity<-sumdat

# alex updated to make consistent with ABCD/HCP
sumdat$inc<-(sumdat$num.ov.acc>0.55 & sumdat$num.fast.guess.p<0.10 & sumdat$num.n>=150 & sumdat$num.omit.p<0.25)
write.csv(sumdat[sumdat$inc,]$s,"num_inc_final.csv")


# final data frame following exclusions
dat.dmc<-dat.c[dat.c$subject_id%in%sumdat[sumdat$inc,"s"],c("s","S","C","R","RT","acc")]

# make s a factor with only existing people
dat.dmc$s<-as.factor(as.character(dat.dmc$s))

######################################
# broad prior fits, traditional DDM ##
######################################

####################################
source ("dmc/dmc.R")
load_model ("DDM","ddm_omit.R") 
# ddm_omit is a custom function, if not modeling omissions, use ddm.R

model <- model.dmc(
  p.map = list(a="1",v=c("S","C"),z="1",d="1",sz="1",
               sv="1",t0="1",st0="1",
               censor="1",gf="1"), 
  responses = c("FEW","MANY"),
  factors=list(S=c("few","many"),C=c("easy","hard")),
  match.map = list(M=list(many.easy="MANY",many.hard="MANY",
                          few.easy="FEW",few.hard="FEW")),
  constants = c(sv=0,sz=0,d=0,censor=3.00), 
  type="rd")


p1 <- c(a=1,v.many.easy=3, v.few.easy = 3, v.many.hard= 3, v.few.hard= 3, 
        z=0.5,t0=0.3, st0=0.1,gf=0)

p.prior <- prior.p.dmc(
  dists = rep("tnorm",9),
  p1=p1,                           
  p2=c(.5,1,1,1,1,.1,.1,.05,1),
  lower=c(0,NA,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

mdi <- data.model.dmc(dat.dmc, model)

samples <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = mdi , thin = 5)
save(samples ,file="num_ddm_broad.RData")

# run on greatlakes

source("greatlakes.R")
run.greatlakes.dmc("num_ddm_broad","DDM","ddm_omit.R","asweigar",
                   "asweigar0",40,wall.hours=10)

######################################
# LBA ################################
######################################

#### Design and contaminant (gf) ----
load_model ("LBA","lba_BpvGF.R")

model <- model.dmc(
  p.map = list(A="1",B="R",
               mean_v=c("S","C","M"),
               sd_v="M",t0="1",censor="1",gf="1"), 
  responses = c("FEW","MANY"),
  match.map = list(M=list(many.easy="MANY",many.hard="MANY",
                          few.easy="FEW",few.hard="FEW")),
  factors=list(S=c("few","many"),C=c("easy","hard")),
  constants = c(sd_v.true=1,censor=3.00), 
  type="norm")



p1 <-     c(A=1,B.FEW=1,B.MANY=1,
            mean_v.few.easy.true=2,mean_v.many.easy.true=2,mean_v.few.hard.true=2, mean_v.many.hard.true=2,
            mean_v.few.easy.false=1,mean_v.many.easy.false=1,mean_v.few.hard.false=1, mean_v.many.hard.false=1,
            sd_v.false=1,
            t0=0.3,gf=0)

p.prior <- prior.p.dmc(
  dists = rep("tnorm",14),
  p1=p1,                           
  p2=c(rep(1,12),.1,1),
  lower=c(rep(0,3),rep(NA,8),0,.1,NA),upper=c(rep(NA,12),2,NA)
)


mdi <- data.model.dmc(dat.dmc, model)

samples <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = mdi , thin = 5)
save(samples ,file="num_lba_broad.RData")

# run on greatlakes

source("greatlakes.R")
run.greatlakes.dmc("num_lba_broad","LBA","lba_BpvGF.R","asweigar",
                   "asweigar0",40,wall.hours=10)


