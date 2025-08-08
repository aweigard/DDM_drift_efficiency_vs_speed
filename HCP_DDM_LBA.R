
### load data
rm(list=ls())

# load data
zero.dat<-read.csv("HCP_0bk_trials_full.csv")
two.dat<-read.csv("HCP_2bk_trials_full.csv")

# exclude session 2 data
zero.dat<-zero.dat[zero.dat$session==1,]
two.dat<-two.dat[two.dat$session==1,]


# make and format DMC columns

zero.dat$s<-zero.dat$subject_num
zero.dat$S<-zero.dat$targettype
zero.dat$R<-zero.dat$resp

zero.dat[zero.dat$targettype=="target","S"]<-"tar"
zero.dat[zero.dat$targettype=="lure","S"]<-"lur"
zero.dat[zero.dat$targettype=="nonlure","S"]<-"non"

zero.dat[zero.dat$resp==2 & !is.na(zero.dat$resp),"R"]<-"YES"
zero.dat[zero.dat$resp==3 & !is.na(zero.dat$resp),"R"]<-"NO"

zero.dat[is.na(zero.dat$resp),"RT"]<-NA

two.dat$s<-two.dat$subject_num
two.dat$S<-two.dat$targettype
two.dat$R<-two.dat$resp

two.dat[two.dat$targettype=="target","S"]<-"tar"
two.dat[two.dat$targettype=="lure","S"]<-"lur"
two.dat[two.dat$targettype=="nonlure","S"]<-"non"

two.dat[two.dat$resp==2 & !is.na(two.dat$resp),"R"]<-"YES"
two.dat[two.dat$resp==3 & !is.na(two.dat$resp),"R"]<-"NO"

two.dat[is.na(two.dat$resp),"RT"]<-NA

# accuracy
zero.dat$acc<-NA
zero.dat[!is.na(zero.dat$R),]$acc<-zero.dat[!is.na(zero.dat$R),]$accuracy

two.dat$acc<-NA
two.dat[!is.na(two.dat$R),]$acc<-two.dat[!is.na(two.dat$R),]$accuracy


#exclude some responses on an extra key
two.dat<-two.dat[two.dat$resp!=4 | is.na(two.dat$resp),]

# summary accuracy stats

sum.base<-data.frame(unique(c(zero.dat$s,two.dat$s)));colnames(sum.base)<-"s"

sum.base$mrt.0.tar<-NA
sum.base$mrt.0.lur<-NA
sum.base$mrt.0.non<-NA
sum.base$sdrt.0.tar<-NA
sum.base$sdrt.0.lur<-NA
sum.base$sdrt.0.non<-NA
sum.base$acc.0.tar<-NA
sum.base$acc.0.lur<-NA
sum.base$acc.0.non<-NA

sum.base$acc.0.total<-NA
sum.base$p.omit.0<-NA

sum.base$mrt.2.tar<-NA
sum.base$mrt.2.lur<-NA
sum.base$mrt.2.non<-NA
sum.base$sdrt.2.tar<-NA
sum.base$sdrt.2.lur<-NA
sum.base$sdrt.2.non<-NA
sum.base$acc.2.tar<-NA
sum.base$acc.2.lur<-NA
sum.base$acc.2.non<-NA

sum.base$acc.2.total<-NA
sum.base$p.omit.2<-NA

for (r in 1:length(sum.base$s)){
  tmp0<-zero.dat[zero.dat$s==sum.base$s[r],]
  tmp2<-two.dat[two.dat$s==sum.base$s[r],]
  
  sum.base$mrt.0.tar[r]<-mean(tmp0[tmp0$S=="tar",]$RT,na.rm=TRUE)
  sum.base$mrt.0.lur[r]<-mean(tmp0[tmp0$S=="lur",]$RT,na.rm=TRUE)
  sum.base$mrt.0.non[r]<-mean(tmp0[tmp0$S=="non",]$RT,na.rm=TRUE)
  sum.base$sdrt.0.tar[r]<-sd(tmp0[tmp0$S=="tar",]$RT,na.rm=TRUE)
  sum.base$sdrt.0.lur[r]<-sd(tmp0[tmp0$S=="lur",]$RT,na.rm=TRUE)
  sum.base$sdrt.0.non[r]<-sd(tmp0[tmp0$S=="non",]$RT,na.rm=TRUE)
  sum.base$acc.0.tar[r]<-mean(tmp0[tmp0$S=="tar",]$acc,na.rm=TRUE)
  sum.base$acc.0.lur[r]<-mean(tmp0[tmp0$S=="lur",]$acc,na.rm=TRUE)
  sum.base$acc.0.non[r]<-mean(tmp0[tmp0$S=="non",]$acc,na.rm=TRUE)
  
  sum.base$acc.0.total[r]<-mean(tmp0$acc,na.rm=TRUE)
  sum.base$p.omit.0[r]<-mean(is.na(tmp0$acc))
  
  sum.base$mrt.2.tar[r]<-mean(tmp2[tmp2$S=="tar",]$RT,na.rm=TRUE)
  sum.base$mrt.2.lur[r]<-mean(tmp2[tmp2$S=="lur",]$RT,na.rm=TRUE)
  sum.base$mrt.2.non[r]<-mean(tmp2[tmp2$S=="non",]$RT,na.rm=TRUE)
  sum.base$sdrt.2.tar[r]<-sd(tmp2[tmp2$S=="tar",]$RT,na.rm=TRUE)
  sum.base$sdrt.2.lur[r]<-sd(tmp2[tmp2$S=="lur",]$RT,na.rm=TRUE)
  sum.base$sdrt.2.non[r]<-sd(tmp2[tmp2$S=="non",]$RT,na.rm=TRUE)
  sum.base$acc.2.tar[r]<-mean(tmp2[tmp2$S=="tar",]$acc,na.rm=TRUE)
  sum.base$acc.2.lur[r]<-mean(tmp2[tmp2$S=="lur",]$acc,na.rm=TRUE)
  sum.base$acc.2.non[r]<-mean(tmp2[tmp2$S=="non",]$acc,na.rm=TRUE)
  
  sum.base$acc.2.total[r]<-mean(tmp2$acc,na.rm=TRUE)
  sum.base$p.omit.2[r]<-mean(is.na(tmp2$acc))
  
}

inc.0b<-sum.base[sum.base$acc.0.total>=.55 & sum.base$p.omit.0<=.25,"s"]
inc.2b<-sum.base[sum.base$acc.2.total>=.55 & sum.base$p.omit.2<=.25,"s"]
exc.2b<-sum.base[sum.base$acc.2.total<.55 | sum.base$p.omit.2>.25,"s"]

save(sum.base,file="sum_HCP.RData")
load("sum_HCP.RData")

# inclusion for each task

inc.0b<-sum.base[sum.base$acc.0.total>=.55 & sum.base$p.omit.0<=.25,"s"]
inc.2b<-sum.base[sum.base$acc.2.total>=.55 & sum.base$p.omit.2<=.25,"s"]

# final groups

base.0.dmc<-zero.dat[zero.dat$s%in%inc.0b,]
base.0.dmc<-base.0.dmc[,c("s","S","R","RT")]
base.0.dmc$s<-factor(base.0.dmc$s)

base.2.dmc<-two.dat[two.dat$s%in%inc.2b,]
base.2.dmc<-base.2.dmc[,c("s","S","R","RT")]
base.2.dmc$s<-factor(base.2.dmc$s)

# exclude fast guesses

base.0.dmc<-base.0.dmc[base.0.dmc$RT>=0.200 | is.na(base.0.dmc$RT),]
base.2.dmc<-base.2.dmc[base.2.dmc$RT>=0.200 | is.na(base.2.dmc$RT),]

######################################
# broad prior fits, traditional DDM ##
######################################

####################################
source ("dmc/dmc.R")
load_model ("DDM","ddm_omit.R") 
# ddm_omit is a custom function, if not modeling omissions, use ddm.R

model <- model.dmc(
  p.map = list(a="1",v=c("S"),z="1",d="1",sz="1",
               sv="1",t0="1",st0="1",
               censor="1",gf="1"), 
  responses = c("NO","YES"),
  match.map = list(M=list(tar="YES",lur="NO",non="NO")),
  factors=list(S=c("lur","non","tar")),
  constants = c(sv=0,sz=0,d=0,censor=2.00), 
  type="rd")

p1 <- c(a=1,v.lur=3, v.non = 3, v.tar= 3,
        z=0.5,t0=0.3, st0=0.1,gf=0)

p.prior <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1,                           
  p2=c(.5,1,1,1,.1,.1,.05,1),
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)

zero.mdi <- data.model.dmc(base.0.dmc, model)
two.mdi <- data.model.dmc(base.2.dmc, model)

sZero  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = zero.mdi, thin = 5)
save(sZero,file="HCP_zero_ddm.RData")

sTwo  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = two.mdi, thin = 5)
save(sTwo,file="HCP_two_ddm.RData")


# run on greatlakes

source("greatlakes.R")
run.greatlakes.dmc("HCP_zero_ddm","DDM","ddm_omit.R","asweigar",
                   "asweigar0",40,wall.hours=10)


source("greatlakes.R")
run.greatlakes.dmc("HCP_two_ddm","DDM","ddm_omit.R","asweigar",
                   "asweigar0",40,wall.hours=10)


######################################
# LBA ################################
######################################

#### Design and contaminant (gf) ----
load_model ("LBA","lba_BpvGF.R")

model <- model.dmc(
  p.map = list(A="1",B="1",
               mean_v=c("S","M"),
               sd_v="M",t0="1",censor="1",gf="1"), 
  responses = c("NO","YES"),
  match.map = list(M=list(tar="YES",lur="NO",non="NO")),
  factors=list(S=c("lur","non","tar")),
  constants = c(sd_v.true=1,censor=2.00), 
  type="norm")

p1 <-     c(A=1,B=1,
            mean_v.tar.true=2,mean_v.lur.true=2,mean_v.non.true=2,
            mean_v.tar.false=1,mean_v.lur.false=1,mean_v.non.false=1,
            sd_v.false=1,
            t0=0.3,gf=0)

p.prior <- prior.p.dmc(
  dists = rep("tnorm",11),
  p1=p1,                           
  p2=c(rep(1,9),.1,1),
  lower=c(rep(0,2),rep(NA,6),0,.1,NA),upper=c(rep(NA,9),2,NA)
)
# par(mfcol=c(2,6)); for (i in names(p.prior)) plot.prior(i,p.prior)


zero.mdi <- data.model.dmc(base.0.dmc, model)
two.mdi <- data.model.dmc(base.2.dmc, model)

sZero  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = zero.mdi, thin = 5)
save(sZero,file="HCP_zero_lba.RData")

sTwo  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = two.mdi, thin = 5)
save(sTwo,file="HCP_two_lba.RData")


# run on greatlakes

source("greatlakes.R")
run.greatlakes.dmc("HCP_zero_lba","LBA","lba_BpvGF.R","asweigar",
                   "asweigar0",40,wall.hours=10)

source("greatlakes.R")
run.greatlakes.dmc("HCP_two_lba","LBA","lba_BpvGF.R","asweigar",
                   "asweigar0",40,wall.hours=10)



