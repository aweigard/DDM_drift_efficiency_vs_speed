rm(list=ls())
source ("dmc/dmc.R")

# This file contains code for fitting the linear ballistic accumulator (LBA)
# to ABCD n-back data. Code for fitting the diffusion decision model (DDM) to the same 
# data are available in a repository from one of our past projects: 
# https://github.com/aweigard/EEA_network_adaptation

################################################################
####### format n-back data and apply exclusion criteria ########
################################################################

# load trial-level n-back data, merged across all 
# baseline participants (once downloaded from aws)
abcd.nback<-read.csv("nback_revised.csv")

# make DMC-based columns

abcd.nback$S<-as.character(abcd.nback$enback_targettype)
abcd.nback[abcd.nback$S=="target",]$S<-"tar"
abcd.nback[abcd.nback$S=="lure",]$S<-"lur"
abcd.nback[abcd.nback$S=="nonlure",]$S<-"non"
abcd.nback$S<-factor(abcd.nback$S)

abcd.nback$R<-NA
abcd.nback[abcd.nback$S%in%c("tar") & abcd.nback$enback_stim_acc==1 & 
             !is.na(abcd.nback$enback_stim_resp),]$R<-"YES"
abcd.nback[abcd.nback$S%in%c("tar") & abcd.nback$enback_stim_acc==0 & 
             !is.na(abcd.nback$enback_stim_resp),]$R<-"NO"
abcd.nback[abcd.nback$S%in%c("lur","non") & abcd.nback$enback_stim_acc==1 & 
             !is.na(abcd.nback$enback_stim_resp),]$R<-"NO"
abcd.nback[abcd.nback$S%in%c("lur","non") & abcd.nback$enback_stim_acc==0 & 
             !is.na(abcd.nback$enback_stim_resp),]$R<-"YES"
abcd.nback$R<-factor(abcd.nback$R)

abcd.nback$RT<-NA
abcd.nback[!is.na(abcd.nback$R),]$RT<-abcd.nback[!is.na(abcd.nback$R),]$enback_stim_rt/1000

abcd.nback$acc<-NA
abcd.nback[!is.na(abcd.nback$R),]$acc<-abcd.nback[!is.na(abcd.nback$R),]$enback_stim_acc

# separate into baseline and year 2

nback.base<-abcd.nback[abcd.nback$eventname=="baseline_year_1_arm_1",]
nback.y2<-abcd.nback[abcd.nback$eventname=="2_year_follow_up_y_arm_1",]

#compute summary stats

sum.base<-data.frame(unique(nback.base$subject));colnames(sum.base)<-"s"
sum.y2<-data.frame(unique(nback.y2$subject));colnames(sum.y2)<-"s"

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
  tmp0<-nback.base[nback.base$subject==sum.base$s[r] & nback.base$enback_loadcon=="0-Back",]
  tmp2<-nback.base[nback.base$subject==sum.base$s[r] & nback.base$enback_loadcon=="2-Back",]
  
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

save(sum.base,file="sum_base_fits.RData")
load("sum_base_fits.RData")

# inclusion for each task

inc.0b<-sum.base[sum.base$acc.0.total>=.55 & sum.base$p.omit.0<=.25,"s"]
inc.2b<-sum.base[sum.base$acc.2.total>=.55 & sum.base$p.omit.2<=.25,"s"]

# final groups

base.0.dmc<-nback.base[nback.base$enback_loadcon=="0-Back" & nback.base$subject%in%inc.0b,]
base.0.dmc<-base.0.dmc[,c("subject","S","R","RT")]
colnames(base.0.dmc)<-c("s","S","R","RT")
base.0.dmc$s<-factor(base.0.dmc$s)

base.2.dmc<-nback.base[nback.base$enback_loadcon=="2-Back" & nback.base$subject%in%inc.2b,]
base.2.dmc<-base.2.dmc[,c("subject","S","R","RT")]
colnames(base.2.dmc)<-c("s","S","R","RT")
base.2.dmc$s<-factor(base.2.dmc$s)

# exclude fast guesses

base.0.dmc<-base.0.dmc[base.0.dmc$RT>=0.200 | is.na(base.0.dmc$RT),]
base.2.dmc<-base.2.dmc[base.2.dmc$RT>=0.200 | is.na(base.2.dmc$RT),]

# save out

save(base.0.dmc,file="base_0_dat_DDM.RData")
save(base.2.dmc,file="base_2_dat_DDM.RData")

#####################################################################
###### LBA estimation with broad priors #############################
#####################################################################

### load data ###

load("base_0_dat_DDM.RData")

load("base_2_dat_DDM.RData")

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

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"                "B"                "mean_v.lur.true" 
# [4] "mean_v.non.true"  "mean_v.tar.true"  "mean_v.lur.false"
# [7] "mean_v.non.false" "mean_v.tar.false" "sd_v.false"      
# [10] "t0"               "gf"    



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


zero.base.mdi <- data.model.dmc(base.0.dmc, model)

two.base.mdi <- data.model.dmc(base.2.dmc, model)

sZero.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = zero.base.mdi, thin = 5)
sZero.base.first<-sZero.base[1:4500]
save(sZero.base.first,file="zero_lba_broad_first.RData")
sZero.base.second<-sZero.base[4501:length(sZero.base)]
save(sZero.base.second,file="zero_lba_broad_second.RData")

sTwo.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = two.base.mdi, thin = 5)
sTwo.base.first<-sTwo.base[1:4500]
save(sTwo.base.first,file="two_lba_broad_first.RData")
sTwo.base.second<-sTwo.base[4501:length(sTwo.base)]
save(sTwo.base.second,file="two_lba_broad_second.RData")

# run on greatlakes 

source("greatlakes.R")
run.greatlakes.dmc("zero_lba_broad_first","LBA","lba_BpvGF.R","asweigar",
                   "asweigar0",40,wall.hours=5)

source("greatlakes.R")
run.greatlakes.dmc("zero_lba_broad_second","LBA","lba_BpvGF.R","asweigar",
                   "asweigar0",40,wall.hours=10)

source("greatlakes.R")
run.greatlakes.dmc("two_lba_broad_first","LBA","lba_BpvGF.R","asweigar",
                   "asweigar0",40,wall.hours=10)

source("greatlakes.R")
run.greatlakes.dmc("two_lba_broad_second","LBA","lba_BpvGF.R","asweigar",
                   "asweigar0",40,wall.hours=10)

#### parameter estimates

# loop to make summary parameters (EEA and SEA)


sZero.base.pars<-sZero.base

for (s in 1:length(sZero.base.pars)){
  t<-sZero.base.pars[[s]]$theta; d<-dim(t)
  gf.natural<-pnorm(t[,"gf",])
  EEA<-( (t[,"mean_v.lur.true",]-t[,"mean_v.lur.false",]) + 
           (t[,"mean_v.non.true",]-t[,"mean_v.non.false",]) +
           (t[,"mean_v.tar.true",]-t[,"mean_v.tar.false",]) )/3
  SEA<-( (t[,"mean_v.lur.true",]+t[,"mean_v.lur.false",]) + 
           (t[,"mean_v.non.true",]+t[,"mean_v.non.false",]) +
           (t[,"mean_v.tar.true",]+t[,"mean_v.tar.false",]) )/6
  D <- d + c(0, 3, 0)
  t2<-array(sapply(1:D[3], function(x) cbind(t[,,x], gf.natural[,x],
                                             EEA[,x],
                                             SEA[,x])), D)
  dimnames(t2)[[2]]<-c(dimnames(t)[[2]],"gf.natural","EEA",
                       "SEA")
  sZero.base.pars[[s]]$theta<-t2
}

sTwo.base.pars<-sTwo.base

for (s in 1:length(sTwo.base)){
  t<-sTwo.base[[s]]$theta; d<-dim(t)
  gf.natural<-pnorm(t[,"gf",])
  EEA<-( (t[,"mean_v.lur.true",]-t[,"mean_v.lur.false",]) + 
           (t[,"mean_v.non.true",]-t[,"mean_v.non.false",]) +
           (t[,"mean_v.tar.true",]-t[,"mean_v.tar.false",]) )/3
  SEA<-( (t[,"mean_v.lur.true",]+t[,"mean_v.lur.false",]) + 
           (t[,"mean_v.non.true",]+t[,"mean_v.non.false",]) +
           (t[,"mean_v.tar.true",]+t[,"mean_v.tar.false",]) )/6
  D <- d + c(0, 3, 0)
  t2<-array(sapply(1:D[3], function(x) cbind(t[,,x], gf.natural[,x],
                                             EEA[,x],
                                             SEA[,x])), D)
  dimnames(t2)[[2]]<-c(dimnames(t)[[2]],"gf.natural","EEA",
                       "SEA")
  sTwo.base.pars[[s]]$theta<-t2
}

# save out point estimates 

sZero.base_medians<-lapply(sZero.base.pars,FUN=function(x) apply(x$theta,2,median))
sZero.base_medians<-as.data.frame(t(as.data.frame(sZero.base_medians)))

sTwo.base_medians<-lapply(sTwo.base.pars,FUN=function(x) apply(x$theta,2,median))
sTwo.base_medians<-as.data.frame(t(as.data.frame(sTwo.base_medians)))

write.csv(sZero.base_medians,file="ABCD_zero_LBA_postmedians.csv")

write.csv(sTwo.base_medians,file="ABCD_two_LBA_postmedians.csv")

### simulations for model fit

sim.Zero.base <- h.post.predict.dmc(sZero.base,cores=16)
save(sim.Zero.base,file="nback_lba_pp_Zero.RData")

sim.Two.base <- h.post.predict.dmc(sTwo.base,cores=16)
save(sim.Two.base,file="nback_lba_pp_Two.RData")

### model fit plots ###

# load pp data

load("nback_lba_pp_Zero.RData")

jpeg(filename = "ABCD_0_lba.jpg",width = 9.2,height = 3,
     units = "in",res=400)
plot.pp.dmc(sim.Zero.base,layout = c(1,3),
            model.legend = F,show.fits=FALSE,mar=c(4,5,3,1),pos=NA,
            data.col = "red")
dev.off()

load("nback_lba_pp_Two.RData")

jpeg(filename = "ABCD_2_lba.jpg",width = 9.2,height = 3,
     units = "in",res=400)
plot.pp.dmc(sim.Two.base,layout = c(1,3),
            model.legend = F,show.fits=FALSE,mar=c(4,5,3,1),pos=NA,
            data.col = "red")
dev.off()

