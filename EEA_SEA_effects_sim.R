
##################
### Overview #####
##################


# The first part of this script goes through multiple simulations 
# across the diffusion model, linear ballistic accumulator model, and 
# racing diffusion model to assess the impact of changes in EEA and SEA 
# on observed summary statistics. We encourage readers to explore this 
# code by generating their own parameter values. Overall, these simulations 
# indicate that increases in both v from the diffusion model and EEA from accumulator
# models always lead to higher accuracy and have subtle effects on RT. In  
# contrast, increases in SEA lead to large decreases in RT but have subtle
# effects on accuracy.

# The second part of the script conducts a more extensive LBA simulation
# study that varies EEA and SEA across wide ranges and plots the results.
# This study confirms that varying EEA has a larger magnitude of impact on 
# accuracy than on RT while varying SEA has a larger magnitude of impact on
# RT than on accuracy.


### first, load EMC2 package:

#remotes::install_github("ampl-psych/EMC2",ref="dev")

library(EMC2)

#######################################################################
### Part 1: explore drift rate effects on RT/accuracy across models ###
#######################################################################

#### simple (wiener) diffusion model ####

Smat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,"d"))
dWDM <- design(model=DDM,contrasts=list(S=Smat),
  factors=list(subjects="1",S=c("left","right"),D=c("easy","hard")),Rlevels=c("left","right"),
  formula=list(a~1,v~D/S,Z~1,t0~1),constants=c(v=0,v_Dhard=0))

pvWDM <- sampled_pars(dWDM)

# parameter values:
pvWDM[1:5] <- c(log(1),1.5,1,0,log(.3))

# A nice way to see the design:
plot(dWDM,pvWDM,factors=list(v=c("S","D")))

# simulate data:
datWDM <- make_data(pvWDM,dWDM,n_trials=1e4)

#  around 77% accuracy
tapply(datWDM$S==datWDM$R,datWDM[,c("D","S")],mean)

# Error and correct equal, both slowing with lower rates by 20ms
tapply(datWDM$rt,cbind(datWDM[,c("D","S")],C=datWDM$S==datWDM$R),mean)


#### full diffusion model ####

# induce slow errors by adding drift variability (sv)
dDDM <- design(model=DDM,contrasts=list(S=Smat),
  factors=list(subjects="1",S=c("left","right"),D=c("easy","hard")),Rlevels=c("left","right"),
  formula=list(a~1,v~D/S,Z~1,t0~1,sv~1),constants=c(v=0,v_Dhard=0))

pvDDM <- sampled_pars(dDDM)

# parameter values:
# Higher rates to roughly compensate for sv
pvDDM[1:6] <- c(log(1),1.75,1,0,log(.3),log(1))

plot(dDDM,pvDDM,factors=list(v=c("S","D")))

# simulate data
datDDM <- make_data(pvDDM,dDDM,n_trials=1e4)

tapply(datDDM$S==datDDM$R,datDDM[,c("D","S")],mean)

# Slower errors overall, but there is consistent slowing due to difficulty 
# across both corrects and errors
tapply(datDDM$rt,cbind(datDDM[,c("D","S")],C=datDDM$S==datDDM$R),mean)


# induce fast errors by adding start point variability (sz)
dDDM <- design(model=DDM,contrasts=list(S=Smat),
  factors=list(subjects="1",S=c("left","right"),D=c("easy","hard")),Rlevels=c("left","right"),
  formula=list(a~1,v~D/S,Z~1,t0~1,SZ~1),constants=c(v=0,v_Dhard=0))

pvDDM <- sampled_pars(dDDM)

# parameter values:
pvDDM[1:6] <- c(log(1),1.5,.9,0,log(.3),qnorm(.5))

plot(dDDM,pvDDM,factors=list(v=c("S","D")))

# simulate data
datDDM <- make_data(pvDDM,dDDM,n_trials=1e4)

tapply(datDDM$S==datDDM$R,datDDM[,c("D","S")],mean)

# Fast errors than corrects, still slowing with hard/lower rates for both
tapply(datDDM$rt,cbind(datDDM[,c("D","S")],C=datDDM$S==datDDM$R),mean)

# Slow errors, but now hard driven mainly by more variable sv gets fast hard  
dDDM <- design(model=DDM,contrasts=list(S=Smat),
  factors=list(subjects="1",S=c("left","right"),D=c("easy","hard")),Rlevels=c("left","right"),
  formula=list(a~1,v~D/S,Z~1,t0~1,sv~D),constants=c(v=0,v_Dhard=0))

pvDDM <- sampled_pars(dDDM)

# parameter values:
# Higher rates to roughly compensate for sv
pvDDM[1:7] <- c(log(1),1.75,1.5,0,log(.3),log(1),log(2))

plot(dDDM,pvDDM,factors=list(v=c("S","D")))
mapped_pars(dDDM,pvDDM)

#simulate data
datDDM <- make_data(pvDDM,dDDM,n_trials=1e4)

tapply(datDDM$S==datDDM$R,datDDM[,c("D","S")],mean)

# Slow errors, slight speeding with lower rates for correct and error
tapply(datDDM$rt,cbind(datDDM[,c("D","S")],C=datDDM$S==datDDM$R),mean)

# So overall a bit faster for hard.
tapply(datDDM$rt,datDDM[,c("D","S")],mean)




#### linear ballistic accumulator model (LBA) ####

ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))

# Only EEA no SEA difference
dLBA <- design(model=LBA,contrasts=list(lM=ADmat),matchfun=\(d)d$S==d$lR,
  factors=list(subjects="1",S=c("left","right"),D=c("easy","hard")),Rlevels=c("left","right"),
  formula=list(v~D/lM,B~1,t0~1,A~1,sv~1),constants=c(sv=0,v_Dhard=0))

pvLBA <- sampled_pars(dLBA)

# parameter values
pvLBA[1:6] <- c(1,1,.4,log(.5),log(.3),log(.15))

mapped_pars(dLBA,pvLBA)

plot(dLBA,pvLBA,factors=list('v'=c("D", "lM")),plot_factor="D")

# simulate data
datLBA <- make_data(pvLBA,dLBA,n_trials=1e4)

# ~10% accuracy difference, around 65% overall
tapply(datLBA$S==datLBA$R,datLBA[,c("D","S")],mean)

# Slow errors, slowing with lower rates for correct but OPPOSITE for errors
tapply(datLBA$rt,cbind(datLBA[,c("D","S")],C=datLBA$S==datLBA$R),median)

# Overall balances out a little so fairly small speeding for easy, but still there
tapply(datLBA$rt,datLBA[,c("D","S")],median)

# Only SEA and no EEA difference, so "hard" has overall lower rate
dLBA <- design(model=LBA,contrasts=list(lM=ADmat),matchfun=\(d)d$S==d$lR,
  factors=list(subjects="1",S=c("left","right"),D=c("easy","hard")),Rlevels=c("left","right"),
  formula=list(v~D*lM,B~1,t0~1,A~1,sv~1),constants=c(sv=0,'v_Dhard:lMd'=0))

pvLBA <- sampled_pars(dLBA)

# parameter values
pvLBA[1:6] <- c(4,-1,1,log(1),log(.3),log(.25))

mapped_pars(dLBA,pvLBA)

plot(dLBA,pvLBA,factors=list('v'=c("D", "lM")),plot_factor="D")

# simulate data
datLBA <- make_data(pvLBA,dLBA,n_trials=1e4)

# Virtually no accuracy difference
tapply(datLBA$S==datLBA$R,datLBA[,c("D","S")],mean)

# Slow errors and slowing with difficulty across both corrects and errors
tapply(datLBA$rt,cbind(datLBA[,c("D","S")],C=datLBA$S==datLBA$R),median)

# large amount of speeding for easy
tapply(datLBA$rt,datLBA[,c("D","S")],median)


# Can get faster hard than easy with variability
# Only EEA no SEA difference, but add sv(hard) > sv(easy)
dLBA <- design(model=LBA,contrasts=list(lM=ADmat),matchfun=\(d)d$S==d$lR,
  factors=list(subjects="1",S=c("left","right"),D=c("easy","hard")),Rlevels=c("left","right"),
  formula=list(v~D/lM,B~1,t0~1,A~1,sv~D),constants=c(sv=0,v_Dhard=0))

pvLBA <- sampled_pars(dLBA)

# parameter values
pvLBA[1:7] <- c(1,1,.5,log(.5),log(.3),log(.15),log(1.5))
mapped_pars(dLBA,pvLBA)

# simulate data
datLBA <- make_data(pvLBA,dLBA,n_trials=1e4)

# ~15% accuracy difference, around 64% overall
tapply(datLBA$S==datLBA$R,datLBA[,c("D","S")],mean)

# Overall now hard is fast!
tapply(datLBA$rt,cbind(datLBA[,c("D","S")],C=datLBA$S==datLBA$R),median)
tapply(datLBA$rt,datLBA[,c("D","S")],median)


#### racing diffusion model (RDM) ####

ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))

# Only EEA no SEA difference
dRDM <- design(model=RDM,contrasts=list(lM=ADmat),matchfun=\(d)d$S==d$lR,
  factors=list(subjects="1",S=c("left","right"),D=c("easy","hard")),Rlevels=c("left","right"),
  formula=list(v~D*lM,B~1,t0~1,s~1),constants=c(s=0))

# A bit harder to wrangle because of the log on v
mapped_pars(dRDM)

dSEA <- function(v,dve,dvh,du) {
  eT <- exp(v + 0.5 * dve)
  eF <- exp(v - 0.5 * dve) 
  hT <- exp(v + du + 0.5 * dve + 0.5 * dvh)
  hF <- exp(v + du - 0.5 * dve - 0.5 * dvh)
  c(eT=eT,eF=eF,hT=hT,hF=hF,eEEA=eT-eF,hEEA=hT-hF,eSEA=eT+eF,hSEA=hT+hF)
}

round(dSEA(v=.5,du=.1,dve=1,dvh=-.5),2)

pvRDM <- sampled_pars(dRDM)

# parameter values
pvRDM[1:6] <- c(.5,.1,1,-.5,log(1),log(.3))

mapped_pars(dRDM,pvRDM)
plot(dRDM,pvRDM,factors=list('v'=c("D", "lM")),plot_factor="D")

# simulate data:
datRDM <- make_data(pvRDM,dRDM,n_trials=1e4)

# 10% accuracy differnece around 70%
tapply(datRDM$S==datRDM$R,datRDM[,c("D","S")],mean)

# Errors slightly slow, slowing with lower rates for correct and errors
tapply(datRDM$rt,cbind(datRDM[,c("D","S")],C=datRDM$S==datRDM$R),median)
# Overall slight hard slowing
tapply(datRDM$rt,datRDM[,c("D","S")],median)


# Only SEA, no EEA, difference, so "hard" has lower overall rate
dRDM <- design(model=RDM,contrasts=list(lM=ADmat),matchfun=\(d)d$S==d$lR,
  factors=list(subjects="1",S=c("left","right"),D=c("easy","hard")),Rlevels=c("left","right"),
  formula=list(v~D*lM,B~1,t0~1,s~1),constants=c(s=0))

round(dSEA(v=.5,du=-.5,dve=1,dvh=.56),2)

pvRDM <- sampled_pars(dRDM)

# parameter values
pvRDM[1:6] <- c(.5,-.5,1,.56,log(1),log(.3))

mapped_pars(dRDM,pvRDM)

plot(dRDM,pvRDM,factors=list('v'=c("D", "lM")),plot_factor="D")

# simulate data:
datRDM <- make_data(pvRDM,dRDM,n_trials=1e4)

# Slighlty MORE accurate for hard
tapply(datRDM$S==datRDM$R,datRDM[,c("D","S")],mean)

# Errors and corrects similar , slowing with lower rates for correct and errors
tapply(datRDM$rt,cbind(datRDM[,c("D","S")],C=datRDM$S==datRDM$R),median)

# Overall clear hard slowing
tapply(datRDM$rt,datRDM[,c("D","S")],median)


#######################################################################
### Part 2: Parametric effects of EEA and SEA on RT/ accuracy #########
#######################################################################

### define base LBA model with only EEA and SEA values that vary

ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))

dLBA <- design(model=LBA,contrasts=list(lM=ADmat),matchfun=\(d)d$S==d$lR,
               factors=list(subjects="1",S=c("left","right")),Rlevels=c("left","right"),
               formula=list(v~lM,B~1,t0~1,A~1,sv~1),
               constants=c(B=log(.5),t0=log(.3),A=log(.15),sv=1))


### matrices of simulated parameter values, fix other parameter at 2 

EEA.par.mat<-data.frame(v=rep(2,41),v_lMd=seq(0,4,.1))
SEA.par.mat<-data.frame(v=seq(0,4,.1),v_lMd=rep(2,41))

### simulate data across all values

EEA.sim<-list()
for (i in 1:length(EEA.par.mat$v_lMd)){
  tmp<-make_data(EEA.par.mat[i,],dLBA,n_trials=1e4)
  tmp$C<-tmp$S==tmp$R
  tmp$EEA<-EEA.par.mat[i,"v_lMd"]
  EEA.sim[[i]]<-tmp
}
EEA.sim<-do.call("rbind",EEA.sim)


SEA.sim<-list()
for (i in 1:length(SEA.par.mat$v_lMd)){
  tmp<-make_data(SEA.par.mat[i,],dLBA,n_trials=1e4)
  tmp$C<-tmp$S==tmp$R
  tmp$SEA<-SEA.par.mat[i,"v"]
  SEA.sim[[i]]<-tmp
}
SEA.sim<-do.call("rbind",SEA.sim)

### plot out 

# accuracy: goes up dramtically with EEA, but also up slightly with SEA (not sure why?)
tapply(EEA.sim$C,EEA.sim$EEA,mean)
tapply(SEA.sim$C,SEA.sim$SEA,mean)

# overall RT: goes way down with SEA, slight decrease with EEA
tapply(EEA.sim$rt,EEA.sim$EEA,mean)
tapply(SEA.sim$rt,SEA.sim$SEA,mean)

# correct RT: goes way down with SEA, slight decrease with EEA
tapply(EEA.sim[EEA.sim$C,]$rt,EEA.sim[EEA.sim$C,]$EEA,mean)
tapply(SEA.sim[SEA.sim$C,]$rt,SEA.sim[SEA.sim$C,]$SEA,mean)

# incorrect RT: way down with SEA, small reverse effect with EEA
tapply(EEA.sim[!EEA.sim$C,]$rt,EEA.sim[!EEA.sim$C,]$EEA,mean)
tapply(SEA.sim[!SEA.sim$C,]$rt,SEA.sim[!SEA.sim$C,]$SEA,mean)

# plot all 

jpeg("EEA_SEA_supp.jpg",7,8,units = 'in',res=300)

par(mfrow=c(2,2))

plot(EEA.par.mat$v_lMd,tapply(EEA.sim$C,EEA.sim$EEA,mean),col="red",pch=16,
     xlab="EEA (red) / SEA (blue)",ylab="accuracy")
points(SEA.par.mat$v,tapply(SEA.sim$C,SEA.sim$SEA,mean),col="blue",pch=16)

plot(EEA.par.mat$v_lMd,tapply(EEA.sim$rt,EEA.sim$EEA,mean),col="red",pch=16,
     xlab="EEA (red) / SEA (blue)",ylab="overall RT",ylim=c(0.4,0.6))
points(SEA.par.mat$v,tapply(SEA.sim[SEA.sim$C,]$rt,SEA.sim[SEA.sim$C,]$SEA,mean),col="blue",pch=16)

plot(EEA.par.mat$v_lMd,tapply(EEA.sim[EEA.sim$C,]$rt,EEA.sim[EEA.sim$C,]$EEA,mean),col="red",pch=16,
     xlab="EEA (red) / SEA (blue)",ylab="correct RT",ylim=c(0.4,0.6))
points(SEA.par.mat$v,tapply(SEA.sim[SEA.sim$C,]$rt,SEA.sim[SEA.sim$C,]$SEA,mean),col="blue",pch=16)

plot(EEA.par.mat$v_lMd,tapply(EEA.sim[!EEA.sim$C,]$rt,EEA.sim[!EEA.sim$C,]$EEA,mean),col="red",pch=16,
     xlab="EEA (red) / SEA (blue)",ylab="incorrect RT",ylim=c(0.4,0.6))
points(SEA.par.mat$v,tapply(SEA.sim[SEA.sim$C,]$rt,SEA.sim[SEA.sim$C,]$SEA,mean),col="blue",pch=16)

dev.off()

