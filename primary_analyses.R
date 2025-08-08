rm(list=ls())

#### load required packages ####
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)
library(ggpointdensity)
library(lme4)
library(partR2)

#### load HCP ####

tmp<-read.csv("HCP_zero_DDM_postmedians.csv")
colnames(tmp)<-c("X",paste0("ddm.zero.",colnames(tmp)[-1]))
dat.hcp<-tmp

tmp<-read.csv("HCP_zero_LBA_postmedians.csv")
colnames(tmp)<-c("X",paste0("lba.zero.",colnames(tmp)[-1]))
dat.hcp<-merge(dat.hcp,tmp,by="X",all=TRUE)

tmp<-read.csv("HCP_two_DDM_postmedians.csv")
colnames(tmp)<-c("X",paste0("ddm.two.",colnames(tmp)[-1]))
dat.hcp<-merge(dat.hcp,tmp,by="X",all=TRUE)

tmp<-read.csv("HCP_two_LBA_postmedians.csv")
colnames(tmp)<-c("X",paste0("lba.two.",colnames(tmp)[-1]))
dat.hcp<-merge(dat.hcp,tmp,by="X",all=TRUE)

# add Subject column
dat.hcp$Subject<-gsub("X","",dat.hcp$X)

# add information about family nesting

hcp_demo_covs<-read.csv("RESTRICTED_HCP.csv")

dat.hcp<-merge(dat.hcp,hcp_demo_covs[,c("Subject","Family_ID")],all.x=TRUE)

#### load ABCD ####

tmp<-read.csv("ABCD_zero_DDM_postmedians.csv")
colnames(tmp)<-c("X",paste0("ddm.zero.",colnames(tmp)[-1]))
dat.abcd<-tmp

tmp<-read.csv("ABCD_zero_LBA_postmedians.csv")
colnames(tmp)<-c("X",paste0("lba.zero.",colnames(tmp)[-1]))
dat.abcd<-merge(dat.abcd,tmp,by="X",all=TRUE)

tmp<-read.csv("ABCD_two_DDM_postmedians.csv")
colnames(tmp)<-c("X",paste0("ddm.two.",colnames(tmp)[-1]))
dat.abcd<-merge(dat.abcd,tmp,by="X",all=TRUE)

tmp<-read.csv("ABCD_two_LBA_postmedians.csv")
colnames(tmp)<-c("X",paste0("lba.two.",colnames(tmp)[-1]))
dat.abcd<-merge(dat.abcd,tmp,by="X",all=TRUE)

# add src_subject_id column
dat.abcd$src_subject_id<-dat.abcd$X

# add information about family nesting

abcd_y_lt<-read.csv("abcd_y_lt.csv")

dat.abcd<-merge(dat.abcd,
                abcd_y_lt[abcd_y_lt$eventname=="baseline_year_1_arm_1",
                            c("src_subject_id","site_id_l","rel_family_id")],
                all.x=TRUE)

#### load numerosity ####

tmp<-read.csv("num_DDM_postmedians.csv")
colnames(tmp)<-c("X",paste0("ddm.",colnames(tmp)[-1]))
dat.num<-tmp

tmp<-read.csv("num_LBA_postmedians.csv")
colnames(tmp)<-c("X",paste0("lba.",colnames(tmp)[-1]))
dat.num<-merge(dat.num,tmp,by="X",all=TRUE)


#### estimate unique r^2 for EEA and SEA relations with v ####

# relevant functions for computing r-squared change

rsq<-function(v,vars,dat){
  out<-c(NA,NA)
  out[1]<-(summary(lm(as.formula(paste0(v,"~",vars[1],"+",vars[2])),data = dat))$r.squared-
      summary(lm(as.formula(paste0(v,"~",vars[2])) ,data = dat))$r.squared)
  out[2]<-(summary(lm(as.formula(paste0(v,"~",vars[1],"+",vars[2])),data = dat))$r.squared-
             summary(lm(as.formula(paste0(v,"~",vars[1])) ,data = dat))$r.squared)
  out
}

rsq.mlm<-function(v,vars,dat,rand){
  tmp<-lmer(as.formula(paste(v," ~ ",paste(vars,collapse=' + '),"+",rand)),
            data = dat)
  tmp<-partR2(tmp,data=dat.abcd,partvars = vars,R2_type = "marginal")
  out<-tmp$R2[tmp$R2$term%in%vars,]$estimate
  out
}

# function for bootstrapped r^2 estimates
sum.mlm<-function(dv,vars,dat,rand){
  dat<-dat[rowSums(is.na(dat[,c(dv,vars)]))==0,]
  mod<-lmer(as.formula(paste(dv," ~ ",paste(vars,collapse=' + '),"+",rand)),
            data = dat)
  r2<-partR2(mod,data=dat,partvars = vars,R2_type = "marginal",nboot=1000)
  out<-list(r2=r2,mod=mod,n=length(dat[,1]))
  out
}



# compute r-squared with standard regression and multi-level models 

u_rsq<-data.frame(task=rep(c("HCP 0-back","HCP 2-back",
                  "ABCD 0-back", "ABCD 2-back",
                  "Num"),2),var=factor(c(rep("EEA",5),rep("SEA",5))),
                  rs=NA)

u_rsq[u_rsq$task=="HCP 0-back","rs"]<-rsq(v="ddm.zero.overall.v",vars=c("lba.zero.EEA","lba.zero.SEA"),dat=dat.hcp)
u_rsq[u_rsq$task=="HCP 2-back","rs"]<-rsq(v="ddm.two.overall.v",vars=c("lba.two.EEA","lba.two.SEA"),dat=dat.hcp)
u_rsq[u_rsq$task=="ABCD 0-back","rs"]<-rsq(v="ddm.zero.overall.v",vars=c("lba.zero.EEA","lba.zero.SEA"),dat=dat.abcd)
u_rsq[u_rsq$task=="ABCD 2-back","rs"]<-rsq(v="ddm.two.overall.v",vars=c("lba.two.EEA","lba.two.SEA"),dat=dat.abcd)
u_rsq[u_rsq$task=="Num","rs"]<-rsq(v="ddm.overall.v",vars=c("lba.EEA","lba.SEA"),dat=dat.num)

u_rsq$rs<-round(u_rsq$rs,3)


u_rsq_f<-u_rsq

u_rsq_f[u_rsq_f$task=="HCP 0-back","rs"]<-rsq.mlm(v="ddm.zero.overall.v",vars=c("lba.zero.EEA","lba.zero.SEA"),
                                                  dat=dat.hcp,rand="(1|Family_ID)")
u_rsq_f[u_rsq_f$task=="HCP 2-back","rs"]<-rsq.mlm(v="ddm.two.overall.v",vars=c("lba.two.EEA","lba.two.SEA"),
                                                  dat=dat.hcp,rand="(1|Family_ID)")
u_rsq_f[u_rsq_f$task=="ABCD 0-back","rs"]<-rsq.mlm(v="ddm.zero.overall.v",vars=c("lba.zero.EEA","lba.zero.SEA"),
                                                   dat=dat.abcd,rand="(1|site_id_l/rel_family_id)")
u_rsq_f[u_rsq_f$task=="ABCD 2-back","rs"]<-rsq.mlm(v="ddm.two.overall.v",vars=c("lba.two.EEA","lba.two.SEA"),
                                                   dat=dat.abcd,rand="(1|site_id_l/rel_family_id)")

u_rsq_f$rs<-round(u_rsq_f$rs,3)


#### main figure #####

hcp.zero.EEA <- ggplot(data = dat.hcp, mapping = aes(x = lba.zero.EEA,dat, y = ddm.zero.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA EEA") + ylab("DDM drift rate") + xlim(c(0,4)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=.5,y=3.5,label=round(cor(dat.hcp[,c("lba.zero.EEA","ddm.zero.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

hcp.zero.SEA <- ggplot(data = dat.hcp, mapping = aes(x = lba.zero.SEA,dat, y = ddm.zero.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA SEA") + ylab("DDM drift rate") + xlim(c(1.5,6)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=2,y=3.5,label=round(cor(dat.hcp[,c("lba.zero.SEA","ddm.zero.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

hcp.zero.rsq <-ggplot(data=u_rsq_f[u_rsq_f$task=="HCP 0-back",],aes(x=var,y=rs, fill = var)) + geom_bar(stat = "identity",show.legend = FALSE) + 
  theme_classic(base_size = 14) + xlab("") + ylab("variance explained") +
  ylim(c(0,1)) +  geom_text(aes(y = rs, label = rs), vjust = -1) 


hcp.two.EEA <- ggplot(data = dat.hcp, mapping = aes(x = lba.two.EEA,dat, y = ddm.two.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA EEA") + ylab("DDM drift rate") + xlim(c(0,4)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=.5,y=3.5,label=round(cor(dat.hcp[,c("lba.two.EEA","ddm.two.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

hcp.two.SEA <- ggplot(data = dat.hcp, mapping = aes(x = lba.two.SEA,dat, y = ddm.two.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA SEA") + ylab("DDM drift rate") + xlim(c(1.5,6)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=2,y=3.5,label=round(cor(dat.hcp[,c("lba.two.SEA","ddm.two.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

hcp.two.rsq <-ggplot(data=u_rsq_f[u_rsq_f$task=="HCP 2-back",],aes(x=var,y=rs, fill = var)) + geom_bar(stat = "identity",show.legend = FALSE) + 
  theme_classic(base_size = 14) + xlab("") + ylab("variance explained") +
  ylim(c(0,1)) +  geom_text(aes(y = rs, label = rs), vjust = -1) 


abcd.zero.EEA <- ggplot(data = dat.abcd, mapping = aes(x = lba.zero.EEA,dat, y = ddm.zero.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA EEA") + ylab("DDM drift rate") + xlim(c(0,4)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=.5,y=3.5,label=round(cor(dat.abcd[,c("lba.zero.EEA","ddm.zero.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

abcd.zero.SEA <- ggplot(data = dat.abcd, mapping = aes(x = lba.zero.SEA,dat, y = ddm.zero.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA SEA") + ylab("DDM drift rate") + xlim(c(1.5,6)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=2,y=3.5,label=round(cor(dat.abcd[,c("lba.zero.SEA","ddm.zero.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

abcd.zero.rsq <-ggplot(data=u_rsq_f[u_rsq_f$task=="ABCD 0-back",],aes(x=var,y=rs, fill = var)) + geom_bar(stat = "identity",show.legend = FALSE) + 
  theme_classic(base_size = 14) + xlab("") + ylab("variance explained") +
  ylim(c(0,1)) +  geom_text(aes(y = rs, label = rs), vjust = -1) 


abcd.two.EEA <- ggplot(data = dat.abcd, mapping = aes(x = lba.two.EEA,dat, y = ddm.two.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA EEA") + ylab("DDM drift rate") + xlim(c(0,4)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=.5,y=3.5,label=round(cor(dat.abcd[,c("lba.two.EEA","ddm.two.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

abcd.two.SEA <- ggplot(data = dat.abcd, mapping = aes(x = lba.two.SEA,dat, y = ddm.two.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA SEA") + ylab("DDM drift rate") + xlim(c(1.5,6)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=2,y=3.5,label=round(cor(dat.abcd[,c("lba.two.SEA","ddm.two.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

abcd.two.rsq <-ggplot(data=u_rsq_f[u_rsq_f$task=="ABCD 2-back",],aes(x=var,y=rs, fill = var)) + geom_bar(stat = "identity",show.legend = FALSE) + 
  theme_classic(base_size = 14) + xlab("") + ylab("variance explained") +
  ylim(c(0,1)) +  geom_text(aes(y = rs, label = rs), vjust = -1) 



num.EEA <- ggplot(data = dat.num, mapping = aes(x = lba.EEA,dat, y = ddm.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA EEA") + ylab("DDM drift rate") + xlim(c(0,4)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=.5,y=3.5,label=round(cor(dat.num[,c("lba.EEA","ddm.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

num.SEA <- ggplot(data = dat.num, mapping = aes(x = lba.SEA,dat, y = ddm.overall.v)) +
  geom_pointdensity(show.legend = FALSE) + xlab("LBA SEA") + ylab("DDM drift rate") + xlim(c(1.5,6)) + ylim(0,4) +
  scale_color_viridis() + theme_classic() + 
  geom_text(x=2,y=3.5,label=round(cor(dat.num[,c("lba.SEA","ddm.overall.v")],use="complete")[1,2],2),cex=5) +
  geom_smooth(method='lm', formula= y~x,col="red", se=FALSE)

num.rsq <-ggplot(data=u_rsq_f[u_rsq_f$task=="Num",],aes(x=var,y=rs, fill = var)) + geom_bar(stat = "identity",show.legend = FALSE) + 
  theme_classic(base_size = 14) + xlab("") + ylab("variance explained") +
  ylim(c(0,1)) +  geom_text(aes(y = rs, label = rs), vjust = -1) 

# merge together for figure 

Fig2<-grid.arrange(hcp.zero.EEA,hcp.zero.SEA,hcp.zero.rsq,
             hcp.two.EEA,hcp.two.SEA,hcp.two.rsq,
             abcd.zero.EEA,abcd.zero.SEA,abcd.zero.rsq,
             abcd.two.EEA,abcd.two.SEA,abcd.two.rsq,
             num.EEA,num.SEA,num.rsq,
             nrow=5,
                left= textGrob(paste("Numerosity",
                                     "ABCD 2-back","ABCD 0-back",
                                     "HCP 2-back", "HCP 0-back",
                                sep = "            "),
               gp = gpar(fontsize = 18),x=.4,y=.50,rot=90))


ggsave(filename = "Figure_2.tif",Fig1,units = "in",dpi=300,
       width = 7,height = 11.5,device="tiff")


#### HCP criterion measure analyses ####

# read in data files 

hcp_cog<-read.csv("unrestricted_HCP.csv")


hcp_covs<-c("CogTotalComp_Unadj" , "CogTotalComp_AgeAdj",
            "ListSort_Unadj" , "ListSort_AgeAdj",
            "CardSort_Unadj" , "CardSort_AgeAdj" ,
            "Flanker_Unadj", "Flanker_AgeAdj" )

dat.hcp<-merge(dat.hcp,hcp_cog[,c("Subject",hcp_covs)],all.x=TRUE)


Total.hcp.lba.0<-sum.mlm(dv="CogTotalComp_Unadj",dat=dat.hcp,rand="(1|Family_ID)",
                         vars=c("lba.zero.EEA","lba.zero.SEA"))
Total.hcp.lba.2<-sum.mlm(dv="CogTotalComp_Unadj",dat=dat.hcp,rand="(1|Family_ID)",
                         vars=c("lba.two.EEA","lba.two.SEA"))



List.hcp.lba.0<-sum.mlm(dv="ListSort_Unadj",dat=dat.hcp,rand="(1|Family_ID)",
                         vars=c("lba.zero.EEA","lba.zero.SEA"))
List.hcp.lba.2<-sum.mlm(dv="ListSort_Unadj",dat=dat.hcp,rand="(1|Family_ID)",
                         vars=c("lba.two.EEA","lba.two.SEA"))



Card.hcp.lba.0<-sum.mlm(dv="CardSort_Unadj",dat=dat.hcp,rand="(1|Family_ID)",
                        vars=c("lba.zero.EEA","lba.zero.SEA"))
Card.hcp.lba.2<-sum.mlm(dv="CardSort_Unadj",dat=dat.hcp,rand="(1|Family_ID)",
                        vars=c("lba.two.EEA","lba.two.SEA"))


Flank.hcp.lba.0<-sum.mlm(dv="Flanker_Unadj",dat=dat.hcp,rand="(1|Family_ID)",
                        vars=c("lba.zero.EEA","lba.zero.SEA"))
Flank.hcp.lba.2<-sum.mlm(dv="Flanker_Unadj",dat=dat.hcp,rand="(1|Family_ID)",
                        vars=c("lba.two.EEA","lba.two.SEA"))


#### ABCD criterion measure analyses ####

# read in data files 

nc_y_nihtb<-read.csv("nc_y_nihtb.csv")


abcd_covs<-c("nihtbx_totalcomp_uncorrected" , "nihtbx_totalcomp_agecorrected",
             "nihtbx_list_uncorrected" , "nihtbx_list_agecorrected",
             "nihtbx_cardsort_uncorrected" , "nihtbx_cardsort_agecorrected" ,
             "nihtbx_flanker_uncorrected", "nihtbx_flanker_agecorrected")

dat.abcd<-merge(dat.abcd,nc_y_nihtb[nc_y_nihtb$eventname=="baseline_year_1_arm_1",c("src_subject_id",abcd_covs)],all.x=TRUE)



Total.abcd.lba.0<-sum.mlm(dv="nihtbx_totalcomp_uncorrected",dat=dat.abcd,rand="(1|site_id_l/rel_family_id)",
                          vars=c("lba.zero.EEA","lba.zero.SEA"))
Total.abcd.lba.2<-sum.mlm(dv="nihtbx_totalcomp_uncorrected",dat=dat.abcd,rand="(1|site_id_l/rel_family_id)",
                          vars=c("lba.two.EEA","lba.two.SEA"))

List.abcd.lba.0<-sum.mlm(dv="nihtbx_list_uncorrected",dat=dat.abcd,rand="(1|site_id_l/rel_family_id)",
                         vars=c("lba.zero.EEA","lba.zero.SEA"))
List.abcd.lba.2<-sum.mlm(dv="nihtbx_list_uncorrected",dat=dat.abcd,rand="(1|site_id_l/rel_family_id)",
                         vars=c("lba.two.EEA","lba.two.SEA"))

Card.abcd.lba.0<-sum.mlm(dv="nihtbx_cardsort_uncorrected",dat=dat.abcd,rand="(1|site_id_l/rel_family_id)",
                         vars=c("lba.zero.EEA","lba.zero.SEA"))
Card.abcd.lba.2<-sum.mlm(dv="nihtbx_cardsort_uncorrected",dat=dat.abcd,rand="(1|site_id_l/rel_family_id)",
                         vars=c("lba.two.EEA","lba.two.SEA"))

Flank.abcd.lba.0<-sum.mlm(dv="nihtbx_flanker_uncorrected",dat=dat.abcd,rand="(1|site_id_l/rel_family_id)",
                          vars=c("lba.zero.EEA","lba.zero.SEA"))
Flank.abcd.lba.2<-sum.mlm(dv="nihtbx_flanker_uncorrected",dat=dat.abcd,rand="(1|site_id_l/rel_family_id)",
                          vars=c("lba.two.EEA","lba.two.SEA"))


### summary table ###

sum_r2<-rbind(
  t(as.numeric(c(Total.hcp.lba.0$r2$R2[2,2:4],Total.hcp.lba.0$r2$R2[3,2:4]))),
  t(as.numeric(c(Total.hcp.lba.2$r2$R2[2,2:4],Total.hcp.lba.2$r2$R2[3,2:4]))),
  t(as.numeric(c(List.hcp.lba.0$r2$R2[2,2:4],List.hcp.lba.0$r2$R2[3,2:4]))),
  t(as.numeric(c(List.hcp.lba.2$r2$R2[2,2:4],List.hcp.lba.2$r2$R2[3,2:4]))),
  t(as.numeric(c(Card.hcp.lba.0$r2$R2[2,2:4],Card.hcp.lba.0$r2$R2[3,2:4]))),
  t(as.numeric(c(Card.hcp.lba.2$r2$R2[2,2:4],Card.hcp.lba.2$r2$R2[3,2:4]))),
  t(as.numeric(c(Flank.hcp.lba.0$r2$R2[2,2:4],Flank.hcp.lba.0$r2$R2[3,2:4]))),
  t(as.numeric(c(Flank.hcp.lba.2$r2$R2[2,2:4],Flank.hcp.lba.2$r2$R2[3,2:4]))),
  t(as.numeric(c(Total.abcd.lba.0$r2$R2[2,2:4],Total.abcd.lba.0$r2$R2[3,2:4]))),
  t(as.numeric(c(Total.abcd.lba.2$r2$R2[2,2:4],Total.abcd.lba.2$r2$R2[3,2:4]))),
  t(as.numeric(c(List.abcd.lba.0$r2$R2[2,2:4],List.abcd.lba.0$r2$R2[3,2:4]))),
  t(as.numeric(c(List.abcd.lba.2$r2$R2[2,2:4],List.abcd.lba.2$r2$R2[3,2:4]))),
  t(as.numeric(c(Card.abcd.lba.0$r2$R2[2,2:4],Card.abcd.lba.0$r2$R2[3,2:4]))),
  t(as.numeric(c(Card.abcd.lba.2$r2$R2[2,2:4],Card.abcd.lba.2$r2$R2[3,2:4]))),
  t(as.numeric(c(Flank.abcd.lba.0$r2$R2[2,2:4],Flank.abcd.lba.0$r2$R2[3,2:4]))),
  t(as.numeric(c(Flank.abcd.lba.2$r2$R2[2,2:4],Flank.abcd.lba.2$r2$R2[3,2:4]))) )

sum_r2<-as.data.frame(round(sum_r2,3))

sum_r2<-cbind(data.frame(sample=c(rep("HCP",8),rep("ABCD",8)),
                         measure=rep(c("Total Cognition", "Total Cognition",
                                       "List (working memory)", "List (working memory)",
                                       "Card Sort (flexibility)","Card Sort (flexibility)",
                                       "Flanker (inhibition)", "Flanker (inhibition)"),2),
                         task=rep(c("0-back", "2-back"),8)),
              sum_r2)

write.csv(sum_r2,file="Table_1_final.csv",row.names=F)