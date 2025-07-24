##############################################################################/
##############################################################################/
#ACCase : Test of the effect of several herbicides on an aphid species
##############################################################################/
##############################################################################/

#loading the data and packages
source("CrossKingdom_load.R")


##############################################################################/
#Effect on fertility by lm####
##############################################################################/

#analysis of the effect of SA on the mean number of larvae
temp<-MyzHerbi[MyzHerbi$Dose!="N/2",]
temp<-drop.levels(temp,reorder=FALSE)
mod_nblarv<-aov(Total~Active_substance*Dose*Clone
                +Error(Repetition),data=temp)
summary(mod_nblarv)
#we can remove the three way interaction
mod_nblarv<-aov(Total~(Active_substance+Dose+Clone)^2
                +Error(Repetition),data=temp)
summary(mod_nblarv)

emm<-emmeans(mod_nblarv,~(Active_substance+Dose+Clone)^2)
plot(emm)
emm<-emmeans(mod_nblarv,~(Clone+Dose+Active_substance)^2)
plot(emm)

#using the afex package
#first we need an ID for the different individual
temp$ID<-paste(temp$Clone,temp$Active_substance,temp$Dose,sep="_")
mod_afex<-aov_car(Total~((Active_substance+Dose+Clone)^2+
                           Error(ID/Repetition)),
                  data=temp)
mod_afex
mod_afex<-aov_car(Total~(Active_substance+Dose+Clone+
                           Active_substance:Dose+Dose:Clone+
                           Error(ID/Repetition)),
                  data=temp)
mod_afex
afex_plot(mod_afex,"Dose","Active_substance","Clone")
afex_plot(mod_afex,"Dose",panel=~Clone+Active_substance)
#no strong Clone effect
afex_plot(mod_afex,"Dose",panel=~Active_substance)


emmAS<-emmeans(mod_afex,~Active_substance)
pairs(emmAS)
afex_plot(mod_afex,"Active_substance")
summary(as.glht(pairs(emmAS)),test=adjusted("BH"))
emmDos<-emmeans(mod_afex,~Dose)
pairs(emmDos)
afex_plot(mod_afex,"Dose")
summary(as.glht(pairs(emmDos)),test=adjusted("BH"))
emmClone<-emmeans(mod_afex,~Clone)
pairs(emmClone)
afex_plot(mod_afex,"Clone")
summary(as.glht(pairs(emmClone)),test=adjusted("BH"))

#investigating the effect of Dose by active substance
emmDoByAS<-emmeans(mod_afex,"Dose",by="Active_substance")
emmDoByAS
pairs(emmDoByAS)


##############################################################################/
#Effect on fertility using glmm####
##############################################################################/

#same kind of analysis but in a glmm framework
mmod_nblarv<-lme(Total~Active_substance*Dose*Clone,
                 random= ~1|Repetition,
                 data=temp,method="ML")
summary(mmod_nblarv)
mmod_nblarv.1<-update(mmod_nblarv,~. -Active_substance:Dose:Clone)
summary(mmod_nblarv.1)
anova(mmod_nblarv.1,mmod_nblarv)
mmod_nblarv.2<-update(mmod_nblarv.1,~. -Active_substance:Clone)
anova(mmod_nblarv.2,mmod_nblarv.1)
#the interaction Active_substance:Clone could be removed
mmod_nblarv.3<-update(mmod_nblarv.2,~. -Dose:Clone)
anova(mmod_nblarv.3,mmod_nblarv.2)
mmod_nblarv.4<-update(mmod_nblarv.2,~. -Active_substance:Dose)
anova(mmod_nblarv.4,mmod_nblarv.2)

#here is the model we keep
summary(mmod_nblarv.1)
#checking the residuals
plot(mmod_nblarv.1)
#response variable vs the fitted value
plot(mmod_nblarv.1,Total~fitted(.))
#normal distribution of errors in the different repetition
qqnorm(mmod_nblarv.1,~resid(.)|Repetition)
Anova(mmod_nblarv.1,type=c("III")) #type III anova
AIC(mmod_nblarv.1)


#same model but using lmer
mmod_nblarv.1bis<-lmer(Total~(Active_substance+Dose+Clone)^2+(1|Repetition),
                       data=temp,REML=FALSE)
plot(mmod_nblarv.1bis)
plot(mmod_nblarv.1bis,Total~fitted(.))
Anova(mmod_nblarv.1bis,type="III")
testDispersion(mmod_nblarv.1bis) #no overdispersion
plotResiduals(mmod_nblarv.1bis)
resid<-simulateResiduals(fittedModel=mmod_nblarv.1bis)
plot(resid)
AIC(mmod_nblarv.1bis)
#checking for multicollinearity
vif(mmod_nblarv.1bis) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10))

#same analysis with package glmmTMB in order to use
#generalized model: glmm with poisson distribution
mmod_nblarv.1ter<-glmmTMB(Total~(Active_substance+Dose+Clone)^2+
                            (1|Repetition),
                          data=temp,REML=FALSE,family=poisson)
summary(mmod_nblarv.1ter)
Anova(mmod_nblarv.1ter,type="III")
testDispersion(mmod_nblarv.1ter) #borderline overdispersion
plotResiduals(mmod_nblarv.1ter)
AIC(mmod_nblarv.1ter)

#same analysis with negative binomial
mmod_nblarv.1qat<-glmmTMB(Total~(Active_substance+Dose+Clone)^2+
                            (1|Repetition),
                          data=temp,REML=FALSE,
                          family=nbinom1) #similar to a quasipoisson
summary(mmod_nblarv.1qat)
Anova(mmod_nblarv.1qat,type="III")
testDispersion(mmod_nblarv.1qat) #no overdispersion
plotResiduals(mmod_nblarv.1qat)
AIC(mmod_nblarv.1qat)

plot(mmod_nblarv.2qat) #doesn't work
#response variable vs the fitted value
plot(mmod_nblarv.2qat,Total~fitted(.)) #doesn't work
#normal distribution of errors in the different repetition
qqnorm(mmod_nblarv.2qat,~resid(.)|Repetition) #doesn't work

emmAS<-emmeans(mmod_nblarv.1qat,~Active_substance)
pairs(emmAS)
afex_plot(mmod_nblarv.1qat,"Active_substance")
summary(as.glht(pairs(emmAS)),test=adjusted("BH"))
emmDos<-emmeans(mmod_nblarv.1qat,~Dose)
pairs(emmDos)
afex_plot(mmod_nblarv.1qat,"Dose")
summary(as.glht(pairs(emmDos)),test=adjusted("BH"))
emmClone<-emmeans(mmod_nblarv.1qat,~Clone)
pairs(emmClone)
afex_plot(mmod_nblarv.1qat,"Clone")
summary(as.glht(pairs(emmClone)),test=adjusted("BH"))

#investigating the effect of Dose by active substance
emmDoByAS<-emmeans(mmod_nblarv.1qat,"Dose",by="Active_substance")
emmDoByAS
pairs(emmDoByAS)

afex_plot(mmod_nblarv.1qat,"Dose","Active_substance","Clone")
afex_plot(mmod_nblarv.1qat,"Dose",panel=~Clone+Active_substance)
#no strong clone effect
afex_plot(mmod_nblarv.1qat,"Dose",panel=~Active_substance)

#figures for the effect of AS on the mean number of larvae
interaction.plot(temp$Dose,temp$Active_substance,temp$Total,las=1)
boxplot(Total~Clone+Dose+Active_substance,data=temp,las=1,
        at=c(((1+6*0):(6*1))+2*0,
             ((1+6*1):(6*2))+2*1,
             ((1+6*2):(6*3))+2*2,
             ((1+6*3):(6*4))+2*3,
             ((1+6*4):(6*5))+2*4,
             ((1+6*5):(6*6))+2*5),
        col=rep(c(1,2),times=18),
        border=rep(c(3,4),each=18),
        names=rep(c("NT","NT","N","N","2N","2N"),times=6))
text(c((0:5)*6+(1:6)*2+1.5),y=99,
     labels=levels(temp$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))
#since there is no significant effect of the clone, they can be pooled 
#together
boxplot(Total~Dose+Active_substance,data=temp,las=1,
        at=c(((1+3*0):(3*1))+2*0,
             ((1+3*1):(3*2))+2*1,
             ((1+3*2):(3*3))+2*2,
             ((1+3*3):(3*4))+2*3,
             ((1+3*4):(3*5))+2*4,
             ((1+3*5):(3*6))+2*5),
        border=rep(c(3,4),each=9),
        names=rep(c("NT","N","2N"),times=6))
text(c((0:5)*3+(1:6)*2),y=99,
     labels=levels(temp$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))


#analysis of the effect of SA on the mean number of larvae per mother
temp<-MyzHerbi[MyzHerbi$Dose!="N/2",]
temp<-drop.levels(temp,reorder=FALSE)
mod_nbparMo<-aov(Total/Laying_females~Active_substance*Dose*Clone
                 +Error(Repetition/Dose/Clone),data=temp)
summary(mod_nbparMo)
mod_nbparMo<-aov(Total/Laying_females~(Active_substance+Dose+Clone)^2
                 +Error(Repetition/Dose/Clone),data=temp)
summary(mod_nbparMo)

#figures for the effect of SA on the mean number of larvae per mother
interaction.plot(temp$Dose,temp$Active_substance,
                 temp$Total/temp$Laying_females,las=1)
boxplot(Total/Laying_females~Clone+Dose+Active_substance,
        data=temp,las=1,
        at=c(((1+6*0):(6*1))+2*0,
             ((1+6*1):(6*2))+2*1,
             ((1+6*2):(6*3))+2*2,
             ((1+6*3):(6*4))+2*3,
             ((1+6*4):(6*5))+2*4,
             ((1+6*5):(6*6))+2*5),
        col=rep(c(1,2),times=18),
        border=rep(c(3,4),each=18),
        names=rep(c("NT","NT","N","N","2N","2N"),times=6))
text(c((0:5)*6+(1:6)*2+1.5),y=5.8,
     labels=levels(temp$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))
#since there is no significant effect of the clone, they can be pooled 
#together
boxplot(Total/Laying_females~Dose+Active_substance,data=temp,las=1,
        at=c(((1+3*0):(3*1))+2*0,
             ((1+3*1):(3*2))+2*1,
             ((1+3*2):(3*3))+2*2,
             ((1+3*3):(3*4))+2*3,
             ((1+3*4):(3*5))+2*4,
             ((1+3*5):(3*6))+2*5),
        border=rep(c(3,4),each=9),
        names=rep(c("NT","N","2N"),times=6))
text(c((0:5)*3+(1:6)*2),y=5.8,
     labels=levels(temp$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))



##############################################################################/
#Effect of the herbicides AS on survival####
##############################################################################/

#dose response analyses by replicate
aggregate(cbind(Live,Total_death,Total)~Clone+Active_substance+
            Dose,data=MyzHerbi,"sum")

MortModl<-glm(cbind(Live,Total_death)~Active_substance*Dose*Clone,
              binomial,data=temp)
summary(MortModl)
anova(MortModl)

r<-temp$Repetition
ra<-temp$Repetition:temp$Active_substance
rad<-temp$Repetition:temp$Active_substance:temp$Dose

mmod_Death.2bis<-glmmTMB(cbind(Live,Total_death)~(Active_substance+Dose+Clone)^2+
                        (1|r)+(1|ra)+(1|rad),
                      data=temp,REML=FALSE,family=binomial)

#figures for the effect of SA on the death rate of clones
interaction.plot(temp$Dose,temp$Active_substance,
                 temp$Total/temp$Laying_females,las=1)
boxplot(Total_death/Total~Clone+Dose+Active_substance,
        data=temp,las=1,
        at=c(((1+6*0):(6*1))+2*0,
             ((1+6*1):(6*2))+2*1,
             ((1+6*2):(6*3))+2*2,
             ((1+6*3):(6*4))+2*3,
             ((1+6*4):(6*5))+2*4,
             ((1+6*5):(6*6))+2*5),
        col=rep(c(1,2),times=18),
        border=rep(c(3,4),each=18),
        names=rep(c("NT","NT","N","N","2N","2N"),times=6))
text(c((0:5)*6+(1:6)*2+1.5),y=1.02,
     labels=levels(temp$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))
#figures for the effect of SA on the number of death for clone (in order
#to remove the illusion of high death rate computed on a very low number 
#of individuals)
boxplot(Total_death~Clone+Dose+Active_substance,
        data=temp,las=1,
        at=c(((1+6*0):(6*1))+2*0,
             ((1+6*1):(6*2))+2*1,
             ((1+6*2):(6*3))+2*2,
             ((1+6*3):(6*4))+2*3,
             ((1+6*4):(6*5))+2*4,
             ((1+6*5):(6*6))+2*5),
        col=rep(c(1,2),times=18),
        border=rep(c(3,4),each=18),
        names=rep(c("NT","NT","N","N","2N","2N"),times=6))
text(c((0:5)*6+(1:6)*2+1.5),y=65,
     labels=levels(temp$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))


##############################################################################/
#END
##############################################################################/