##############################################################################/
##############################################################################/
#ACCase : Test of the effect of two insecticides on several weed species
##############################################################################/
##############################################################################/

#loading the data and packages
source("CrossKingdom_load.R")

#first we inspect the data
str(WeedInsc)
head(WeedInsc)
tail(WeedInsc)
descr(WeedInsc)
skim(WeedInsc) #the design is not complete for all species
summary(WeedInsc) #there are only Resistant phenotype for 2 species
dfSummary(WeedInsc) 

#we reduced the different phenotype categories to 2 categories:
#dead and alive individuals
WeedInsc$dead<-WeedInsc$r..1+WeedInsc$r...+WeedInsc$S
WeedInsc$alive<-WeedInsc$R+WeedInsc$r.+WeedInsc$r
str(WeedInsc)
head(WeedInsc)
skim(WeedInsc)

#we first investigate the effect at dose N only 
WeedInscN<-WeedInsc[WeedInsc$Dose=="N"|
                      WeedInsc$Dose=="NT",]
WeedInscN<-drop.levels(WeedInscN,reorder=FALSE)
str(WeedInscN)
summary(WeedInscN)
skim(WeedInscN)


##############################################################################/
#Effect of the insecticides AS on survival by plant species####
##############################################################################/

#Alopecurus####
WeedAlo<-WeedInscN[WeedInscN$Weed=="Alopecurus",]
WeedAlo<-drop.levels(WeedAlo,reorder=FALSE)
str(WeedAlo)
summary(WeedAlo)
skim(WeedAlo) #there is only susceptible phenotypes for this species
boxplot(dead/Total~Dose+Active_substance,
        data=WeedAlo,las=1,ylim=c(0,1))
#simple generalized linear model
MortAlo<-glm(cbind(alive,dead)~Active_substance*Dose,
              binomial,data=WeedAlo)
summary(MortAlo)
anova(MortAlo) #no interaction effect
MortAlo<-glm(cbind(alive,dead)~Active_substance+Dose,
             binomial,data=WeedAlo)
summary(MortAlo)
anova(MortAlo)
#mixed generalized linear model
mmod_Alo<-glmer(cbind(alive,dead)~Active_substance*Dose+
                  (1|Repetition),
                  data=WeedAlo,family=binomial)
Anova(mmod_Alo,type="III") #no interaction effect
mmod_Alo<-glmer(cbind(alive,dead)~Active_substance+Dose+
                  (1|Repetition),
                data=WeedAlo,family=binomial)
Anova(mmod_Alo,type="III") #no effect of both
#checking for multicollinearity
vif(mmod_Alo) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10)) 
testDispersion(mmod_Alo) #no overdispersion
plotResiduals(mmod_Alo)
resid<-simulateResiduals(fittedModel=mmod_Alo)
plot(resid) #no dramatic deviation of the residuals
AIC(mmod_Alo)
plot(mmod_Alo)
plot(mmod_Alo,Total~fitted(.))
Anova(mmod_Alo,type="III") #no Dose main effect but strong interaction
summary(mmod_Alo)
memm<-emmeans(mmod_Alo,~(Dose+Active_substance)^2)
plot(memm)
#all simple main effect comparisons
pairs(memm,simple="each")

#Digitaria####
WeedDig<-WeedInscN[WeedInscN$Weed=="Digitaria",]
WeedDig<-drop.levels(WeedDig,reorder=FALSE)
str(WeedDig)
summary(WeedDig)
skim(WeedDig) #there is only susceptible phenotypes for this species
boxplot(dead/Total~Dose+Active_substance,
        data=WeedDig,las=1,ylim=c(0,1))
#simple generalized linear model
MortDig<-glm(cbind(alive,dead)~Active_substance*Dose,
             binomial,data=WeedDig)
summary(MortDig)
anova(MortDig) #no interaction effect
MortDig<-glm(cbind(alive,dead)~Active_substance+Dose,
             binomial,data=WeedDig)
summary(MortDig)
anova(MortDig)
#mixed generalized linear model
mmod_Dig<-glmer(cbind(alive,dead)~Active_substance*Dose+
                  (1|Repetition),
                data=WeedDig,family=binomial)
Anova(mmod_Dig,type="III") #no interaction effect
mmod_Dig<-glmer(cbind(alive,dead)~Active_substance+Dose+
                  (1|Repetition),
                data=WeedDig,family=binomial)
Anova(mmod_Dig,type="III") #no effect of both
#checking for multicollinearity
vif(mmod_Dig) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10)) 
testDispersion(mmod_Dig) #no overdispersion
plotResiduals(mmod_Dig)
resid<-simulateResiduals(fittedModel=mmod_Dig)
plot(resid) #no dramatic deviation of the residuals
AIC(mmod_Dig)
plot(mmod_Dig)
plot(mmod_Dig,Total~fitted(.))
Anova(mmod_Dig,type="III") #no Dose main effect but strong interaction
summary(mmod_Dig)
memm<-emmeans(mmod_Dig,~Dose+Active_substance)
plot(memm)
#all simple main effect comparisons
pairs(memm,simple="each")

#Lolium####
WeedLol<-WeedInscN[WeedInscN$Weed=="Lolium",]
WeedLol<-drop.levels(WeedLol,reorder=FALSE)
str(WeedLol)
summary(WeedLol)
skim(WeedLol) #there is only susceptible phenotypes for this species
boxplot(dead/Total~Dose+Active_substance,
        data=WeedLol,las=1,ylim=c(0,1)) #nothing happening
#simple generalized linear model
MortLol<-glm(cbind(alive,dead)~Active_substance*Dose,
             binomial,data=WeedLol)
summary(MortLol)
anova(MortLol) #no interaction effect
MortLol<-glm(cbind(alive,dead)~Active_substance+Dose,
             binomial,data=WeedLol)
summary(MortLol)
anova(MortLol)


#Setaria####
WeedSet<-WeedInscN[WeedInscN$Weed=="Setaria",]
WeedSet<-drop.levels(WeedSet,reorder=FALSE)
str(WeedSet)
summary(WeedSet)
skim(WeedSet) #there is both susceptible and resistant phenotypes
boxplot(dead/Total~Phenotype+Dose+Active_substance,
        data=WeedSet,las=1,ylim=c(0,1))
#simple generalized linear model
MortSet<-glm(cbind(alive,dead)~Active_substance*Dose*Phenotype,
             binomial,data=WeedSet)
summary(MortSet)
anova(MortSet) #no three way interaction effect
MortSet<-glm(cbind(alive,dead)~(Active_substance+Dose+Phenotype)^2,
             binomial,data=WeedSet)
summary(MortSet)
anova(MortSet)
#mixed generalized linear model
mmod_Set<-glmer(cbind(alive,dead)~Active_substance*Dose*Phenotype+
                  (1|Repetition),
                data=WeedSet,family=binomial)
Anova(mmod_Set,type="III") #no three way interaction effect
mmod_Set<-glmer(cbind(alive,dead)~(Active_substance+Dose+Phenotype)^2+
                  (1|Repetition),
                data=WeedSet,family=binomial)
Anova(mmod_Set,type="III") #no effect of both
mmod_Set<-glmer(cbind(alive,dead)~(Active_substance+Dose+Phenotype
                                   +Active_substance:Dose+
                                     Active_substance:Phenotype)+
                  (1|Repetition),
                data=WeedSet,family=binomial)
Anova(mmod_Set,type="III") #no effect of both
mmod_Set<-glmer(cbind(alive,dead)~(Active_substance+Dose+Phenotype
                                   +Active_substance:Phenotype)+
                  (1|Repetition),
                data=WeedSet,family=binomial)
Anova(mmod_Set,type="III") #no effect of both
mmod_Set<-glmer(cbind(alive,dead)~(Active_substance+Dose+Phenotype)+
                  (1|Repetition),
                data=WeedSet,family=binomial)
Anova(mmod_Set,type="III") #no effect of both
#checking for multicollinearity
vif(mmod_Set) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10)) 
testDispersion(mmod_Set) #no overdispersion
plotResiduals(mmod_Set)
resid<-simulateResiduals(fittedModel=mmod_Set)
plot(resid) #no dramatic deviation of the residuals
AIC(mmod_Set)
plot(mmod_Set)
plot(mmod_Set,Total~fitted(.))
Anova(mmod_Set,type="III") #no Dose main effect but strong interaction
summary(mmod_Set)
memm<-emmeans(mmod_Set,~Dose+Active_substance)
plot(memm)
#all simple main effect comparisons
pairs(memm,simple="each")

#Zea####
WeedZea<-WeedInscN[WeedInscN$Weed=="Zea",]
WeedZea<-drop.levels(WeedZea,reorder=FALSE)
str(WeedZea)
summary(WeedZea)
skim(WeedZea) #there is both susceptible and resistant phenotypes
boxplot(dead/Total~Phenotype+Dose+Active_substance,
        data=WeedZea,las=1,ylim=c(0,1))
#simple generalized linear model
MortZea<-glm(cbind(alive,dead)~Active_substance*Dose*Phenotype,
             binomial,data=WeedZea)
summary(MortZea)
anova(MortZea) #no three way interaction effect
MortZea<-glm(cbind(alive,dead)~(Active_substance+Dose+Phenotype)^2,
             binomial,data=WeedZea)
summary(MortZea)
anova(MortZea)
#mixed generalized linear model
mmod_Zea<-glmer(cbind(alive,dead)~Active_substance*Dose*Phenotype+
                  (1|Repetition),
                data=WeedZea,family=binomial)
Anova(mmod_Zea,type="III") #no three way interaction effect
mmod_Zea<-glmer(cbind(alive,dead)~(Active_substance+Dose+Phenotype)^2+
                  (1|Repetition),
                data=WeedZea,family=binomial)
Anova(mmod_Zea,type="III") #no effect of both
mmod_Zea<-glmer(cbind(alive,dead)~(Active_substance+Dose+Phenotype
                                   +Active_substance:Dose+
                                     Active_substance:Phenotype)+
                  (1|Repetition),
                data=WeedZea,family=binomial)
Anova(mmod_Zea,type="III") #no effect of both
mmod_Zea<-glmer(cbind(alive,dead)~(Active_substance+Dose+Phenotype
                                   +Active_substance:Phenotype)+
                  (1|Repetition),
                data=WeedZea,family=binomial)
Anova(mmod_Zea,type="III") #no effect of both
mmod_Zea<-glmer(cbind(alive,dead)~(Active_substance+Dose+Phenotype)+
                  (1|Repetition),
                data=WeedZea,family=binomial)
Anova(mmod_Zea,type="III") #no effect of both
#checking for multicollinearity
vif(mmod_Zea) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10)) 
testDispersion(mmod_Zea) #no overdispersion
plotResiduals(mmod_Zea)
resid<-simulateResiduals(fittedModel=mmod_Zea)
plot(resid) #no dramatic deviation of the residuals
AIC(mmod_Zea)
plot(mmod_Zea)
plot(mmod_Zea,Total~fitted(.))
Anova(mmod_Zea,type="III") #no Dose main effect but strong interaction
summary(mmod_Zea)
memm<-emmeans(mmod_Zea,~Dose+Active_substance)
plot(memm)
#all simple main effect comparisons
pairs(memm,simple="each")




##############################################################################/
#Estimating DL50 of the two insecticide on Setaria####
##############################################################################/

#preparing the dataset
SetDosRep<-WeedInsc[WeedInsc$Weed=="Setaria",]
SetDosRep$DoseQ<-SetDosRep$Dose

#dose response analyses by population for spiromesifen
SetDosSpiM<-SetDosRep[SetDosRep$Active_substance=="Spiromesifen",]
levels(SetDosSpiM$DoseQ)<-c("0","0.50","1","2") #g/L or 150-300-600g/ha
SetDosSpiM$DoseQ<-as.numeric(as.character(SetDosSpiM$DoseQ))
SpiromesifenS.m1<-drm(dead/Total~DoseQ,
                  curveid=Phenotype,
                  weights=Total,
                  data=SetDosSpiM,
                  fct=LN.3u(),type="binomial")
plot(SpiromesifenS.m1,ylim=c(0,1.1),xlim=c(0,10),
     main="Setaria - Spiromesifen",
     ylab="mortality rate",col=TRUE,
     legendPos=c(0.5,1.1))
#estimating EC50
ED50SpimS<-ED(SpiromesifenS.m1,0.50,type="absolute")
#FR for the "resistant" population
ED50SpimS[2]/ED50SpimS[1]

#dose response analyses by population for spirotetramat
SetDosSpiT<-SetDosRep[SetDosRep$Active_substance=="Spirotetramat",]
levels(SetDosSpiT$DoseQ)<-c("0","0.50","1","2") #g/L or 150-300-600g/ha
SetDosSpiT$DoseQ<-as.numeric(as.character(SetDosSpiT$DoseQ))
SpirotetramatS.m1<-drm(dead/Total~DoseQ,
                     curveid=Phenotype,
                     weights=Total,
                     data=SetDosSpiT,
                     fct=LN.3u(),type="binomial")
plot(SpirotetramatS.m1,ylim=c(0,1.1),xlim=c(0,5),
     main="Setaria - Spirotetramat",
     ylab="mortality rate",col=TRUE,
     legendPos=c(0.5,1.1))
#estimating EC50
ED50SpitS<-ED(SpirotetramatS.m1,0.50,type="absolute")
#FR for the "resistant" population
ED50SpitS[2]/ED50SpitS[1]


##############################################################################/
#Estimating DL50 of the two insecticide on Zea####
##############################################################################/

#preparing the dataset
ZeaDosRep<-WeedInsc[WeedInsc$Weed=="Zea",]
ZeaDosRep$DoseQ<-ZeaDosRep$Dose

#dose response analyses by population for spirotetramat
ZeaDosSpiT<-ZeaDosRep[ZeaDosRep$Active_substance=="Spirotetramat",]
levels(ZeaDosSpiT$DoseQ)<-c("0","0.50","1","2") #g/L or 150-300-600g/ha
ZeaDosSpiT$DoseQ<-as.numeric(as.character(ZeaDosSpiT$DoseQ))
SpirotetramatZ.m1<-drm(dead/Total~DoseQ,
                      curveid=Phenotype,
                      weights=Total,
                      data=ZeaDosSpiT,
                      fct=LN.3u(),type="binomial")
plot(SpirotetramatZ.m1,ylim=c(0,1.1),xlim=c(0,5),
     main="Zea - Spirotetramat",
     ylab="mortality rate",col=TRUE,
     legendPos=c(0.5,1.1))
#estimating EC50
ED50SpitZ<-ED(SpirotetramatZ.m1,0.50,type="absolute")
#FR for the "resistant" population
ED50SpitZ[2]/ED50SpitZ[1]


##############################################################################/
#END
##############################################################################/











#we first investigate the effect at dose N only 
WeedInscN<-WeedInsc[WeedInsc$Dose=="N"|
                      WeedInsc$Dose=="NT",]
WeedInscN<-drop.levels(WeedInscN,reorder=FALSE)
str(WeedInscN)
summary(WeedInscN)
skim(WeedInscN)
#and on sensitive populations only
WeedInscNS<-WeedInscN[WeedInscN$Phenotype=="S",]
WeedInscNS<-drop.levels(WeedInscNS,reorder=FALSE)
str(WeedInscNS)
summary(WeedInscNS)
skim(WeedInscNS)

#Mixed effect model with repetition as random effect
WeedSen.mmod<-glmer(cbind(dead,alive)~Active_substance*Dose*Weed+
                       (1|Repetition),
                     data=WeedInscNS,family=binomial)
Anova(WeedSen.mmod,type="III") #no three way interaction effect
WeedSen.mmod<-glmer(cbind(dead,alive)~(Active_substance+Dose+Weed)^2+
                       (1|Repetition),
                     data=WeedInscNS,family=binomial)
Anova(WeedSen.mmod,type="III") #removing Dose:Weed interaction
WeedSen.mmod<-glmer(cbind(dead,alive)~Active_substance+Dose+Weed+
                       Active_substance:Dose+Active_substance:Weed+
                       (1|Repetition),
                     data=WeedInscNS,family=binomial)
Anova(WeedSen.mmod,type="III") #removing Active_substance:Weed interaction
WeedSen.mmod<-glmer(cbind(dead,alive)~Active_substance+Dose+Weed+
                      Active_substance:Dose+
                      (1|Repetition),
                     data=WeedInscNS,family=binomial)
Anova(WeedSen.mmod,type="III") #removing Active_substance:Dose interaction
WeedSen.mmod<-glmer(cbind(dead,alive)~Active_substance+Dose+Weed+
                       (1|Repetition),
                     data=WeedInscNS,family=binomial)
Anova(WeedSen.mmod,type="III") #removing Phenotype:Weed interaction
WeedSen.mmod<-glmer(cbind(dead,alive)~Active_substance+Dose+Weed+
                       Active_substance:Dose+
                       (1|Repetition),
                     data=WeedInscNS,family=binomial)
Anova(WeedSen.mmod,type="III") #final model


##############################################################################/
#Effect of the insecticides AS on survival####
##############################################################################/

#we first investigate the effect at dose N only
WeedInscN<-WeedInsc[WeedInsc$Dose=="N"|
                      WeedInsc$Dose=="NT",]
WeedInscN<-drop.levels(WeedInscN,reorder=FALSE)
str(WeedInscN)
summary(WeedInscN)
skim(WeedInscN)


WeedSurv.mod<-glm(cbind(dead,alive)~Active_substance*Dose*Phenotype*Weed,
              binomial,data=WeedInscN)
summary(WeedSurv.mod)
anova(WeedSurv.mod)
#we can remove the four way and three way interactions
WeedSurv.mod<-glm(cbind(dead,alive)~(Active_substance+Dose+Phenotype+Weed)^2,
                  binomial,data=WeedInscN)
summary(WeedSurv.mod)
anova(WeedSurv.mod) #removing Active_substance:Weed interaction
WeedSurv.mod<-glm(cbind(dead,alive)~Active_substance+Dose+Phenotype+Weed+
                    Active_substance:Dose+Active_substance:Phenotype+
                    Dose:Phenotype+Dose:Weed+Phenotype:Weed,
                  binomial,data=WeedInscN)
summary(WeedSurv.mod)
anova(WeedSurv.mod) #removing Dose:Phenotype interaction
WeedSurv.mod<-glm(cbind(dead,alive)~Active_substance+Dose+Phenotype+Weed+
                    Active_substance:Dose+Active_substance:Phenotype+
                    Dose:Weed+Phenotype:Weed,
                  binomial,data=WeedInscN)
summary(WeedSurv.mod)
anova(WeedSurv.mod) #removing Phenotype:Weed interaction
WeedSurv.mod<-glm(cbind(dead,alive)~Active_substance+Dose+Phenotype+Weed+
                    Active_substance:Dose+Active_substance:Phenotype+
                    Dose:Weed,
                  binomial,data=WeedInscN)
summary(WeedSurv.mod)
anova(WeedSurv.mod) #removing Phenotype:Weed interaction


#with repetition as a random effect
WeedSurv.mmod<-glmmTMB(cbind(dead,alive)~(Active_substance+Dose+Phenotype)^2+
                      (1|Repetition),
                    data=WeedInscN,REML=FALSE,family=binomial)
WeedSurv.mmod<-glmer(cbind(dead,alive)~Active_substance*Dose*Phenotype*Weed+
                    (1|Repetition),
                  data=WeedInscN,family=binomial)
Anova(WeedSurv.mmod,type="III") #no three way interaction effect
WeedSurv.mmod<-glmer(cbind(dead,alive)~(Active_substance+Dose+
                                          Phenotype+Weed)^2+
                       (1|Repetition),
                     data=WeedInscN,family=binomial)
Anova(WeedSurv.mmod,type="III") #removing Active_substance:Weed interaction
WeedSurv.mmod<-glmer(cbind(dead,alive)~Active_substance+Dose+Phenotype+Weed+
                       Active_substance:Dose+Active_substance:Phenotype+
                       Dose:Phenotype+Dose:Weed+Phenotype:Weed+
                       (1|Repetition),
                     data=WeedInscN,family=binomial)
Anova(WeedSurv.mmod,type="III") #removing Dose:Weed interaction
WeedSurv.mmod<-glmer(cbind(dead,alive)~Active_substance+Dose+Phenotype+Weed+
                       Active_substance:Dose+Active_substance:Phenotype+
                       Dose:Phenotype+Phenotype:Weed+
                       (1|Repetition),
                     data=WeedInscN,family=binomial)
Anova(WeedSurv.mmod,type="III") #removing Dose:Phenotype interaction
WeedSurv.mmod<-glmer(cbind(dead,alive)~Active_substance+Dose+Phenotype+Weed+
                       Active_substance:Dose+Active_substance:Phenotype+
                       Phenotype:Weed+
                       (1|Repetition),
                     data=WeedInscN,family=binomial)
Anova(WeedSurv.mmod,type="III") #removing Phenotype:Weed interaction
WeedSurv.mmod<-glmer(cbind(dead,alive)~Active_substance+Dose+Phenotype+Weed+
                       Active_substance:Dose+Active_substance:Phenotype+
                       (1|Repetition),
                     data=WeedInscN,family=binomial)
Anova(WeedSurv.mmod,type="III") #final model

#checking for multicollinearity
vif(WeedSurv.mmod) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10)) 
                   #Doesn't work
#because no main effect and some multicollinearity, we remove the Clone effect
testDispersion(WeedSurv.mmod) #some overdispersion
plotResiduals(WeedSurv.mmod)
resid<-simulateResiduals(fittedModel=WeedSurv.mmod)
plot(resid) #not dramatic but some problems
AIC(WeedSurv.mmod)
plot(WeedSurv.mmod)
plot(WeedSurv.mmod,Total~fitted(.))
qqnorm(WeedSurv.mmod,~resid(.)|Repetition)
Anova(WeedSurv.mmod,type="III") #no Dose main effect but strong interaction
summary(WeedSurv.mmod)
Wemm<-emmeans(WeedSurv.mmod,~Dose+Active_substance+Weed)
plot(Wemm)
#all simple main effect comparisons
pairs(Wemm,simple="each")


#figures for the effect of SA on the death rate of clones
interaction.plot(WeedInscN$Dose,WeedInscN$Active_substance,
                 WeedInscN$dead/WeedInscN$Total,las=1)
boxplot(dead/Total~Phenotype+Dose+Active_substance+Weed,
        data=WeedInscN,las=1,
        at=c(((1+4*0):(4*1))+2*0,
             ((1+4*1):(4*2))+2*1,
             ((1+4*2):(4*3))+2*2,
             ((1+4*3):(4*4))+2*3),
        col=rep(c(1,2),times=8),
        border=rep(c(3,4),each=8),
        names=rep(c("NT","NT","N","N"),times=4))
text(c((0:3)*4+(1:4)*2+1),y=1.02,
     labels=levels(MyzHerbiM$Active_substance),cex=1.5,
     col=c(3,3,4,4))
#figures for the effect of SA on the number of death for clone (in order
#to remove the illusion of high death rate computed on a very low number 
#of individuals)
boxplot(Total_death~Clone+Dose+Active_substance,
        data=MyzHerbiS,las=1,
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
     labels=levels(MyzHerbiS$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))

