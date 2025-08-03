##############################################################################/
##############################################################################/
#ACCase : Test of the effect of several herbicides on an aphid species
##############################################################################/
##############################################################################/

#loading the data and packages
source("CrossKingdom_load.R")

#first we inspect the data
str(MyzHerbi)
head(MyzHerbi)
tail(MyzHerbi)
descr(MyzHerbi)
skim(MyzHerbi)
summary(MyzHerbi) #the design is not complete for the Dose N/2
dfSummary(MyzHerbi) #the N/2 dose has not been 
#we remove the N/2 to streamline the analyses
MyzHerbiS<-MyzHerbi[MyzHerbi$Dose!="N/2",]
MyzHerbiS<-drop.levels(MyzHerbiS,reorder=FALSE)
str(MyzHerbiS)
summary(MyzHerbiS)
skim(MyzHerbiS)


##############################################################################/
#Effect on total fertility by lm####
##############################################################################/

#first we investigate the effect of the different factors on total
#fertility
hist(MyzHerbiS$Total)
#analysis of the effect of SA on the mean number of larvae
mod_nblarv<-aov(Total~Active_substance*Dose*Clone
                +Error(Repetition),data=MyzHerbiS)
summary(mod_nblarv)
#we can remove the three way interaction
mod_nblarv<-aov(Total~(Active_substance+Dose+Clone)^2
                +Error(Repetition),data=MyzHerbiS)
summary(mod_nblarv)

emm<-emmeans(mod_nblarv,~(Clone+Dose+Active_substance)^2)
plot(emm) #the magnitude of the differences between clones are small
#all simple main effect comparisons
pairs(emm,simple="each")
pairs(emm,simple="each")[1] #very few comparison between clones significant
emmip(mod_nblarv,~Dose|Active_substance*Clone)
emmip(mod_nblarv,Active_substance~Dose|Clone)
joint_tests(mod_nblarv,by="Clone")
#same model but without the Clone variable as a main effect but as a 
#random effect instead
mod_nblarv<-aov(Total~(Active_substance+Dose)^2
                +Error(Repetition+Clone),data=MyzHerbiS)
summary(mod_nblarv)
emm<-emmeans(mod_nblarv,~(Dose+Active_substance)^2)
plot(emm)
#all simple main effect comparisons
pairs(emm,simple="each")
summary(as.glht(pairs(emm,simple="Dose")),test=adjusted("BH"))
summary(as.glht(pairs(emm,simple="Active_substance")),test=adjusted("BH"))
#no significant difference between active substance at dose NT
#At dose N, only Quizalofop is significantly different from all other AS
#At dose 2N, in addition to Quizalofop, differences also for Cycloxydim
#Within AS, only Quizalofop and Cycloxydim display a significant 
#difference between doses NT and N. There are also additional marginal 
#differences for Clodinafop and Fluazifop between doses NT and 2N
#In conclusion, it seems that at dose N, Quizalofop and Cycloxydim
#affect the fertility whatever the clone  -> these AS should be removed 
#for survival analyses


##############################################################################/
#Effect on fertility using glmm####
##############################################################################/

#same model but using lmer
mmod_nblarv.1<-lmer(Total~Active_substance*Dose*Clone+
                         (1|Repetition),
                       data=MyzHerbiS,REML=FALSE)
mmod_nblarv.2<-lmer(Total~(Active_substance+Dose+Clone)^2+
                         (1|Repetition),
                       data=MyzHerbiS,REML=FALSE)
anova(mmod_nblarv.2,mmod_nblarv.1) #three way interaction can be removed

plot(mmod_nblarv.2)
plot(mmod_nblarv.2,Total~fitted(.))
Anova(mmod_nblarv.2,type="III") #no Clone main effect
#checking for multicollinearity
vif(mmod_nblarv.2) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10))
#because no main effect and some multicollinearity, we remove the Clone 
mmod_nblarv<-lmer(Total~(Active_substance+Dose)^2+
                    (1|Repetition),
                  data=MyzHerbiS,REML=FALSE)
anova(mmod_nblarv,mmod_nblarv.2) #marginal improvement of loglik
vif(mmod_nblarv) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10))
#less colinearity
testDispersion(mmod_nblarv) #no overdispersion
plotResiduals(mmod_nblarv)
resid<-simulateResiduals(fittedModel=mmod_nblarv)
plot(resid)
AIC(mmod_nblarv)
plot(mmod_nblarv)
plot(mmod_nblarv,Total~fitted(.))
qqnorm(mmod_nblarv,~resid(.)|Repetition) #not working
Anova(mmod_nblarv,type="III") #no Dose main effect but strong interaction
summary(mmod_nblarv)
memm<-emmeans(mmod_nblarv,~(Dose+Active_substance)^2)
plot(memm)
#all simple main effect comparisons
pairs(memm,simple="each")
summary(as.glht(pairs(memm,simple="Dose")),test=adjusted("BH"))
summary(as.glht(pairs(memm,simple="Active_substance")),test=adjusted("BH"))
#no significant difference between active substance at dose NT
#At dose N, only Quizalofop is significantly different from all other AS
#At dose 2N, in addition to Quizalofop, differences also for Cycloxydim
#Within AS, only Quizalofop and Cycloxydim display a significant 
#difference between doses NT and N. There are also additional marginal 
#differences for Clodinafop and Fluazifop between doses NT and 2N
#In conclusion, it seems that at dose N, Quizalofop and Cycloxydim
#affect the fertility whatever the clone  -> these AS should be removed 
#for survival analyses



#figures for the effect of AS on the mean number of larvae
interaction.plot(MyzHerbiS$Dose,MyzHerbiS$Active_substance,
                 MyzHerbiS$Total,las=1)
boxplot(Total~Clone+Dose+Active_substance,data=MyzHerbiS,las=1,
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
     labels=levels(MyzHerbiS$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))
#since there is no significant effect of the clone, they can be pooled 
#together
boxplot(Total~Dose+Active_substance,data=MyzHerbiS,las=1,
        at=c(((1+3*0):(3*1))+2*0,
             ((1+3*1):(3*2))+2*1,
             ((1+3*2):(3*3))+2*2,
             ((1+3*3):(3*4))+2*3,
             ((1+3*4):(3*5))+2*4,
             ((1+3*5):(3*6))+2*5),
        border=rep(c(3,4),each=9),
        names=rep(c("NT","N","2N"),times=6))
text(c((0:5)*3+(1:6)*2),y=99,
     labels=levels(MyzHerbiS$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))



##############################################################################/
#Effect of the herbicides AS on survival####
##############################################################################/

#removing unnecessary active substances and dose
#because of the induces and unspecific mortality we remove 
#Quizalofop and Cycloxydim
#also in a first step and because we are primarily interested by the dose N
#we remove the dose 2N to simplify the first modelling
MyzHerbiM<-MyzHerbiS[MyzHerbiS$Dose!="2N",]
MyzHerbiM<-MyzHerbiM[MyzHerbiM$Active_substance!="Quizalofop",]
MyzHerbiM<-MyzHerbiM[MyzHerbiM$Active_substance!="Cycloxydim",]
MyzHerbiM<-drop.levels(MyzHerbiM,reorder=FALSE)
str(MyzHerbiM)
summary(MyzHerbiM)
skim(MyzHerbiM)


MortModl<-glm(cbind(Live,Total_death)~Active_substance*Dose*Clone,
              binomial,data=MyzHerbiM)
summary(MortModl)
anova(MortModl)

r<-MyzHerbiM$Repetition
ra<-MyzHerbiM$Repetition:MyzHerbiM$Active_substance
rad<-MyzHerbiM$Repetition:MyzHerbiM$Active_substance:MyzHerbiM$Dose


mmod_Death<-glmer(cbind(Total_death,Live)~(Active_substance+Dose+Clone)^2+
                    (1|r),
                  data=MyzHerbiM,family=binomial)

mmod_Death<-glmmTMB(cbind(Total_death,Live)~(Active_substance+Dose+Clone)^2+
                        (1|Repetition),
                      data=MyzHerbiM,REML=FALSE,family=binomial)
mmod_Death<-glmer(cbind(Total_death,Live)~Active_substance*Dose*Clone+
                    (1|Repetition),
                  data=MyzHerbiM,family=binomial)

Anova(mmod_Death,type="III") #no Clone main effect
#checking for multicollinearity
vif(mmod_Death) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10)) Doesn't work
#because no main effect and some multicollinearity, we remove the Clone 
testDispersion(mmod_Death) #some overdispersion
plotResiduals(mmod_Death)
resid<-simulateResiduals(fittedModel=mmod_Death)
plot(resid) #not dramatic but some problems
AIC(mmod_Death)
plot(mmod_Death)
plot(mmod_Death,Total~fitted(.))
qqnorm(mmod_Death,~resid(.)|Repetition)
Anova(mmod_Death,type="III") #no Dose main effect but strong interaction
summary(mmod_Death)
memm<-emmeans(mmod_Death,~(Clone+Dose+Active_substance)^2)
plot(memm)
#all simple main effect comparisons
pairs(memm,simple="each")


#figures for the effect of SA on the death rate of clones
interaction.plot(MyzHerbiM$Dose,MyzHerbiM$Active_substance,
                 MyzHerbiM$Total_death/MyzHerbiM$Total,las=1)
boxplot(Total_death/Total~Clone+Dose+Active_substance,
        data=MyzHerbiM,las=1,
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



##############################################################################/
#Estimating DL50 of the two herbicides on Myzus####
##############################################################################/

#preparing the dataset
MyzDosRep<-MyzHerbi[MyzHerbi$Active_substance=="Clethodim"|
                      MyzHerbi$Active_substance=="Pinoxaden",]
MyzDosRep$DoseQ<-MyzDosRep$Dose

#dose response analyses by clone for clethodim
MyzDosClet<-MyzDosRep[MyzDosRep$Active_substance=="Clethodim",]
levels(MyzDosClet$DoseQ)<-c("0","1.50","3.00","6.00") #150-300-600g/ha
MyzDosClet$DoseQ<-as.numeric(as.character(MyzDosClet$DoseQ))
Clethodim.m1<-drm(Total_death/Total~DoseQ,
               curveid=Clone,
               weights=Total,
               data=MyzDosClet,
               fct=LN.3u(),type="binomial")
plot(Clethodim.m1,ylim=c(0,1.1),xlim=c(0,10),
     main="Clethodim",
     ylab="mortality rate",col=TRUE,
     legendPos=c(0.5,1.1))
#estimating EC50
ED50cle<-ED(Clethodim.m1,0.50,type="absolute")
#FR for the resistant clone
ED50cle[2]/ED50cle[1]

#dose response analyses by clone for clethodim
MyzDosPino<-MyzDosRep[MyzDosRep$Active_substance=="Pinoxaden",]
levels(MyzDosPino$DoseQ)<-c("0","3.0","6.0","12.0") #30-60-120g/ha
MyzDosPino$DoseQ<-as.numeric(as.character(MyzDosPino$DoseQ))
Pinoxaden.m1<-drm(Total_death/Total~DoseQ,
                  curveid=Clone,
                  weights=Total,
                  data=MyzDosPino,
                  fct=LN.3u(),type="binomial")
plot(Pinoxaden.m1,ylim=c(0,1.1),xlim=c(0,50),
     main="Pinoxaden",
     ylab="mortality rate",col=TRUE,
     legendPos=c(0.9,1.1))
#estimating EC50
ED50pin<-ED(Pinoxaden.m1,0.50,type="absolute")
#FR for the resistant clone
ED50pin[2]/ED50pin[1]


##############################################################################/
#END
##############################################################################/




#using the afex package
#first we need an ID for the different individual
MyzHerbiS$ID<-paste(MyzHerbiS$Clone,MyzHerbiS$Active_substance,
                    MyzHerbiS$Dose,sep="_")
mod_afex<-aov_car(Total~((Active_substance+Dose+Clone)^2+
                           Error(ID/Repetition)),
                  data=MyzHerbiS)
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




#same kind of analysis but in a glmm framework
mmod_nblarv<-lme(Total~Active_substance*Dose*Clone,
                 random= ~1|Repetition,
                 data=MyzHerbiS,method="ML")
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
                       data=MyzHerbiS,REML=FALSE)
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
                          data=MyzHerbiS,REML=FALSE,family=poisson)
summary(mmod_nblarv.1ter)
Anova(mmod_nblarv.1ter,type="III")
testDispersion(mmod_nblarv.1ter) #borderline overdispersion
plotResiduals(mmod_nblarv.1ter)
AIC(mmod_nblarv.1ter)


#same analysis with negative binomial
mmod_nblarv.1qat<-glmmTMB(Total~(Active_substance+Dose+Clone)^2+
                            (1|Repetition),
                          data=MyzHerbiS,REML=FALSE,
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






##############################################################################/
#Mean number of larvae as response variable
##############################################################################/


#analysis of the effect of SA on the mean number of larvae per mother
MyzHerbiS<-MyzHerbi[MyzHerbi$Dose!="N/2",]
MyzHerbiS<-drop.levels(MyzHerbiS,reorder=FALSE)
mod_nbparMo<-aov(Total/Laying_females~Active_substance*Dose*Clone
                 +Error(Repetition/Dose/Clone),data=MyzHerbiS)
summary(mod_nbparMo)
mod_nbparMo<-aov(Total/Laying_females~(Active_substance+Dose+Clone)^2
                 +Error(Repetition/Dose/Clone),data=MyzHerbiS)
summary(mod_nbparMo)

#figures for the effect of SA on the mean number of larvae per mother
interaction.plot(MyzHerbiS$Dose,MyzHerbiS$Active_substance,
                 MyzHerbiS$Total/MyzHerbiS$Laying_females,las=1)
boxplot(Total/Laying_females~Clone+Dose+Active_substance,
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
text(c((0:5)*6+(1:6)*2+1.5),y=5.8,
     labels=levels(MyzHerbiS$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))
#since there is no significant effect of the clone, they can be pooled 
#together
boxplot(Total/Laying_females~Dose+Active_substance,data=MyzHerbiS,las=1,
        at=c(((1+3*0):(3*1))+2*0,
             ((1+3*1):(3*2))+2*1,
             ((1+3*2):(3*3))+2*2,
             ((1+3*3):(3*4))+2*3,
             ((1+3*4):(3*5))+2*4,
             ((1+3*5):(3*6))+2*5),
        border=rep(c(3,4),each=9),
        names=rep(c("NT","N","2N"),times=6))
text(c((0:5)*3+(1:6)*2),y=5.8,
     labels=levels(MyzHerbiS$Active_substance),cex=1.5,
     col=c(3,3,3,4,4,4))



