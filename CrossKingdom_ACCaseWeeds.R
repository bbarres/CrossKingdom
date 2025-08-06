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
skim(WeedInsc)
summary(WeedInsc) #the design is not complete for the Dose N/2
dfSummary(WeedInsc) #the N/2 dose has not been 

#we reduced the different phenotype categories to 2 categories:
#dead and alive individuals
WeedInsc$dead<-WeedInsc$r..1+WeedInsc$r...+WeedInsc$S
WeedInsc$alive<-WeedInsc$R+WeedInsc$r.+WeedInsc$r
str(WeedInsc)
head(WeedInsc)
skim(WeedInsc)


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
WeedSurv.mmod<-glmer(cbind(dead,alive)~(Active_substance+Dose+Phenotype+Weed)^2+
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
vif(WeedSurv.mmod) #should be <2.2 (sqrt(5)) or at least <3.2 (sqrt(10)) Doesn't work
#because no main effect and some multicollinearity, we remove the Clone 
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



##############################################################################/
#END
##############################################################################/





