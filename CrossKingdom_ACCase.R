##############################################################################/
##############################################################################/
#ACCase : Test of the effect of several herbicides on an aphid species
##############################################################################/
##############################################################################/

#loading the data and packages
source("CrossKingdom_load.R")


##############################################################################/
#Effect on fertility####
##############################################################################/

#analysis of the effect of SA on the mean number of larvae
temp<-MyzHerbi[MyzHerbi$Dose!="N/2",]
temp<-drop.levels(temp,reorder=FALSE)
mod_nblarv<-aov(Total~Active_substance*Dose*Clone
                +Error(Repetition/Dose/Clone),data=temp)
summary(mod_nblarv)
mod_nblarv<-aov(Total~(Active_substance+Dose+Clone)^2
                +Error(Repetition/Dose/Clone),data=temp)
summary(mod_nblarv)
TukeyHSD(mod_nblarv) #doesn't work with an Error term in the aov

emm<-emmeans(mod_nblarv,~(Active_substance+Dose+Clone)^2)
plot(emm)
emm<-emmeans(mod_nblarv,~(Clone+Dose+Active_substance)^2)
plot(emm)
#alternatively, you can use mixed model and glht function from the 
#multcomp package other possibilité afex package function mixed + 
#lsmeans package


#same kind of analysis but in a glmm framework
mmod_nblarv<-lme(Total~Active_substance*Dose*Clone,
                 random= ~1|Repetition/Dose/Clone,
                 data=temp,method="ML")
summary(mmod_nblarv)
mmod_nblarv.1<-update(mmod_nblarv,~. -Active_substance:Dose:Clone)
summary(mmod_nblarv.1)
anova(mmod_nblarv.1,mmod_nblarv)
mmod_nblarv.2<-update(mmod_nblarv.1,~. -Active_substance:Clone)
anova(mmod_nblarv.2,mmod_nblarv.1)
#the interaction Active_substance:Clone can be removed
mmod_nblarv.3<-update(mmod_nblarv.2,~. -Dose:Clone)
anova(mmod_nblarv.3,mmod_nblarv.2)
mmod_nblarv.4<-update(mmod_nblarv.2,~. -Active_substance:Dose)
anova(mmod_nblarv.4,mmod_nblarv.2)

#here is the "minimal" model
summary(mmod_nblarv.2)
#checking the residuals
plot(mmod_nblarv.2)
#response variable vs the fitted value
plot(mmod_nblarv.2,Total~fitted(.))
#normal distribution of errors in the different repetition
qqnorm(mmod_nblarv.2,~resid(.)|Repetition)
Anova(mmod_nblarv.2,type=c("III"))
#checking for overdispersion
overdisp_fun = function(model) {
  sum(residuals(model,type="pearson")^2)/df.residual(model)
}

overdisp_fun(mmod_nblarv.2)


#same model but using lmer
r<-temp$Repetition
rd<-temp$Repetition:temp$Dose
rdc<-temp$Repetition:temp$Dose:temp$Clone
mmod_nblarv.2bis<-lmer(Total~Active_substance+Dose+Clone+
                         Active_substance:Dose+Dose:Clone+(1|r)+(1|rd),
                       data=temp,REML=FALSE)
print(mmod_nblarv.2bis,cor=FALSE)

mmod_nblarv.2ter<-glmmTMB(Total~Active_substance+Dose+Clone+
                            Active_substance:Dose+Dose:Clone+(1|r)+(1|rd),
                          data=temp,REML=FALSE)
summary(mmod_nblarv.2ter)
overdisp_fun(mmod_nblarv.2ter)
Anova(mmod_nblarv.2ter)

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