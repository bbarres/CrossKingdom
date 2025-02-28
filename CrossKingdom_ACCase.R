##############################################################################/
##############################################################################/
#ACCase : Test of the effect of several herbicides on insect
##############################################################################/
##############################################################################/

#loading the data and packages
source("CrossKingdom_load.R")


##############################################################################/
#Effect on fertility####
##############################################################################/

#first we compare the repetitions

#dose response analyses by replicate
aggregate(cbind(Live,Total_death,Total)~Clone+Active_substance+
            Dose,data=MyzHerbi,"sum")


glm(c(V,M)~SA*Clone,binomial)
glm(c(Pond,PasPond)~SA*Clone,binomial)





#analysis of the effect of SA on the mean number of larvae
temp<-MyzHerbi[MyzHerbi$Dose!="N/2",]
temp<-drop.levels(temp,reorder=FALSE)
mod_nblarv<-aov(Total~Active_substance*Dose*Clone
                +Error(Repetition/Dose/Clone),data=temp)
summary(mod_nblarv)
mod_nblarv<-aov(Total~(Active_substance+Dose+Clone)^2
                +Error(Repetition/Dose/Clone),data=temp)
summary(mod_nblarv)
#same kind of analysis but in a glmm framework
mmod_nblarv<-lme(Total~Active_substance*Dose*Clone,
                 random= ~1|Repetition/Dose/Clone,data=temp)
summary(mmod_nblarv)
mmod_nblarv<-lme(Total~(Active_substance+Dose+Clone)^2,
                 random= ~1|Repetition/Dose/Clone,data=temp)
summary(mmod_nblarv)

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
#since there is no significant effect of the clone, they can be pooled together
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
#since there is no significant effect of the clone, they can be pooled together
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
#END
##############################################################################/