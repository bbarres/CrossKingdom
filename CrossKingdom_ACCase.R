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
aggregate(cbind(Vivants,Morts_H72,Total)~Clone+SA+Dose,data=MyzHerbi,"sum")


glm(c(V,M)~SA*Clone,binomial)
glm(c(Pond,PasPond)~SA*Clone,binomial)





#analysis of the effect on number of larvae / mean number of larvae
temp<-MyzHerbi[MyzHerbi$Dose!="N/2",]
temp<-drop.levels(temp)
temp$Dose<-factor(temp$Dose,levels=c("NT","N","2N"))
mod_nblarv<-aov(Total~SA*Dose*Clone+Error(Repetition/Dose/Clone),data=temp)
summary(mod_nblarv)
interaction.plot(temp$Dose,temp$SA,temp$Total,las=1)


##############################################################################/
#END
##############################################################################/