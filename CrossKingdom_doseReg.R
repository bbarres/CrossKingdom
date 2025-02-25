##############################################################################/
##############################################################################/
#Regression of dose response experiment
##############################################################################/
##############################################################################/

#loading the data and packages
source("CrossKingdom_load.R")


##############################################################################/
#Dose response analyses of clone 24-0116-0001####
##############################################################################/

#dose response analyses by replicate
Clo116.m1<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             curveid=dat_test,
             weights=(nb_mtot+nb_vi),
             data=dataReg116,
             fct=LN.3u(),type="binomial")
#dose-response plot for the 24-0116-0001 replicate
plot(Clo116.m1,ylim=c(0,1.1),xlim=c(0,200),
     main="Clone 24-0116-0001 spirotetramat by replicate
     (72h after removing female)",
     ylab="mortality rate",col=TRUE,
     legendPos=c(0.5,1.1))
#LD50 comparison between replicate
compParm(Clo116.m1,"e")
EDcomp(Clo116.m1,c(50,50))
#LD50 for each replicate
ED50m1<-ED(Clo116.m1,0.50,type="absolute")

#dose response analysis with pooled replicates
Clo116.m2<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataRegRep116,
             fct=LN.3u(),type="binomial")
#LD50 of the 24-0116-0001 clone
ED50_116<-ED(Clo116.m2,0.50,type="absolute")
#dose response plot
plot(Clo116.m2,ylim=c(0,1.1),xlim=c(0,200),
     main="Clone 24-0116-0001 spirotetramat",
     type="confidence",ylab="mortality rate",
     las=1)
#adding the observed values for each dose
plot(Clo116.m2,ylim=c(0,1.1),xlim=c(0,200),
     main=names(table(dataReg116$ech_id))[1],
     type="all",add=TRUE)
#adding the average mortality rate for each dose
points((nb_mtot/(nb_mtot+nb_vi))~dose,data=dataRegDos116,
       col="red",pch=19)
#adding the LD50 value
text(x=20,y=0.3,labels=paste("ED50=",round(ED50_116[1],digits=3)))


##############################################################################/
#Dose response analyses of clone 17-0153-0020####
##############################################################################/

#dose response analyses by replicate
Clo153.m1<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataReg153,
             curveid=dat_test,
             fct=LN.3u(),type="binomial")
#dose-response plot for the 17-0116-0001 replicate
plot(Clo153.m1,ylim=c(0,1.1),xlim=c(0,100),
     main="Clone 17-0153-0020 spirotetramat by replicate
     (72h after removing female)",
     ylab="mortality rate",legendPos=c(100,0.35),col = T )
#LD50 comparison between replicate
compParm(Clo153.m1,"e")
EDcomp(Clo153.m1,c(50,50))
#LD50 for each replicate
ED50m1<-ED(Clo153.m1,0.50,type="absolute")

#dose response analysis with pooled replicates
Clo153.m2<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataRegRep153,
             fct=LN.3u(),type="binomial")
#LD50 of the 17-0153-0020 clone
ED50v<-ED(Clo153.m2,0.50,type="absolute")
#dose response plot
plot(Clo153.m2,ylim=c(0,1.1),xlim=c(0,100),
     main="Clone 17-0153-0020 spirotetramat",
     type="confidence",ylab="mortality rate",
     las=1)
#adding the observed values for each dose
plot(Clo153.m2,ylim=c(0,1.1),xlim=c(0,100),
     main=names(table(dataReg153$ech_id))[1],
     type="all",add=TRUE)
#adding the average mortality rate for each dose
points((nb_mtot/(nb_mtot+nb_vi))~dose,data=dataRegDos153,
       col="red",pch=19)
#adding the LD50 value
text(x=20,y=0.3,labels=paste("ED50=",round(ED50v[1],digits=3)))


##############################################################################/
#Comparison between sensitive and resistant clone####
##############################################################################/

#combining the dataset 
data_rep<-rbind(cbind(dataRegRep116,clone="17-0116"),
                cbind(dataRegRep153,clone="17-0153"))

Compclo<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             curveid=clone,
             weights=(nb_mtot+nb_vi),
             data=data_rep,
             fct=LN.3u(),type="binomial")

#plotting the dose-response
plot(Compclo,main="spirotetramat dose response curves 
for both clones",
     ylab="mortality rate",legendPos=c(0.05,0.8))


##############################################################################/
#END
##############################################################################/