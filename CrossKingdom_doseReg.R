##############################################################################/
##############################################################################/
#Loading the data sets and packages needed for CrossKingdom analyses
##############################################################################/
##############################################################################/

#loading the data and packages
source("CrossKingdom_load.R")


##############################################################################/
#Loading and preparing the main data set####
##############################################################################/




##############################################################################/
#Writing session information for reproducibility####
##############################################################################/

#dose response analyses by date
Clo116.m1<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             curveid=dat_test,
             weights=(nb_mtot+nb_vi),
             data=dataReg116,
             fct=LN.3u(),type="binomial")
#plot de la courbe dose-réponse
plot(Clo116.m1,ylim=c(0,1.1),xlim=c(0,100),
     main="Clone 24-0116-0001 spirotetramat by date 
     (72h after removing female)",
     ylab="mortality rate",col=TRUE,
     )
#comparaison des courbes entre dates
compParm(Clo116.m1,"e")
#comparaison DL50 entre dates
EDcomp(Clo116.m1,c(50,50))
#DL 50 pour chaque date
ED50m1<-ED(Clo116.m1,0.50,type="absolute")


# Création de la courbe dose réponse avec
# le data où les vivants et les morts sont agrégés par date et dose
Clo116.m2<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataRegRep,
             fct=LN.3u(),type="binomial")
#DL50 de la courbe dose réponse avec cumul de dates
ED50v<-ED(Clo116.m2,0.50,type="absolute")

#plot de la courbe "moyenne" avec intervalles de confiance
plot(Clo116.m2,ylim=c(0,1.1),xlim=c(0,100),
     main="Courbe moyenne dose-réponse du clone 24-0116-0001 \nau spirotétramate",
     type="confidence",ylab="mortality rate",
     las=1)

# plot des taux de mortalité par dose et par date
plot(Clo116.m2,ylim=c(0,1.1),xlim=c(0,100),
     main=names(table(dataReg116$ech_id))[1],
     type="all",add=TRUE)

# plot des taux de mortalité par dose
points((nb_mtot/(nb_mtot+nb_vi))~dose,data=dataRegDos116,
       col="red",pch=19)

#plot de la valeur de la DL50
text(x=20,y=0.3,labels=paste("ED50=",round(ED50v[1],digits=4)))






#Idem pour le clone de référence du projet : 17-0153-0020

setwd("I:/LYON/R4P/TRANS REGNE/InhibiteursACCase/Test_Myzus_Spirotétramate/Data")
dataReg1<-read.table("ExtractionR_MyzusSpiro_1701530020_20241025.txt",
                     header=TRUE,stringsAsFactors=TRUE,sep=";")

#aggregating the number of individuals per rep and dose
dataRegRep1<-as.data.frame(aggregate(cbind(nb_vi,nb_mtot)~dose + dat_test,
                                     data=dataReg1,"sum"))

#aggregating the number of individuals per dose
dataRegDos1<-as.data.frame(aggregate(cbind(nb_vi,nb_mtot)~dose,
                                     data=dataReg1,"sum"))

# Création de la courbe dose réponse par dates
temp.m3<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataReg1,
             curveid=dat_test,
             fct=LN.3u(),type="binomial")
#plot de la courbe dose-réponse
plot(temp.m3,ylim=c(0,1.1),xlim=c(0,100),
     main="Courbe dose réponse du clone 17-0153-0020 \nau spriotétramate par date \n(72h après retrait des femelles)",
     ylab="mortality rate",legendPos=c(100,0.35),col = T )
#comparaison des courbes entre dates
compParm(temp.m3,"e")
#comparaison DL50 entre dates
EDcomp(temp.m3,c(50,50))
#DL 50 pour chaque date
ED50m1<-ED(temp.m3,0.50,type="absolute")


# Création de la courbe dose réponse avec
# le data où les vivants et les morts sont agrégés par date et dose
temp.m4<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataRegRep1,
             fct=LN.3u(),type="binomial")
#DL50 de la courbe dose réponse avec cumul de dates
ED50v<-ED(temp.m4,0.50,type="absolute")

#plot de la courbe "moyenne" avec intervalles de confiance
plot(temp.m4,ylim=c(0,1.1),xlim=c(0,100),
     main="Courbe moyenne dose-réponse du clone 17-0153-0020 \nau spirotétramate",
     type="confidence",ylab="mortality rate",
     las=1)

# plot des taux de mortalité par dose et par date
plot(temp.m4,ylim=c(0,1.1),xlim=c(0,100),
     main=names(table(dataReg1$ech_id))[1],
     type="all",add=TRUE)

# plot des taux de mortalité par dose
points((nb_mtot/(nb_mtot+nb_vi))~dose,data=dataRegDos1,
       col="red",pch=19)

#plot de la valeur de la DL50
text(x=20,y=0.3,labels=paste("ED50=",round(ED50v[1],digits=4)))



# Comparer sur un même graphique les données du clone de référence et du clone porteur de la mutation A2226V

#Création d'un jeu de données avec l'ensemble des données des deux clones  
data_tot <- rbind(dataReg,dataReg1)

dataRegRepTot<-as.data.frame(aggregate(cbind(nb_vi,nb_mtot)~dose + dat_test + ech_id,
                                       data=data_tot,"sum"))
dataRegRepTot$date_puceron <- paste0(dataRegRepTot$ech_id,"_",dataRegRepTot$dat_test)
temp.x<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
            weights=(nb_mtot+nb_vi),
            data=dataRegRepTot,
            curveid=date_puceron,
            fct=LN.3u(),type="binomial")

#plot de la courbe dose-réponse
plot(temp.m1, col = "red",
     main="Courbe dose réponse des clones\nau spirotétramate pour chaque date",
     ylab="mortality rate",legendPos=c(10,0.45))
plot(temp.m3, add = TRUE, col = "darkolivegreen2")


##############################################################################/
#END
##############################################################################/