##############################################################################/
##############################################################################/
#Loading the data sets and packages needed for CrossKingdom analyses
##############################################################################/
##############################################################################/

##############################################################################/
#Loading the libraries####
##############################################################################/

library(drc)
library(plotrix)
library(gdata)
library(ggplot2)
library(dplyr)


##############################################################################/
#Loading and preparing the main data set####
##############################################################################/

#load the dose response dataset
dataReg<-read.table("data/ExtractionR_MyzusSpiro_2401160001_20241025.txt",
                    header=TRUE,stringsAsFactors=TRUE,sep=";")

#aggregating the number of individuals per rep and dose
dataRegRep<-as.data.frame(aggregate(cbind(nb_vi,nb_mtot)~dose + dat_test,
                                    data=dataReg,"sum"))

#aggregating the number of individuals per dose
dataRegDos<-as.data.frame(aggregate(cbind(nb_vi,nb_mtot)~dose,
                                    data=dataReg,"sum"))


##############################################################################/
#Writing session information for reproducibility####
##############################################################################/

sink("session_info.txt")
print(sessioninfo::session_info())
sink()
#inspired by an R gist of François Briatte: 
#https://gist.github.com/briatte/14e47fb0cfb8801f25c889edea3fcd9b


##############################################################################/
#END
##############################################################################/