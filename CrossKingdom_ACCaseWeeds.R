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






##############################################################################/
#END
##############################################################################/





