##############################################################################/
##############################################################################/
#Loading the data sets and packages needed for CrossKingdom analyses
##############################################################################/
##############################################################################/

##############################################################################/
#Loading the libraries####
##############################################################################/

library(drc)
library(ape)


##############################################################################/
#Loading and preparing the main data set####
##############################################################################/

#data by individuals, including all individuals (hd, families, parents and CC)
RTdata<-read.table("data/datatot.txt",sep="\t",stringsAsFactors=FALSE,
                   header=TRUE)
#because some functions do not like "." within colnames, we replace them
colnames(RTdata)<-str_replace_all(colnames(RTdata),"[.]","_")


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