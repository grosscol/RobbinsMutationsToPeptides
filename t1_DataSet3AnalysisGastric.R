################################################################################
# Title:  Examine Data Set Three: Gastric                                      #
# Date:   2012-04-04                                                           #
# Author: Colin A. Gross                                                       #
# Desc:   Sort out issues transcribing and translating mutation data into      #
#           21-mer peptide sequences for robbins                               #
################################################################################



################################################################################
#        SETUP                                                                 #
################################################################################

#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#biocLite("GenomicFeatures")

#source("http://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg18")
#biocLite("BSgenome")


### SETUP AN OUTPUT DIRECTORY VARIABLE ###
myOutDir <- "S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata/"

### SETUP OUTPUT SINK ### uncomment to enable
#sink(paste(myOutDir,"Rawoutput.txt", sep=""))

### SET WORKING DIRECTORY ###
setwd("C:/Users/grossco/Documents/Devel/Rwork/working")


library("Biostrings")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg18")
library("GenomicFeatures")
library("plyr")
library("stringr")
library("reshape")


################################################################################
#        MAIN                                                                  #
################################################################################

#############
# Function  #
#   Defs    #
#############

################################################
###  Import Data                        #######
##############################################

###########################
# Transcript + Mutations  #
###########################

#Import Directory
impdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata/'

















