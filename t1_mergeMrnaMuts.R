################################################################################
# Title:  Merge mRNA and Mutations Data                                        #
# Date:   2012-03-26                                                           #
# Author: Colin A. Gross                                                       #
# Desc:   Read in initial transcript + mutation data and combine with the      #
#           querried mRNA data. Write combined data                            #
################################################################################



################################################################################
#        SETUP                                                                 #
################################################################################

#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

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
library("BSgenome.Hsapiens.UCSC.hg18")
library("plyr")

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

### Data Set One from Dr. Robbins
#Transcript + mRNA 
infile<-paste(impdir,'mrna2012-03-26_172605.txt',sep='')
cls<- c( "character","character" )
df.in1<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                     comment.char = "#", encoding="UTF-8", stringsAsFactors=FALSE)

#Transcript + mutations 
infile<-paste(impdir,'muts2012-03-26_172605.txt',sep='')
cls<- c( "character","character" )
df.in2<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                    comment.char = "#", encoding="UTF-8", stringsAsFactors=FALSE)

### Data Set Two from Dr. Robbins
#Transcript + annotation + mRNA
infile<-paste(impdir,'supplementarytable3.txt',sep='')
cls<- c( rep("character", 8) )
df.in3<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                     comment.char = "#", encoding="UTF-8", stringsAsFactors=FALSE)

#Transcript + mutations
infile<-paste(impdir,'supplementarytable3.txt',sep='')
cls<- c( rep("character", 8) )
df.in4<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                   comment.char = "#", encoding="UTF-8", stringsAsFactors=FALSE)


################################################
###  Summary of Strings                 #######
##############################################

#Names of transcripts given vs names of transcripts given were also in UCSC DB
length(df.in2$Transcript)
sum(df.in2$Transcript %in% df.in1$name)

#Number of mutations that occur in the same transcript name
sum(duplicated(df.in2$Transcript))

#get base transcript names (first 8 characters)
btn<-sapply(df.in2$Transcript, FUN=substr, start=1, stop=8 )
length(unique(btn))
sum(duplicated(btn))
btn[duplicated(btn)]

################################################
###  Parse Mutation Strings             #######
##############################################

#Regular expressions for getting specific parts
regoldbase <- '^c\\.(\\w)\\d+\\w$' #will always be the 3 character
regnewbase <- '^c\\.\\w\\d+(\\w)$' #will always be the 4 to len - 1 chars
regbaseloc <- '^c\\.\\w(\\d+)\\w$' #will always be the last char

#Instead, just use known positional format to get values.
pgdxparsed<-data.frame(
    name=df.in2$Transcript,
    basescript=substr(df.in2$Transcript, start=1,stop=8),
    cdnachange=df.in2$cDNA.change,    
    oldbase=substr(df.in2$cDNA.change, start=3,stop=3),
    loc=substr(df.in2$cDNA.change,start=4,stop=nchar(df.in2$cDNA.change)-1),
    newbase=substr(df.in2$cDNA.change,start=nchar(df.in2$cDNA.change),stop=nchar(df.in2$cDNA.change))
    )

#Column bind parsed mutation info to orginal input and merge with mrna data
df.muts<-merge(cbind(df.in2, pgdxparsed), df.in1)

#Get genomic DNA from UCSC hg18 database
seqnames(Hsapiens)

################################################################################
#        Clean Up                                                       #######
##############################################################################

### REMOVE ALL OBJECTS ### 
rm(list=ls(all=TRUE))

### Detach Packages (reverse order from load) ###
detach("package:plyr")
detach("BSgenome.Hsapiens.UCSC.hg18")
detach("package:Biostrings")






