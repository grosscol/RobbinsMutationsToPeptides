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

#Transcript + mRNA 
infile<-paste(impdir,'mrna2012-03-26_172605.txt',sep='')
cls<- c( "character","character" )
df.mrnas<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                     comment.char = "#", encoding="UTF-8", stringsAsFactors=FALSE)

#Transcript + mutations 
infile<-paste(impdir,'muts2012-03-26_172605.txt',sep='')
cls<- c( "character","character" )
pgdxmuts<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                    comment.char = "#", encoding="UTF-8", stringsAsFactors=FALSE)


################################################
###  Summary of Strings                 #######
##############################################

#Names of transcripts given vs names of transcripts given were also in UCSC DB
length(pgdxmuts$Transcript)
sum(pgdxmuts$Transcript %in% df.mrnas$name)

#Number of mutations that occur in the same transcript name
sum(duplicated(pgdxmuts$Transcript))

sapply(pgdxmuts$Transcript, FUN=nchar)

#get bast transcript names (first 8 characters)
btn<-sapply(pgdxmuts$Transcript, FUN=substr, start=1, stop=8 )
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
    name=pgdxmuts$Transcript,
    basescript=substr(pgdxmuts$Transcript, start=1,stop=8),
    cdnachange=pgdxmuts$cDNA.change,    
    oldbase=substr(pgdxmuts$cDNA.change, start=3,stop=3),
    loc=substr(pgdxmuts$cDNA.change,start=4,stop=nchar(pgdxmuts$cDNA.change)-1),
    newbase=substr(pgdxmuts$cDNA.change,start=nchar(pgdxmuts$cDNA.change),stop=nchar(pgdxmuts$cDNA.change))
    )

#Column bind parsed mutation info to orginal input and merge with mrna data
df.muts<-merge(cbind(pgdxmuts, pgdxparsed), df.mrnas)

#Get genomic DNA from UCSC hg18 database


seqnames(Hsapiens)

################################################################################
#        Clean Up                                                       #######
##############################################################################
detach("package:plyr")
detach("package:Biostrings")
### REMOVE ALL OBJECTS ### 
rm(list=ls(all=TRUE)) #uncomment to enable






