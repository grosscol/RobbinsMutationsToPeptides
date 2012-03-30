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
library("stringr")

################################################################################
#        MAIN                                                                  #
################################################################################

#############
# Function  #
#   Defs    #
#############

#helper function will be nested inside of an lapply or sapply call
#take string of digitis separated by commas, return list of integers.
digitStringToArray <- function(x){
  ret<-lapply(str_extract_all(x,'\\d+'),FUN=as.integer, names=NULL)
  unname(ret)
}

################################################
###  Import Data                        #######
##############################################

###########################
# Transcript + Mutations  #
###########################

#Import Directory
impdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata/'

### Data Set One from Dr. Robbins ###
#Transcript + mRNA 
infile<-paste(impdir,'mrna2012-03-26_172605.txt',sep='')
cls<- c( "character","character" )
df.in1<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                     comment.char="#",encoding="UTF-8",stringsAsFactors=FALSE)

#Transcript + mutations 
infile<-paste(impdir,'muts2012-03-26_172605.txt',sep='')
cls<- c( "character","character" )
df.in2<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                     comment.char="#",encoding="UTF-8",stringsAsFactors=FALSE)

### Data Set Two from Dr. Robbins ###
#Transcript + annotation + mRNA
infile<-paste(impdir,'mrna2012-03-30_140733.txt',sep='')
cls<- c( rep("character",3), rep("integer",5), rep("character",5) )
df.in3<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                     comment.char="#",encoding="UTF-8",stringsAsFactors=FALSE)

#Transcript + mutations
infile<-paste(impdir,'muts2012-03-30_140733.txt',sep='')
cls<- c( "character","character","integer","integer","character","character" )
df.in4<-read.table(infile, header=TRUE, sep="\t", colClasses=cls, 
                     comment.char="#",encoding="UTF-8",stringsAsFactors=FALSE)

################################################
###  Data Frame Formatting              #######
##############################################
### Using data set one ###
#make column names all lower case
colnames(df.in1) <- tolower(colnames(df.in1))
colnames(df.in2) <- tolower(colnames(df.in2))
#Convert sequence data to DNAString
df.in1$seq <- sapply(df.in1$seq, FUN=DNAString)

### Using data set one ###
#make column names all lower case
colnames(df.in3) <- tolower(colnames(df.in3))
colnames(df.in4) <- tolower(colnames(df.in4))
#Convert sequence data to DNAString
df.in3$seq <- sapply(df.in3$seq, FUN=DNAString)
#Parse strings of digits separated by commas to array of integers
df.in3$exonstarts<-unname(sapply(df.in3$exonstarts, FUN=digitStringToArray))
df.in3$exonends<-unname(sapply(df.in3$exonends, FUN=digitStringToArray))


################################################
###  Summary of Strings                 #######
##############################################

### Using data set one ###
#Names of transcripts given vs names of transcripts given also had cDNA seqs
length(df.in2$transcript)
sum(df.in2$transcript %in% df.in1$name)

#Number of mutations that occur in the same transcript name
sum(duplicated(df.in2$transcript))

#get base transcript names (first 8 characters)
btn1<-sapply(df.in2$transcript, FUN=substr, start=1, stop=8)
#number of uniqe base names, number of duplicated, and which were dupes.
length(unique(btn1))
sum(duplicated(btn1))
btn1[duplicated(btn1)]

### Using data set two ###
length(df.in4$transcript)
sum(df.in4$transcript %in% df.in3$name)

#Number of mutations that occur in the same transcript name
sum(duplicated(df.in4$transcript))

#get base transcript names (first 8 characters)
btn2<-sapply(df.in4$transcript, FUN=substr, start=1, stop=8)
#number of uniqe base names, number of duplicated, and which were dupes.
length(unique(btn2))
sum(duplicated(btn2))
btn2[duplicated(btn2)]

#######################################################
###  Parse Mutation Strings. Data Set One      #######
#####################################################

#Regular expressions for getting specific parts of cdna.change column
regoldbase <- '^c\\.(\\w)\\d+\\w$' #will always be the 3 character
regnewbase <- '^c\\.\\w\\d+(\\w)$' #will always be the 4 to len - 1 chars
regbaseloc <- '^c\\.\\w(\\d+)\\w$' #will always be the last char

#Instead, just use known positional format to get values.
df2.parsed<-data.frame(
    name       =df.in2$Transcript,
    basescript =substr(df.in2$transcript, start=1,stop=8),
    cdnachange =df.in2$cdna.change,    
    oldbase    =substr(df.in2$cdna.change, start=3,stop=3),
    loc        =substr(df.in2$cdna.change,start=4,
                       stop=nchar(df.in2$cdna.change)-1),
    newbase    =substr(df.in2$cdna.change,start=nchar(df.in2$cdna.change),
                       stop=nchar(df.in2$cdna.change))
    )

#Column bind parsed mutation info to orginal input and merge with mrna data
df1.muts<-merge(cbind(df.in2, df2.parsed), df.in1)

# TODO: Finish work on this data set later.

#######################################################
###  Parse Mutation Strings. Data Set Two      #######
#####################################################



################################################################################
#        Clean Up                                                       #######
##############################################################################

### REMOVE ALL OBJECTS ### 
rm(list=ls(all=TRUE))

### Detach Packages (reverse order from load) ###
detach("package:stringr")
detach("package:plyr")
detach("package:BSgenome.Hsapiens.UCSC.hg18")
detach("package:Biostrings")






