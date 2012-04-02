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

source("http://bioconductor.org/biocLite.R")
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

#helper function will be nested inside another mapply
#make sequence from two lists: exonStart, exonEnd
startStopToSequence <- function(x,y){
 as.vector(unlist(mapply(x,y,FUN=seq)))

}

#calculate cDNA position of genomic position.
#requires lists of genomic positions of exonStarts and exonEnds
getMrnaPos <- function(p,exS,exE){
  #calculate array of exon lengths
  exonls <- exE - exS
  #find out which exon the mutation occurs in.
  #highest start positon that is less than mutation positon
  exN <- which(exS == max(exS[p >= exS]) )
  
  #Sum the lengths of the prior exons along with length to mutation from exN.
  #This will be the length into to cDNA where the mutation position is.
  cdnal <- sum(exonls[0:(exN-1)]) + (p - exS[exN])

  return(cdnal) 
} 

#calculate sum of exon lengths
#requires lists of genomic positions of exonStarts and exonEnds
calcExonLength <- function(exS,exE){
  #calculate array of exon lengths
  exonls <- exE - exS
  #sum exon lengths
  return(sum(exonls))
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
#change name column of df.in3 (knownGene) to "transcript"
df.in3 <- plyr::rename(df.in3, c("name" = "transcript") )
#Convert sequence data to DNAString
df.in3$seq <- sapply(df.in3$seq, FUN=DNAString)
#Parse strings of digits separated by commas to array of integers
df.in3$exonstarts<-unname(sapply(df.in3$exonstarts, FUN=digitStringToArray))
df.in3$exonends<-unname(sapply(df.in3$exonends, FUN=digitStringToArray))


################################################
###  Brief Summary of Data              #######
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
sum(df.in4$transcript %in% df.in3$transcript)

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
df1 <-merge(cbind(df.in2, df2.parsed), df.in1)

# TODO: Finish work on this data set later.

#######################################################
###  Parse Mutation Strings. Data Set Two      #######
#####################################################

### Merge knownGene data with mutation data
df2 <- merge(df.in4, df.in3, by.x='transcript', by.y='transcript') 

#This will construct a sequence as long as the mrna sequence for each transcript
#However, the values in this will be the genomic position corresponding to the
#cDNA of the transcript sequence.  Verbose and wasteful.
#d<-mapply(FUN=startStopToSequence,df.in3$exonstarts,df.in3$exonends)

#get cDNA position of start of mutation.
#cdp <- getCDnaPos(df2$leftflank[[5]] +1, df2$exonstarts[[5]], df2$exonends[[5]])

#mutation genome position
mutgbpos <- df2$leftflank[5] +1
#exonstarts and ends
b <- df2$exonstarts[[5]]
c <- df2$exonends[[5]]
#calc mutation position in cDNA
mutcpos <- getMrnaPos(mutgbpos, b, c)
#calc coding sequence start and end (lengths into seq)
cdss <- getMrnaPos(df2$cdsstart[[5]],b,c)
cdse <- getMrnaPos(df2$cdsend[[5]],b,c)
#calc transcription start and end (lengths into seq)
txs <- getMrnaPos(df2$txstart[[5]],b,c)
txe <- getMrnaPos(df2$txen[[5]],b,c)

#make a copy of the genomic sequence
mutseq <- reverseComplement(df2$seq[[5]])
#check that ref_allele matches
mutseq[mutcpos]  == DNAString(df2$ref_allele[5])
#Make mutation
mutseq[mutcpos] <- DNAString(df2$var_allele[5])

#Exon lengths
c-b

#translate coding sequence of seq
mrna<-df2$seq[[5]]
cdmrna <- subseq(mrna,length(mrna)-cdse+1,length(mrna)-cdss)
aas <- (translate(dna2rna(sub1)))



ts <- extractTranscripts(Hsapiens[[df2$chrom[[5]]]], 
                   IntegerList(b+1), IntegerList(c), 
                   "+"
                   )[[1]]
df2$seq[[5]]
reverseComplement(df2$seq[[5]])

calcExonLength(b,c)
length(df2$seq[[5]])


d<-df2$seq[[5]]
e<-dna2rna((d))
f<-translate(reverse(e))



################################################################################
#        Clean Up                                                       #######
##############################################################################

### REMOVE ALL OBJECTS ### 
rm(list=ls(all=TRUE))

### Detach Packages (reverse order from load) ###
detach("package:stringr")
detach("package:plyr")
detach("package:GenomicFeatures")
detach("package:BSgenome.Hsapiens.UCSC.hg18")
detach("package:BSgenome")
detach("package:Biostrings")






