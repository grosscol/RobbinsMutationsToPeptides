################################################################################
# Title:  Examine Data Set Three: Gastric                                      #
# Date:   2012-05-08                                                           #
# Author: Colin A. Gross                                                       #
# Desc:   SnpEff output to text from VCF. Output mutated peptides              #
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
library("BSgenome.Hsapiens.UCSC.hg19")
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

#input Directory
impdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata'

#List of files to process.
# infilenames<-c(
# '2219_tumor-vs-normal_discordantSNPs_annotated_new.txt',
# '2221_tumor-vs-normal_discordantSNPs_annotated_new.txt',
# '2246_tumor-vs-normal_discordantSNPs_annotated_new.txt',
# '2359_tumor-vs-normal_discordantSNPs_annotated_new.txt',
# '2556_tumor-vs-normal_discordantSNPs_annotated_new.txt',
# '3466_tumor-vs-normal_discordantSNPs_annotated_new.txt',
# 'Gastric2_tumor-vs-normal_discordantSNPs_annotated_new.txt'
# )

infilenames <- c(
  'mutsrefs2012-04-30_141449.txt',
  '2219_mutsrefs2012-05-08_183106.txt',
  '2221_mutsrefs2012-05-08_183335.txt'
  )

#import combined data frame
infile<-paste(impdir,infilenames[1],sep='/')
#cls<- c( rep("character",3), rep("integer",5), rep("character",6), integer, rep("character",2) )
dfc <-read.table(infile, header=TRUE, sep="\t",
                   comment.char="#",encoding="UTF-8",stringsAsFactors=FALSE)

################################################
###  Format Data                        #######
##############################################

#Convert column names to lower case
colnames(dfc) <- tolower(colnames(dfc))

#change name column from knownGene table to "transcript"
dfc <- plyr::rename(dfc, c("name" = "transcript") )
#change var and ref column names to match expected column names with _allele
dfc <- plyr::rename(dfc, c("var" = "var_allele", "ref" = "ref_allele") )

#Parse strings of digits separated by commas to array of integers
dfc$exonstarts<-unname(sapply(dfc$exonstarts, FUN=digitStringToArray))
dfc$exonends<-unname(sapply(dfc$exonends, FUN=digitStringToArray))

### CALCULATE LEFTFLANK AND RIGHTFLANK ###

#get typs of mutations as lists of logical values
# var_alleles with '-' are the way deletions are marked in this data
# insertions marked with '+' in var_allele
dels <- grepl(pattern='\\-', x=dfc$var_allele)
ins <- grepl(pattern='\\+',x=dfc$var_allele)
pnts <- dfc$ref_allele != '*'

#For point mutations (given ref allele is not *)
dfc$leftflank[pnts]  <- dfc$pos[pnts] - 1
dfc$rightflank[pnts] <- dfc$pos[pnts] + 1

#For deletions
#copy var allele data to ref allele (where it should have been)
dfc$ref_allele[dels] <- substr(dfc$var_allele[dels],2,nchar(dfc$var_allele[dels]))
dfc$leftflank[dels] <- dfc$pos[dels] - 1
dfc$rightflank[dels] <- dfc$pos[dels] + nchar(dfc$ref_allele[dels])
#Set var_alleles to blank... since they are supposed to be deletions
dfc$var_allele[dels] <- ''


#For insertions (position information contextual with gene strand.)
minus <- dfc$strand=='-'
#Flank calculations for '+' strand genes.
dfc$leftflank[ins] =  dfc$pos[ins] -1
dfc$rightflank[ins] =  dfc$pos[ins]
# #Flanks for '-' Strand Insertions (seems odd)
dfc$leftflank[ins & minus] = dfc$pos[ins & minus] 
dfc$rightflank[ins & minus] = dfc$pos[ins & minus] +1
#remove '+' from the var allele
dfc$var_allele[ins] <- substr(dfc$var_allele[ins],2,nchar(dfc$var_allele[ins]))
#Set ref_allele to blank... since it is supposed to be blank
dfc$ref_allele[ins] <- ''


################################################################################
###  CALCULATIONS                                                       #######
##############################################################################

#dfc is input data frame. 
# Expected initial columns required:
# transcript chrom strand txstart txend cdsstart cdsend exoncount
# exonstarts exonends proteinid alignid seq ref_allele var_allele
# leftflank rightflank

#rename input data frame to "d" for the sake of brevity
d <- dfc
d$seq <- NULL #seq column is superfluous

#InO(Intermediat Output)
print(paste(nrow(d), "Rows total."))

##########
### 1  ### Check if mutation is outside of coding region
##########
d$noncoding <- d$leftflank < d$cdsstart | d$leftflank >= d$cdsend
# InO (Intermediate Output)
print(paste(sum(d$noncoding),"Rows w/ leftflank outside cdsstart and cdsend:"))

##########
### 2  ### Check if mutation is inside an exon segment
##########
#Pass in list of starts, ends, and a position
inExons <- function(exonstarts,exonends,leftflank){
  max(which(leftflank >= unlist(exonstarts))) == min(which(leftflank < unlist(exonends)))
}
#Apply the inExons() function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$inexon<-apply(d[,c('exonstarts','exonends','leftflank')], MARGIN=1, FUN=splat(inExons))
#InO
print(paste(sum(d$inexon),"Rows w/ leftflank in exons."))

##########
### 3  ### Create criteria for further processing
##########
#mutation is in an exon and mutation within coding region
L1crite <- d$inexon & !d$noncoding
#InO
print(paste(sum(L1crite),"Rows w/ leftflank in exons AND leftflank in coding region."))

##########
### 4  ### Calculate sum of exon lengths
##########
#requires lists of genomic positions of exonStarts and exonEnds
sumExonLengths <- function(exonstarts,exonends){
  #calculate array of exon lengths
  exonlens <- exonends - exonstarts
  #sum exon lengths
  return(sum(exonlens))
}
#Apply the sumExonLengths() function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$exonslen <- apply(d[,c('exonstarts','exonends')],MARGIN=1,FUN=splat(sumExonLengths))
#InO
print(paste(sum(is.na(d$exonslen)),"Rows w/ NA result of sum exon lengths."))

##########
### 5  ### Calculate exon in which coding begins
##########
getCodingStartExon <- function(exonstarts,cdsstart){
  #highest start positon that is less than cdsstart positon
  which(exonstarts == max(exonstarts[cdsstart >= exonstarts]) )
}
#Apply the getCodingExonStart() function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$exCDS <- apply(d[,c('exonstarts','cdsstart')],MARGIN=1,FUN=splat(getCodingStartExon))


##########
### 5  ### Calculate exon in which coding begins
##########
getCodingEndExon <- function(exonends,cdsend){
  #lowest end positon that is greater than cdsend positon
  which(exonends == min(exonends[cdsend <= exonends]) )
}
#Apply the getCodingExonStart() function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$exCDE <- apply(d[,c('exonends','cdsend')],MARGIN=1,FUN=splat(getCodingEndExon))

##########
### 6  ### Get starts for coding exons (as 1-based index)
##########
getCodingExonStarts <- function(exonstarts,exCDS,exCDE,cdsstart){
  #exonstarts are 0 based. Add 1 to convert to 1-based position
  ret<-exonstarts[exCDS:exCDE] + 1 
  #replace first coding exon start with coding start (and + 1 for conversion)
  ret[1] <- cdsstart + 1
  ret
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$codingStarts <- apply(d[,c('exonstarts','exCDS','exCDE','cdsstart')],
                        MARGIN=1,FUN=splat(getCodingExonStarts))

##########
### 7  ### Get ends for coding exons.
##########
getCodingExonEnds <- function(exonends,exCDS,exCDE,cdsend){
  #exonends are already 1-based position index.
  ret<-exonends[exCDS:exCDE]
  #replace last coding exon end with the coding end position (alredy 1-based)
  ret[length(ret)]<-cdsend
  ret
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$codingEnds <- apply(d[,c('exonends','exCDS','exCDE','cdsend')],
                      MARGIN=1,FUN=splat(getCodingExonEnds))

###########################################################
### Begin Using Criteria to Only Analyze Coding Mutants ###
###########################################################
#Use L1 criteria

                  MARGIN=1,FUN=splat(modifyCodingEnds))



##########
### 8  ### Get all the coding starts prior to mutation
##########
getPreMutCodingStarts <- function(codingStarts,leftflank){
  codingStarts[codingStarts <= (leftflank+1)]
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$premutstarts[L1crite] <- apply(d[L1crite, c('codingStarts','leftflank')],
                        MARGIN=1, FUN=splat(getPreMutCodingStarts))

##########
### 9  ### Get all the coding ends prior to mutation (append mut here)
##########
getPreMutCodingEnds <- function(codingEnds,leftflank){
  #Do not include the coding end if it is the mutation position.
  ret <- codingEnds[codingEnds < (leftflank+1)]
  #append mutation start position to this list
  ret[length(ret)+1] <- leftflank+1
  ret
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$premutends[L1crite] <- apply(d[L1crite, c('codingEnds','leftflank')],
                                MARGIN=1, FUN=splat(getPreMutCodingEnds))

##########
### 10 ### Check for mismatches between number of pre mutation starts and ends
##########
flagNumStartStops <- function(premutstarts,premutends){
  length(premutstarts) != length(premutends)
}
d$flagMutSE[L1crite] <- apply(d[L1crite, c('premutstarts','premutends')],
                             MARGIN=1, FUN=splat(flagNumStartStops))
#InO
print(paste(sum(d[L1crite,'flagMutSE']),
            "ERROR flag(s) thrown from num of starts & stops mismatch."))

##########
### 11 ### Calculate the mutation position relative to the reference transcript
##########
getMutTrnscrtDNAPosition <- function(premutstarts,premutends){
  #length calculations have to accomodated 1 based start (subtract 1 from each)
  sum(premutends - (premutstarts - 1 ))
}
d$mutTrnscrtDNAPos[L1crite] <- apply(d[L1crite, c('premutstarts','premutends')],
                            MARGIN=1, FUN=splat(getMutTrnscrtDNAPosition))


##########
### 12 ### Get transcript reference DNA from UCSC genome
##########
getTranscriptDNA <- function(codingStarts,codingEnds,chrom){
  dnas <- DNAStringSet(unmasked(Hsapiens[[chrom]]), 
                       start=as.vector(codingStarts), end=as.vector(codingEnds) )
  #unlist the segements to create a single transcript
  dnas <- unlist(dnas)
}
#get start time
start.time <- Sys.time()
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$trnscrtDNA[L1crite] <- apply(d[L1crite,c('codingStarts','codingEnds','chrom')],
                               MARGIN=1,FUN=splat(getTranscriptDNA))
#compute run duration
dur <- Sys.time() - start.time
#InO
print(paste(dur,attr(dur,'units'),"required to do reference DNA lookup."))

##########
### 13 ### Calculate mutatant amino acid position
##########
getMutAAPosition <- function(mutTrnscrtDNAPos,strand,trnscrtDNA){
  if(strand=='+'){
    ret <- ceiling(mutTrnscrtDNAPos / 3)
  }else{
    ret <- (length(trnscrtDNA)/3) - (ceiling(mutTrnscrtDNAPos/3)) + 1
  }
  ret
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$mutAAPos[L1crite] <- apply(d[L1crite,c('mutTrnscrtDNAPos','strand','trnscrtDNA')],
                               MARGIN=1,FUN=splat(getMutAAPosition))

#do adjustment (+1) for insertions at start of codon '+' strand genes
crite <- L1crite & d$strand=='+' & d$ref_allele == '' & (d$mutTrnscrtDNAPos %% 3 == 0)
d$mutAAPos[crite] <- d$mutAAPos[crite] + 1
#InO
print(paste(sum(crite),"insertions at start of codon on '+' strand genes"))

#do adjustment (+1) for insertions at start of codon '-' strand genes
crite <- L1crite & d$strand=='-' & d$ref_allele == '' & 
  ((length(d$trnscrtDNA[crite]) - d$mutTrnscrtDNAPos[crite]+1) %% 3 == 0)
d$mutAAPos[crite] <- d$mutAAPos[crite] + 1
#InO
print(paste(sum(crite),"insertions at start of codon on minus strand genes"))


##########
### 14 ### Make mutation and store mutant sequence
##########
getMutantRefDNA <- function(trnscrtDNA,lpmut,rpmut,var_allele){
  #Concatenation of [left side],[variant allele],[right side]
    xscat(
    substr(trnscrtDNA, 1, lpmut), 
    DNAString(var_allele), 
    substr(trnscrtDNA, rpmut, length(trnscrtDNA))
    )
}
#calculate left and right sides of variant sequence
d$lpmut[L1crite] <- d$mutTrnscrtDNAPos[L1crite] - 1
d$rpmut[L1crite] <- d$mutTrnscrtDNAPos[L1crite] + nchar(d$ref_allele[L1crite])
#get start time
start.time <- Sys.time()
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$mutTrnscrtDNA[L1crite] <- apply(d[L1crite,c('trnscrtDNA','lpmut','rpmut','var_allele')],
                               MARGIN=1,FUN=splat(getMutantRefDNA))
#compute run duration
dur <- Sys.time() - start.time
#InO
print(paste(dur,attr(dur,'units'),"required to do concatenate mutant ref DNA."))

##########
### 15 ### Do translations. Not elegant, but ham handed works.
##########
#get start time
start.time <- Sys.time()
#translate normal '+' strand genes
crite <- L1crite & d$strand=='+'
d$aanorm[crite] <- sapply(d$trnscrtDNA[crite],
                          FUN=function(x){translate(dna2rna(x))})
#translate normal '-' strand genes
crite <- L1crite & d$strand=='-'
d$aanorm[crite] <- sapply(d$trnscrtDNA[crite],
                          FUN=function(x){translate(dna2rna(reverseComplement(x)))})
#translate mutant '+' strand genes
crite <- L1crite & d$strand=='+'
d$aamut[crite] <- sapply(d$mutTrnscrtDNA[crite],
                          FUN=function(x){translate(dna2rna(x))})
#translate mutant '-' strand genes
crite <- L1crite & d$strand=='-'
d$aamut[crite] <- sapply(d$mutTrnscrtDNA[crite],
                          FUN=function(x){translate(dna2rna(reverseComplement(x)))})
#compute run duration
dur <- Sys.time() - start.time
#InO
print(paste(dur,attr(dur,'units'),"required to do transcription & translation of all ref DNA."))


##########
### 16 ### Check if variant amino acid is a stop.
##########
isMutAAstop <- function(aamut, mutAAPos){
  #Check for edge/error cases
  if(mutAAPos > length(aamut)){
    cat(length(aamut),mutAAPos,'; ')
    return(NA)
  }
  aamut[mutAAPos] == AAString("*")
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$isTrunc[L1crite] <- apply(d[L1crite,c('aamut','mutAAPos')],
                                  MARGIN=1,FUN=splat(isMutAAstop))
#InO
print(paste(sum(d$isTrunc[L1crite & !is.na(d$isTrunc)]),
            "cases with mutant AA position is stop. Truncatation"))
print(paste(sum(is.na(d$isTrunc[L1crite])),
            "cases with mutant AA position > length mutant peptide."))
#For the time being, flag these as truncated. (TODO: Handle these cases)
crite <- is.na(d$isTrunc) & L1crite
d$isTrunc[crite] <- TRUE

###########################################################
### Add Truncation to Coding Criteria and Apply to      ###
###  Subsequent Analysis                                ###
###########################################################
L2crite <- L1crite & !d$isTrunc
#Convert NA to false
L2crite[is.na(L2crite)] <- FALSE
#Use L2 criteria


##########
### 17 ### Check variant peptide against normal peptide
##########
isSynonymousMutation <- function(mutAAPos,aamut,aanorm){
  aamut[mutAAPos:length(aamut)] == aanorm[mutAAPos:length(aanorm)]
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$isSynon[L2crite] <- apply(d[L2crite,c('aanorm','aamut','mutAAPos')],
                            MARGIN=1,FUN=splat(isSynonymousMutation))
#InO
print(paste(sum(d$isSynon[L2crite]),"Synonymous mutations flagged."))
print(paste(sum(!d$isSynon[L2crite]),"Non-Synonymous mutations."))

###########################################################
### Switch to Non-Synonymous Mutation Criteria          ###
###  for Subsequent Analysis                            ###
###########################################################
L3crite <- !d$isSynon
#Convert NA to false
L3crite[is.na(L3crite)] <- FALSE
#Use L2 criteria

##########
### 18 ### Calculate report cut length to left of variant AA
##########
#default to position 1
d$lareport[L3crite] <- 1
#for those with a mutation position > 10, report 10 AA to left of variant
crite <- d$mutAAPos > 10 & L3crite
d$lareport[crite] <- d$mutAAPos[crite] - 10

##########
### 19 ### Calculate report cut length to left of variant AA
##########
#Calc length of mutant amino acid sequence (aamut)
d$lenAAMut[L3crite] <- sapply(d$aamut[L3crite],length)
#default to end position
d$rareport[L3crite] <- d$lenAAMut[L3crite]

crite <- (d$lenAAMut - d$mutAAPos) > 10 & L3crite
d$rareport[crite] <- d$mutAAPos[crite] + 10

##########
### 20 ### Check for frame shift
##########
#Create criteria to identify relevant frame shift mutations
FScrite <- (nchar(d$var_allele) != nchar(d$ref_allele)) & 
  ( (nchar(d$var_allele)-nchar(d$ref_allele))%%3 !=0 ) & L3crite
#InO
print(paste(sum(crite),"Non-synonymous frame shift mutations flagged."))

##########
### 21 ### Recalculate right amino acid report cut length
##########
getNewRightAAReport <- function(aamut, rareport){
  #dangerously assuming first stop codon found will be beyond mut position
  stoppos <- regexpr('\\*',aamut)
  if(stoppos == -1){
    #no stop codon found.  TODO:translate into 3'UTR ?
    #FLAG using -1
    return(stoppos)
  }else{
    #Cut at the star (will include the '*')
    return(stoppos[1])
    #set cut position to 1 before first stop codon
    #return(stoppos[1]-1)
  }
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$rareport[FScrite] <- apply(d[FScrite,c('aamut','rareport')],
                            MARGIN=1,FUN=splat(getNewRightAAReport))
#Flag stop loss and set rareport to end of transcript
d$isStopLoss <- FALSE
d$isStopLoss[which(d$rareport==-1)] <- TRUE


##########
### 22 ### Cut the mutant peptide and store result for reporting
##########
getMutantAAReportSequence <- function(aamut,lareport,rareport){
  subseq(aamut,start=lareport,end=rareport)
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$mutaaReport[L3crite] <- apply(d[L3crite,c('aamut','lareport','rareport')],
                            MARGIN=1,FUN=splat(getMutantAAReportSequence))

#InO
print(paste(sum(L3crite),"reported mutant amino acid sequences"))



##########
### 23 ### Make small subset of unique mutant peptides
##########
#dupes <- duplicated(sapply(d$mutaaReport[L3crite],as.character))



################################################
###  Output Result of    Calculations   #######
##############################################
#Copy data frame for output
d.o <- d

#debug compare old and new
dnew <-  sapply(d.o$aanorm[L3crite],as.character)
dold <-  sapply(d$aanorm[L3crite],as.character)


dold <-  sapply(dold, function(x){substr(x, nchar(x)-10,nchar(x)) })
dnew <-  sapply(dnew, function(x){substr(x, nchar(x)-11,nchar(x)) })

comp<-data.frame(oldAAPos=d$mutAAPos[L3crite], newAAPoss=d.o$mutAAPos[L3crite],
                 strand=d$strand[L3crite],aaOld=dold,aaNew=dnew,row.names=NULL)
View(comp)
View(unname(cbind(dold,dnew)))

#Convert Biostrings to regular character strings
d.o$trnscrtDNA <- sapply(d.o$trnscrtDNA,as.character)
d.o$mutTrnscrtDNA <- sapply(d.o$mutTrnscrtDNA,as.character)
d.o$aanorm <- sapply(d.o$aanorm,as.character)
d.o$aamut <- sapply(d.o$aamut,as.character)
d.o$mutaaReport <- sapply(d.o$mutaaReport,as.character)
View(d.o[L3crite & d.o$transcript=='uc003lli.3',c(1,3,15,16,42,35,36)])


################################################################################
#        Clean Up                                                       #######
##############################################################################

### REMOVE ALL OBJECTS ### 
rm(list=ls(all=TRUE))

### Detach Packages (reverse order from load) ###
detach("package:reshape")
detach("package:stringr")
detach("package:plyr")
detach("package:GenomicFeatures")
detach("package:BSgenome.Hsapiens.UCSC.hg19")
detach("package:BSgenome")
detach("package:Biostrings")



################################################################################
#        Scrap Code                                                     #######
##############################################################################

# #Include these just before "Get all the coding starts prior to mutation"  
# #to add 3 nucleotides to the 3' end of the resulting mrna transcript
# ##########
# ### NEW ## Modify coding starts of minus strand genes to include stop codon 
# ##########
# modifyCodingStarts <- function(codingStarts,strand){
#   #replace first coding start with location 3 towards the 5' end.
#   #will force inclusion of stop codon for minus strand genes.
#   codingStarts[1] <- codingStarts[1] - 3
#   codingStarts
# }
# #Need to ensure that normal stop codon is included in sequence for (-) genes
# crite <- d$strand=='-' & L1crite
# #Apply the function to each row of the data frame
# #Use splat() instead of spelling out function arguments
# d$codingStarts[crite] <- apply(d[crite,c('codingStarts','strand')],
#                         MARGIN=1,FUN=splat(modifyCodingStarts))
# 
# ##########
# ### NEW ## Modify coding ends of plus strand genes to include stop codon 
# ##########
# modifyCodingEnds <- function(codingEnds,strand){
#   #replace first coding start with location 3 towards the 5' end.
#   #will force inclusion of stop codon for minus strand genes.
#   codingEnds[length(codingEnds)] <- codingEnds[length(codingEnds)] + 3
#   codingEnds
# }
# #Need to ensure that normal stop codon is included in sequence for (+) genes
# crite <- d$strand=='+' & L1crite
# #Apply the function to each row of the data frame
# #Use splat() instead of spelling out function arguments
# d$codingEnds[crite] <- apply(d[crite,c('codingEnds','strand')],








