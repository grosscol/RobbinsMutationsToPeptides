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
  ret <- lapply(str_extract_all(x,'\\d+'),FUN=as.integer, names=NULL)
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

infilenames <- c(
  '2219_mutsrefs_2012-06-04.txt',
  '2221_mutsrefs_2012-06-04.txt',
  '2246_mutsrefs_2012-06-04.txt',
  '2359_mutsrefs_2012-06-04.txt',
  '2556_mutsrefs_2012-06-04.txt',
  '3466_mutsrefs_2012-06-04.txt',
  'gastric/gastric2_mutsrefs_2012-06-04.txt'
  )

#import combined data frame
infile<-paste(impdir,infilenames[7],sep='/')
#cls<- c( rep("character",3), rep("integer",5), rep("character",6), integer, rep("character",2) )
dfc <-read.table(infile, header=TRUE, sep="\t",
                   comment.char="#",encoding="UTF-8",stringsAsFactors=FALSE)

#output Prefix
outprfx <- ''

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
#begin output sink
sink(paste(myOutDir,outprfx,"mutsToPeps_AnalysisSink.txt", sep=""),type = c("output", "message"))
print("Beginning Mutations to Peptides Analysis")
print(infile)

##########
### 0  ### Record Start Time
##########
analysis.start <- Sys.time()
print(analysis.start)

#dfc is input data frame. 
# Expected initial columns required:
# transcript chrom strand txstart txend cdsstart cdsend exoncount
# exonstarts exonends proteinid alignid seq ref_allele var_allele
# leftflank rightflank

#rename input data frame to "d" for the sake of brevity
d <- dfc
d$seq <- NULL #if sequence is present, drop the column. Superfluous.

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
### 15 ### Store mutant & normal DNA transcript lengths
##########
d$mutTrnscrtDNAlen[L1crite] <- sapply(d$mutTrnscrtDNA[L1crite], length)
d$var_allele_len[L1crite] <- nchar(d$var_allele[L1crite])
d$ref_allele_len[L1crite] <- nchar(d$ref_allele[L1crite])

##########
### 16 ### Do translations. Not elegant, but ham handed works.
##########
#get start time
start.time <- Sys.time()
#translate normal '+' strand genes
crite <- L1crite & d$strand=='+'
d$aaseqnorm[crite] <- sapply(d$trnscrtDNA[crite],
                          FUN=function(x){translate(dna2rna(x))})
#translate normal '-' strand genes
crite <- L1crite & d$strand=='-'
d$aaseqnorm[crite] <- sapply(d$trnscrtDNA[crite],
                          FUN=function(x){translate(dna2rna(reverseComplement(x)))})
#translate mutant '+' strand genes
crite <- L1crite & d$strand=='+'
d$aaseqmut[crite] <- sapply(d$mutTrnscrtDNA[crite],
                          FUN=function(x){translate(dna2rna(x))})
#translate mutant '-' strand genes
crite <- L1crite & d$strand=='-'
d$aaseqmut[crite] <- sapply(d$mutTrnscrtDNA[crite],
                          FUN=function(x){translate(dna2rna(reverseComplement(x)))})
#compute run duration
dur <- Sys.time() - start.time
#InO
print(paste(dur,attr(dur,'units'),"required to do transcription & translation of all ref DNA."))


##########
### 17 ### Check if variant amino acid is a stop.
##########
isMutAAstop <- function(aaseqmut, mutAAPos){
  #Check for edge/error cases
  if(mutAAPos > length(aaseqmut)){
    cat(length(aaseqmut),mutAAPos,'; ')
    return(NA)
  }
  aaseqmut[mutAAPos] == AAString("*")
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$isTrunc[L1crite] <- apply(d[L1crite,c('aaseqmut','mutAAPos')],
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
### 18 ### Check variant peptide against normal peptide
##########
isSynonymousMutation <- function(mutAAPos,aaseqmut,aaseqnorm){
  aaseqmut[mutAAPos:length(aaseqmut)] == aaseqnorm[mutAAPos:length(aaseqnorm)]
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$isSynon[L2crite] <- apply(d[L2crite,c('aaseqnorm','aaseqmut','mutAAPos')],
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
### 19 ### Calculate report cut length to left of variant AA
##########
#default to position 1
d$lareport[L3crite] <- 1
#for those with a mutation position > 10, report 10 AA to left of variant
crite <- d$mutAAPos > 10 & L3crite
d$lareport[crite] <- d$mutAAPos[crite] - 10

##########
### 20 ### Calculate report cut length to left of variant AA
##########
#Calc length of mutant amino acid sequence (aaseqmut)
d$lenAAMut[L3crite] <- sapply(d$aaseqmut[L3crite],length)
#default to end position
d$rareport[L3crite] <- d$lenAAMut[L3crite]
#check of right position is more than three away from end.
crite <- (d$lenAAMut - d$mutAAPos) > 10 & L3crite
d$rareport[crite] <- d$mutAAPos[crite] + 10

##########
### 21 ### Check for frame shift
##########
#Create criteria to identify relevant frame shift mutations
FScrite <- (nchar(d$var_allele) != nchar(d$ref_allele)) & 
  ( (nchar(d$var_allele)-nchar(d$ref_allele))%%3 !=0 ) & L3crite
#InO
print(paste(sum(FScrite),"Non-synonymous frame shift mutations flagged."))

##########
### 22 ### Recalculate right amino acid report cut length for frame shifts
##########
getNewRightAAReport <- function(aaseqmut, rareport){
  #dangerously assuming first stop codon found will be beyond mut position
  stoppos <- regexpr('\\*',aaseqmut)
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
d$rareport[FScrite] <- apply(d[FScrite,c('aaseqmut','rareport')],
                            MARGIN=1,FUN=splat(getNewRightAAReport))
#Flag no stop found
d$noStop <- FALSE
d$noStop[which(d$rareport==-1)] <- TRUE
#Set right reporting end to end of peptide for those without stop codon
d$rareport[d$noStop] <- d$lenAAMut[d$noStop]

##########
### 23 ### Cut the mutant peptide and store result for reporting
##########
getMutantAAReportSequence <- function(aaseqmut,lareport,rareport){
  subseq(aaseqmut,start=lareport,end=rareport)
}
#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$mutaaReport[L3crite] <- apply(d[L3crite,c('aaseqmut','lareport','rareport')],
                            MARGIN=1,FUN=splat(getMutantAAReportSequence))

#InO
print(paste(sum(L3crite),"reported mutant amino acid sequences"))

##########
### 24 ### Get Ref and Var AA @ mutant position
##########
getAAwildtype <- function(aaseqnorm,mutAAPos,ref_allele_len){
  w <- ceiling(ref_allele_len / 3)
  ret<-as.character(subseq(aaseqnorm,start=mutAAPos,width=w))
  ret
}

getAAvarianttype <- function(aaseqmut,mutAAPos,var_allele_len){
  w <- ceiling(var_allele_len / 3)
  ret<-as.character(subseq(aaseqmut,start=mutAAPos,width=w))
  ret
}
#Function to get variants for frame shifts
getFSAAvarianttype <- function(aaseqmut,mutAAPos,rareport){
  if( class(aaseqmut)=='AAString' ){
    ret <- subseq(aaseqmut,start=mutAAPos,end=rareport)
  }
}


d$aawild[L3crite] <- apply(d[L3crite,c('aaseqnorm','mutAAPos','ref_allele_len')],
                                MARGIN=1,FUN=splat(getAAwildtype))

d$aavar[L3crite] <- apply(d[L3crite,c('aaseqmut','mutAAPos','var_allele_len')],
                           MARGIN=1,FUN=splat(getAAvarianttype))

#Make corrections to Frame shifted varints using rareport position from above.
d$aavar[FScrite] <- apply(d[FScrite,c('aaseqmut','mutAAPos','rareport')], 
                          MARGIN=1,FUN=splat(getFSAAvarianttype))

##########
### 25 ### Get +/- 150 nt from the mutation flanks 
##########
getTranscriptLeftOfMut <- function(leftnt, mutTrnscrtDNAPos, mutTrnscrtDNA){
  mutTrnscrtDNA <- as.character(mutTrnscrtDNA)
  if(mutTrnscrtDNAPos==1){
    ret <- ''
    }
  else{
    e <- mutTrnscrtDNAPos - 1
    ret <- substr(mutTrnscrtDNA, start=leftnt, stop=e)
  }
  ret
} 

getTrnascriptRightOfMut <- function(rightnt, mutTrnscrtDNAPos, var_allele_len, mutTrnscrtDNA){
  mutTrnscrtDNA <- as.character(mutTrnscrtDNA)
  if(mutTrnscrtDNAPos+var_allele_len >= nchar(mutTrnscrtDNA)){
    ret <- ''
    }
  else{
    s <- mutTrnscrtDNAPos+var_allele_len
    ret <- substr(mutTrnscrtDNA, start=s, stop=rightnt )    
  }
  ret
}

#calc left and right ends (check if enough space)
d$leftnt[L3crite] <- ifelse(d$mutTrnscrtDNAPos[L3crite] > 150,
                   unlist(d$mutTrnscrtDNAPos[L3crite]) - 150,
                   1 )
d$rightnt[L3crite] <- 
  ifelse(d$mutTrnscrtDNAlen[L3crite] - (d$mutTrnscrtDNAPos[L3crite] + d$var_allele_len[L3crite] ) > 150,
         d$mutTrnscrtDNAPos[L3crite] + d$var_allele_len[L3crite] +  150,
         d$mutTrnscrtDNAlen[L3crite] )  

#Apply the function to each row of the data frame
#Use splat() instead of spelling out function arguments
d$leftTrnscrptDNA[L3crite] <- apply(d[L3crite,c('leftnt','mutTrnscrtDNAPos','mutTrnscrtDNA')],
                                  MARGIN=1,FUN=splat(getTranscriptLeftOfMut))

d$rightTrnscrptDNA[L3crite] <- apply(d[L3crite,c('rightnt','mutTrnscrtDNAPos','var_allele_len','mutTrnscrtDNA')],
                                    MARGIN=1,FUN=splat(getTrnascriptRightOfMut))



################################################
###  Output Result of    Calculations   #######
##############################################
#Fill in 

#Append level criteria to data frame:
d$L1crite <- L1crite
d$L2crite <- L2crite
d$L3crite <- L3crite
d$FScrite <- FScrite

#Copy data frame for output
d.o <- d

#Fill in Nulls

#Convert Biostrings to regular character strings
d.o$trnscrtDNA       <- sapply(d.o$trnscrtDNA,       as.character)
d.o$mutTrnscrtDNA    <- sapply(d.o$mutTrnscrtDNA,    as.character)
d.o$aaseqnorm        <- sapply(d.o$aaseqnorm,        as.character)
d.o$aaseqmut         <- sapply(d.o$aaseqmut,         as.character)
d.o$aanorm           <- sapply(d.o$aanorm,           as.character)
d.o$aamut            <- sapply(d.o$aamut,            as.character)
d.o$mutaaReport      <- sapply(d.o$mutaaReport,      as.character)
d.o$leftTrnscrptDNA  <- sapply(d.o$leftTrnscrptDNA,  as.character)
d.o$rightTrnscrptDNA <- sapply(d.o$rightTrnscrptDNA, as.character)

#Replace character(0) with ''
replaceCharZero <- function(x){
  if( identical(x,character(0)) ){ret <- '' }
  else{ret <- x}
  ret
}

d.o$trnscrtDNA       <- sapply(d.o$trnscrtDNA,       replaceCharZero)
d.o$mutTrnscrtDNA    <- sapply(d.o$mutTrnscrtDNA,    replaceCharZero)
d.o$aaseqnorm        <- sapply(d.o$aaseqnorm,        replaceCharZero)
d.o$aaseqmut         <- sapply(d.o$aaseqmut,         replaceCharZero)
d.o$aanorm           <- sapply(d.o$aanorm,           replaceCharZero)
d.o$aamut            <- sapply(d.o$aamut,            replaceCharZero)
d.o$mutaaReport      <- sapply(d.o$mutaaReport,      replaceCharZero)
d.o$leftTrnscrptDNA  <- sapply(d.o$leftTrnscrptDNA,  replaceCharZero)
d.o$rightTrnscrptDNA <- sapply(d.o$rightTrnscrptDNA, replaceCharZero)

#Convert Exon Starts & Exon Stops back to a character string
d.o$exonstarts   <- sapply(d.o$exonstarts,   paste, sep='',collapse=' ')
d.o$exonends     <- sapply(d.o$exonends,     paste, sep='',collapse=' ')
d.o$codingStarts <- sapply(d.o$codingStarts, paste, sep='',collapse=' ')
d.o$codingEnds   <- sapply(d.o$codingEnds,   paste, sep='',collapse=' ')
d.o$premutstarts <- sapply(d.o$premutstarts, paste, sep='',collapse=' ')
d.o$premutends   <- sapply(d.o$premutends,   paste, sep='',collapse=' ')

#Make smaller set of unique mutant amino acids 21-mers
dupes <- duplicated(d.o$mutaaReport) #will include one instance of character(0)
d.o.u <- d.o[L3crite & !dupes, ]

#InO
print(paste(sum(dupes[L3crite]),"duplicate mutant report (~ 21mer) peptides"))
print(paste(sum(!dupes[L3crite]),"unique mutant report (~ 21mer) peptides"))


### Text description of columns for header of output data frame ###
coldescripts <- c('transcript'='from UCSC', 'chrom'='from UCSC', 'strand'='from UCSC', 
'txstart'='from UCSC', 'txend'='from UCSC', 'cdsstart'='from UCSC', 'cdsend'='From UCSC',
  'exoncount'='from UCSC', 'exonstarts'='from UCSC', 'exonends'='from UCSC',
  'proteinid'='from UCSC', 'alignid'='from UCSC', 
  'chrom.1'='from given data should match UCSC', 
  'geneSymbol'='geneSymbol data from UCSC kgXref table.',
  'pos'='Position of mutation from source data.',
  'ref_allele'='normal nt from source data',
  'var_allele'='altered nt from source data', 
  'leftflank'='calculated ref genome position one nt to the left of mutation',
  'rightflank'='calculated ref genome position one nt to the right of mutation',
  'noncoding'='boolean flag if mutation is within cdsstart and cdsend',
  'inexon'='boolean flag if mutation is with an exon start and end',
  'exonslen'='calculated number of exons',
  'exCDS'='calculated coding exon starts', 
  'exCDE'='calculated coding exon ends', 
  'codingStarts'='calculated coding starts with modified first coding start',
  'codingEnds'='calculated coding ends with modified last coding start',
  'premutstarts'='calculated coding starts that occur before the mutation',
  'premutends'='calculated coding ends that occurr before the mutation',
  'flagMutSE'='boolean flag for error condition where different number of coding starts and ends',
  'mutTrnscrtDNAPos'='calculated mutation position within the transcript',
  'trnscrtDNA'='normal transcript DNA from UCSC between coding starts and ends',
  'mutAAPos'='calculated amino acid position of mutation',
  'lpmut'='calculated transcript position one nt to the left of mutation start',
  'rpmut'='calculated transcript position one nt to the right of mutation end',
  'mutTrnscrtDNA'='mutated transcript DNA',
  'mutTrnscrtDNAlen'='length of mutant transcript DNA',
  'var_allele_len'='length of variant allele',
  'aaseqnorm'='normal peptide sequence w/o mutation',
  'aaseqmut'='mutant peptide sequence with mutation',
  'isTrunc'='boolean flag for if mutant is immediately truncated',
  'isSynon'='boolean flag for if mutant makes a synonymous mutation',
  'lareport'='calculated amino acid position for left side of ~21 mer to report',
  'lenAAMut'='length of the mutant peptide', 
  'rareport'='calculated amino acid position for right side of ~21 mer to report', 
  'noStop'='boolean flag indicating no stop (*) found in mutant amino acid sequence',
  'mutaaReport'='shortened peptide to +/- of mutation',
  'aawild'='wild type amino acid @ mutation position',
  'aavar'='variant amino acid @ mutation position',
  'leftnt'='transcript postion 150 nt to left of mutation',
  'rightnt'='transcript position 150 nt to right of mutation',
  'leftTrnscrptDNA'='150 nt transcript DNA to left of mutation',
  'rightTrnscrptDNA'='150 nt of transcript DNA to right of mutation',
  'L1crite'='boolean for mutation is in an exon and mutation within coding region',
  'L2crite'='boolean for L1crite AND peptide is not truncated immediately',
  'L3crite'='boolean for L2crite AND mutation is not synonymous',
  'FScrite'='boolean for frame shifts. L3crite AND length var and ref not equal AND diff not evenly divisible by 3'
   )

### Write to File Operations ### 
#Make AAStringSet
peplist <- unlist(d.o$mutaaReport[L3crite])
peplist <- gsub(pattern='\\*$',replacement='', x=peplist) #remove * char at end
aaOutputSet <- AAStringSet(peplist)
names(aaOutputSet)<-paste(d.o$transcript[L3crite],"leftflank",d.o$leftflank[L3crite])
#Make Uniques AAStringSet
peplist <- unlist(d.o.u$mutaaReport)
peplist <- gsub(pattern='\\*',replacement='', x=peplist) #remove * characters
aaUniqueSet <- AAStringSet(peplist)
names(aaUniqueSet)<-paste(d.o.u$transcript,"leftflank",d.o.u$leftflank)

#write FASTA format output
outname <- 'shortMutantPeptides_FULL.fa'
outfile <- paste(myOutDir,outprfx,outname,sep='')
write.XStringSet(aaOutputSet, filepath=outfile, append=FALSE, format="fasta")

#write FASTA format output
outname <- 'shortMutantPeptides_Unique.fa'
outfile <- paste(myOutDir,outprfx,outname,sep='')
write.XStringSet(aaUniqueSet, filepath=outfile, append=FALSE, format="fasta")


### Write All Data for Full Set. ###
cls<- unlist( lapply(d.o,class), use.names=FALSE )
df.out<- d.o
#open connection
outfile<-file(description=paste(myOutDir,outprfx,'mutsToPepsAllData_FULL.txt',sep=''),
              open='w', encoding='UTF-8', raw=FALSE)
#write data
write.table(x=df.out, file=outfile,append=FALSE,quote=FALSE,sep='\t',eol='\n',na='NA',
            dec='.', row.names=FALSE, col.names=TRUE, qmethod='escape')
#close connection
close(outfile)

### Write All Data for Uniques Set. ###
df.out<- d.o.u
#open connection
outfile<-file(description=paste(myOutDir,outprfx,'mutsToPepsAllData_Unique.txt',sep=''),
              open='w', encoding='UTF-8', raw=FALSE)
#write data
write.table(x=df.out, file=outfile,append=FALSE,quote=FALSE,sep='\t',eol='\n',na='NA',
            dec='.', row.names=FALSE, col.names=TRUE, qmethod='escape')
#close connection
close(outfile)

### Write Column definitions for All Data ###
df.out <- as.data.frame(coldescripts)
#open connection
outfile<-file(description=paste(myOutDir,outprfx,'mutsToPepsAllData_ColumnDefs.txt',sep=''),
              open='w', encoding='UTF-8', raw=FALSE)
#write data
write.table(x=df.out, file=outfile,append=FALSE,quote=FALSE,sep='\t',eol='\n',na='NA',
            dec='.', row.names=TRUE, col.names=TRUE, qmethod='escape')
#close connection
close(outfile)


#InO
rtt <- Sys.time() - analysis.start
print("Run Time Total:")
print(rtt)
print("END OF LINE")
sink()

################################################################################
#        Clean Up                                                       #######
##############################################################################

### REMOVE ALL OBJECTS ### 
rm(list=ls(all=TRUE))

### Detach Packages (reverse order from load) ###
# detach("package:reshape")
# detach("package:stringr")
# detach("package:plyr")
# detach("package:GenomicFeatures")
# detach("package:BSgenome.Hsapiens.UCSC.hg19")
# detach("package:BSgenome")
# detach("package:Biostrings")



################################################################################
#        Scrap Code                                                     #######
##############################################################################

# arrayToDigitString <- function(x){
#   ret <- paste(x, sep='',collapse=' ')
#   ret
# }


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
#                       MARGIN=1,FUN=splat(modifyCodingEnds))








