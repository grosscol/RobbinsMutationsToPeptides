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
# NOTE: UCSC database stores start positions using 0-based counting
#        They store end positions using 1-based index.
getMrnaPos <- function(p,exS,exE){
  #calculate array of exon lengths
  exonls <- exE - exS
  
  #find out which exon the mutation occurs in.
  #highest start positon that is less than mutation positon
  exN <- which(exS == max(exS[p >= exS]) )
  
  #print(c(exN,p,exS))
  #Sum the lengths of the prior exons along with length to mutation from exN.
  #This will be the length into to cDNA where the mutation position is.
  cdnal <- sum(exonls[0:(exN-1)]) + (p - exS[exN])
  
  #Since mRna is reversed from genomic dna, get 1-based position into mRna by:
  # (length of exons) - (index into exon)
  #Careful with length calcs. The above have been corrected for:
  #The UCSC database start postions are 0 based
  #The UCSC database end positions are 1 based
  #Genome browser is 1 based index
  cdnal <- sum(exonls) - cdnal + 1
  
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


#Function wrapper designed for splat with df2 & verbose
doRowVerbose <- function(transcript,chr,leftflank,rightflank,var_allele,
                         chrom,strand,txstart,txend,cdsstart,cdsend,
                         exoncount,exonstarts,exonends,proteinid,alignid,seq,
                         rettype='rowdata', ...){
  
  crap <- list(...) #make a named list of the crap we are ignoring
  doRow(transcript,chr,leftflank,rightflank,
        ref_allele, var_allele,
        chrom,strand,txstart,txend,cdsstart,cdsend,
        exoncount,exonstarts,exonends,proteinid,alignid,seq,
        rettype='rowdata')
  
}

#Function designed for splat with df2
doRow <- function(transcript,chr,leftflank,rightflank,
                  ref_allele, var_allele,
                  chrom,strand,txstart,txend,cdsstart,cdsend,
                  exoncount,exonstarts,exonends,proteinid,alignid,seq,
                  rettype='mutaa', ...){
  #Ignore unused variables
  crap <- list(...) #make a named list of the crap we are ignoring
  
  #Check return type
  if( rettype %in% c('mutaa','rowdata') ){
    #Valid return type specified
  }else{
    print('rettype not \"mutaa\" or \"rowdata\"')
    return(NA)
  }
  
  #only deal with point mutations for now
  if( (rightflank - leftflank - 1) > nchar(var_allele) ){
    print("Site longer than variant. Deletion. Deal with later.")
    return(NA)
  }else if((rightflank - leftflank - 1) < nchar(var_allele)){
    print("Site shorter than variant. Insertion. Deal with later.")
    return(NA)
  }
  
  #Passed in sequence will be a char string. Convert to DNAString for
  # internal use.
  seq <- DNAString(seq)
  #get mRNA length of data from knownGeneMRna table seq column in hg18 database.
  mrnalen <- length(seq)
  
  #calculate sum of exon lengths
  exonlen <- calcExonLength(exonstarts,exonends)
  
  #get coding start exon
  #highest start positon that is less than cdsstart positon
  exCDS <- which(exonstarts == max(exonstarts[cdsstart >= exonstarts]) )
  
  #get coding end exon
  #lowest end positon that is greater than cdsend positon
  exCDE <- which(exonends == min(exonends[cdsend <= exonends]) )
  
  #Get starts and stops for coding exons
  codingExons  <- c(exCDS:exCDE)
  codingStarts <- exonstarts[codingExons] + 1 #starts are 0 based
  codingEnds   <- exonends[codingExons]
  
  #replace first exon coding start position with cdsstart
  codingStarts[1] <- cdsstart + 1 #starts are 0-based
  #replace last exon coding end position with cdsend
  codingEnds[length(codingEnds)] <-cdsend
  
  #This takes a fucking while. Get all coding exon dna in a set
  trnscrtDNAset <- DNAStringSet(unmasked(Hsapiens[[chrom]]), 
                                start=as.vector(codingStarts), 
                                end=as.vector(codingEnds) )
  #Combine all sequences together.
  trnscrtDNA <- unlist(trnscrtDNAset)
  
  
  #calculate mutation start position (0-based ? or 1-based?) assuming 1-based
  mutStartPos <- leftflank + 1
  
  #get all the start and end positions that occur before the mutation position
  #length = mutStartPos - exonstarts[exMUT] + sum(exonEnds[])
  mutLengthStarts <- codingStarts[codingStarts <= mutStartPos]
  mutLengthEnds <- codingEnds[codingEnds < mutStartPos]
  #append the mutation position to the list of coding ends that occured 
  #before the mutation position.
  mutLengthEnds[length(mutLengthEnds)+1] <- mutStartPos
  
  #a mismatch in the number of starts and ends indicates that the mutation
  #position is not in a coding region.
  if(length(mutLengthEnds) != length(mutLengthStarts)){
    print(paste("Mutation at non-coding site for:",transcript,
                "left flank",leftflank,sep=' '))
    return(NA)
  }
  
  #The position in the trnscrtDNA string = length into the string.
  #length calculations have to accomodated 1 based start (so just subtract 1)
  mutTrnscrtDNAPosition <- sum(mutLengthEnds - (mutLengthStarts - 1 ))
  
  #calculate mutation position in amino acid sequence
  if(strand=='+'){
    mutAAPosition <- ceiling(mutTrnscrtDNAPosition / 3)
  }else{
    mutAAPosition <- floor((length(trnscrtDNA)-mutTrnscrtDNAPosition) / 3 ) + 1
  }
  
  
  ###Make mutation and get mutant DNA!
  #Concatenation of [left side],[variant allele],[right side]
  lpmut <- mutTrnscrtDNAPosition - 1
  rpmut <- mutTrnscrtDNAPosition + nchar(ref_allele)
  mutantTrnscrptDNA <- xscat(
                          substr(trnscrtDNA, 1, lpmut), 
                          DNAString(var_allele), 
                          substr(trnscrtDNA,rpmut, length(trnscrtDNA))
                       )
  #if length of mutation site == length of variant at site, then substitution
#   if((rightflank - leftflank - 1) == nchar(var_allele)){
#     
#     
#     xscat(trnscrt)
#     #make a copy
#     mutantTrnscrptDNA <- trnscrtDNA
#     #make the mutation
#     subseq(mutantTrnscrptDNA, start=mutTrnscrtDNAPosition, 
#            width=nchar(var_allele)) <- DNAString(var_allele)
#     
#   }
  
  #Get coding MRNA and coding Mutant MRNA
  #genome is stored as '-' sense sequence. reverseComplement this strand.
  # if strand is '+' simply do rna substitution.
  if(strand == '+'){
    codingMRNA <- dna2rna(trnscrtDNA)
    mutantMRNA <- dna2rna(mutantTrnscrptDNA)
  }else{
    codingMRNA <- dna2rna(reverseComplement(trnscrtDNA))
    mutantMRNA <- dna2rna(reverseComplement(mutantTrnscrptDNA))
  }
  
  #Do translation
  aaSequence <- translate(codingMRNA)
  aaMutSequence <- translate(mutantMRNA)
  
  #Check if variant amino acid is a stop.
  if(aaMutSequence[mutAAPosition] == AAString("*")){
    print(paste("Mutation truncates peptide.",transcript,leftflank))
    #Subset the mutant AA sequence at the stop codon
    #aaMutSequence <- subseq(aaMutSequence,start=1,end=mutAAPosition)
    return(NA)
  }
  
  #Check Mutant Amino Acid against reference
  if(aaMutSequence[mutAAPosition] == aaSequence[mutAAPosition]){
    print(paste("Synonymous mutation.",
                mutAAPosition, 
                length(trnscrtDNA),
                mutTrnscrtDNAPosition,
                as.character(aaMutSequence[mutAAPosition]),
                as.character(aaSequence[mutAAPosition]),
                as.character(trnscrtDNA[1:6]),
                as.character(codingMRNA[(length(codingMRNA)-6):(length(codingMRNA))])
                ))
    return(NA)
  }
  
  #Calc 21-mer trunc points for AAString to 10 upstream & downstream
  la<-1
  ra<-length(aaMutSequence)
  #Check that there are at least 10 AA to the left
  if(mutAAPosition >= 10){ 
    la <- mutAAPosition-10
  }
  #Check that there are at least 10 AA to the right
  if( (length(aaMutSequence) - mutAAPosition) > 10) {
    ra <- mutAAPosition + 10
  }
  
  #Do truncation
  aaMutTrunc <- subseq(aaMutSequence,start=la,end=ra)
  
  
  #Sanity Checks
  #   mrnalen - exonlen # should be 0
  #   txstart - exonstarts[[1]] # should be 0
  #   txend - exonends[[exoncount]] # should be 0
  #   length(trnscrtDNA) - length(seq) #should be zero or negative
  #   length(trnscrtDNA) / 3 #should be an integer
  #   #Report reference alleles
  #    paste(transcript, leftflank, "refs (NT AA):",
  #          as.character(subseq(trnscrtDNA,start=mutTrnscrtDNAPosition,width=1)),
  #          as.character(subseq(aaSequence,start=mutAAPosition,width=1)), sep=' ')
  #   length(aaMutTrunc) #should be 21 or less
  
  
  if(rettype == 'rowdata'){
    #Return a row of data for rbinding
    retRow <- list(
      'transcript'=transcript,
      'leftflank'=leftflank,
      'mutTrnscPos'=mutTrnscrtDNAPosition,
      'trnscStrand'=strand,
      length(trnscrtDNA),
      'mutAAPos'=mutAAPosition,
      'refAA'=as.character(aaSequence[mutAAPosition]),
      'mutAA'=as.character(aaMutSequence[mutAAPosition]),
      'drunkmer'=as.character(aaMutTrunc)
      )
    return(retRow)
  }
  #Return truncated mutant AA sequence
  return(aaMutTrunc)
}


################################################
###  Import Data                        #######
##############################################

###########################
# Transcript + Mutations  #
###########################

#Import Directory
impdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata/'

#import combined data frame
infile<-paste(impdir,'mutsrefs2012-04-30_141449.txt',sep='')
cls<- c( rep("character",3), rep("integer",5), rep("character",6), integer, rep("character",2) )
dfc <-read.table(infile, header=TRUE, sep="\t",  colClasses=cls,
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
#For point mutations (given ref allele is not *)
pnts <- dfc$ref_allele != '*'
dfc$leftflank[pnts]  <- dfc$pos[pnts] - 1
dfc$rightflank[pnts] <- dfc$pos[pnts] + 1

#For deletions
#Find var alleles with '-' as these are the way deletions are marked in this data
dels <- grep(pattern='\\-', x=dfc$var_allele)
#copy var allele data to ref allele (where it should have been)
dfc$ref_allele[dels] <- substr(dfc$var_allele[dels],2,nchar(dfc$var_allele[dels]))
dfc$leftflank[dels] <- dfc$pos[dels] - 1
dfc$rightflank[dels] <- dfc$pos[dels] + nchar(dfc$ref_allele[dels])
#Set var_alleles to blank... since they are supposed to be deletions
dfc$var_allele[dels] <- ''

#For insertions (do so that first nuc of var sequence occupies position pos)
ins <- grep(pattern='\\+',x=dfc$var_allele)
dfc$leftflank[ins] =  dfc$pos[ins] - 1 
dfc$rightflank[ins] =  dfc$pos[ins]
#remove '+' from the var allele
dfc$var_allele[ins] <- substr(dfc$var_allele[ins],2,nchar(dfc$var_allele[ins]))
#Set ref_allele to blank... since it is supposed to be blank
dfc$ref_allele[ins] <- ''

################################################
###  Analysis                           #######
##############################################

#Small test
dfsm <- dfc[c(1:20),]
dfsm$aaS[1:20] <- apply(dfc[1:20,], MARGIN=1, FUN=function(x){splat(doRow)(x)})


### Full Run ###
dfbig <- dfc[,c(1:16)]
dfbig$aaS <- apply(dfc, MARGIN=1, FUN=function(x){splat(doRow)(x)})

################################################
###  Output Data                        #######
##############################################

### Output sorted and cleaned up data
#reorder data frame so that columns of like classes are next to each other
#character & Date, factor,logical, integer, numeric columns  
# cls<- unlist( lapply(dfc,class), use.names=FALSE )
# ord <- c( grep('character',cls), grep('Date',cls), grep('factor',cls), grep('logical',cls),
#           grep('integer',cls),grep('numeric',cls) )
# dfc.o<-dfc[,ord]
# rle(cls[ord])

















