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
doRowVerbose <- function(transcript,chr,leftflank,rightflank,
                         ref_allele,var_allele,
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
  
  if(leftflank < cdsstart || leftflank > cdsend){
    print(paste("Mutation outside of coding region for:",transcript,
                "left flank",leftflank,sep=' '))
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
  trnscrtDNAset <- tryCatch({
    DNAStringSet(unmasked(Hsapiens[[chrom]]), 
                                start=as.vector(codingStarts), 
                                end=as.vector(codingEnds) );
    }, error= function(ex){
      cat("Genome Lookup Error!")
      return(NA)
    }
  )
  if(class(trnscrtDNAset) != "DNAStringSet")
  {
    cat(" Returning NA.\n")
    return(NA)
  }
  

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
  
  #calculate mutation position in amino acid sequence.
  #also determine if mutation is at a codon start (used later)
  if(strand=='+'){
    mutAAPosition <- ceiling(mutTrnscrtDNAPosition / 3) 
    #if remainder is zero and ref_allele == '' (insertion) add 1 to AA position
    #Only apply for INDELs
    if(mutTrnscrtDNAPosition %% 3 == 0 && nchar(var_allele) != nchar(ref_allele)){
        mutAAPosition <- mutAAPosition +1
    }
  }else{
    mutAAPosition <- (length(trnscrtDNA)/3) - (ceiling(mutTrnscrtDNAPosition/3)) +1
    if((length(trnscrtDNA)-mutTrnscrtDNAPosition+1) %% 3 == 0  && nchar(var_allele) != nchar(ref_allele) ){
        mutAAPosition <- mutAAPosition +1
    }
  }
  
  
  #debug here
  ###Make mutation and get mutant DNA!
  #Concatenation of [left side],[variant allele],[right side]
  lpmut <- mutTrnscrtDNAPosition - 1
  rpmut <- mutTrnscrtDNAPosition + nchar(ref_allele)
  if(lpmut == 0){
    print("Mutation includes start/end of coding sequence. Handle Later.")
    return(NA)
  }
  
  
  mutantTrnscrptDNA <- tryCatch( 
    expr= { 
      xscat(
        substr(trnscrtDNA, 1, lpmut), 
        DNAString(var_allele), 
        substr(trnscrtDNA,rpmut, length(trnscrtDNA))
        );
    }, error = function(ex) {
      cat("Holy crap an error\n")
      return(NA)
    }
  )
  
  if(class(mutantTrnscrptDNA) != 'DNAString' ){
    cat("mutantTranscptProblem. class: ",
        class(mutantTrnscrptDNA),
        " Val: ",
        as.character(mutantTrnscrptDNA),
        "\n")
    return(NA) #return from function NA
  }

  
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
  #aaSequence <- translate(codingMRNA)
  aaSequence <-tryCatch({
    translate(codingMRNA)
  }, error = function(ex) {
    cat("Normal Translation Error in:")
    print(codingMRNA)
    print(ex)
    return(NA)
  }
  )
  
  if(class(aaSequence) != "AAString"){
    cat("AA error Returning NA\n")
    return(NA)
  }
  
  aaMutSequence <-tryCatch({
    translate(mutantMRNA)
    }, error = function(ex) {
      cat("Mutant Translation Error in: ")
      print(mutantMRNA)
      print(ex)
      return(NA)
    }
  )
  
  if(class(aaMutSequence) != "AAString"){
    cat("AA error Returning NA\n")
    return(NA)
  }
  
  #Check if variant amino acid is a stop.
  if(aaMutSequence[mutAAPosition] == AAString("*")){
    print(paste("Mutation truncates peptide.",transcript,leftflank))
    return(NA)
  }
  
  #Check Mutant Amino Acid against reference (from mutAA to end of peptide)
  if(aaMutSequence[mutAAPosition:length(aaMutSequence)] == aaSequence[mutAAPosition:length(aaSequence)]){
    print(paste("Synonymous mutation.",
                mutAAPosition, 
                length(trnscrtDNA),
                mutTrnscrtDNAPosition,
                as.character(aaMutSequence[(mutAAPosition):(mutAAPosition)]),
                as.character(aaSequence[(mutAAPosition):(mutAAPosition)]),
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
  
  #Check for frame shift
  frameShift <- FALSE #false by default
  if( nchar(var_allele) != nchar(ref_allele) ){
    #check if sum of ref_alleles and var_alleles not divisible by 3 evenly
    if( (nchar(var_allele)+nchar(ref_allele)) %%3 !=0 ){
      frameShift <- TRUE
    }
  }
  
  #In the case of a frame shift, find stop codon or end of aa
  if(frameShift){
    stoppos <- regexpr('\\*',aaMutSequence) #find location of *
    if(stoppos == -1){
      #no stop codon?
      #set truncation right side to end of amino
      ra <- length(aaMutSequence)
    }
    else{
      #set truncation right side to 1 before first stop
      ra <- stoppos[1] - 1
    }
  }
  
  #Do truncation
  aaMutTrunc <- tryCatch({
     subseq(aaMutSequence,start=la,end=ra);
  }, error = function(ex){
    print(ex)
    cat("AminoAcid Trunc Error\n")
    
    return(NA)
  })
  #Debug truncation
  aaNormTrunc <- tryCatch({
    subseq(aaSequence,start=la,end=ra);
  }, error = function(ex){
    print(ex)
    cat("AminoAcid Trunc Error\n")
    
    return(NA)
  })
  
  if(rettype == 'rowdata'){
    #Return a row of data for rbinding
    retRow <- list(
      'transcript'=transcript,
      'leftflank'=leftflank,
      'mutTrnscPos'=mutTrnscrtDNAPosition,
      'trnscStrand'=strand,
      'lengthTranscript'=length(trnscrtDNA),
      'mutAAPos'=mutAAPosition,
      'refAA'=as.character(aaSequence[mutAAPosition]),
      'mutAA'=as.character(aaMutSequence[mutAAPosition]),
      'norm_mer'=as.character(aaNormTrunc),
      'drunkmer'=as.character(aaMutTrunc),
      'trnscrptDNA'=as.character(trnscrtDNA),
      'seq'=as.character(seq)
      )
    return(retRow)
  }
  
  #Return truncated mutant AA sequence
  return(as.character(aaMutTrunc))

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
  '2219_mutsrefs2012-05-08_183106.txt',
  '2221_mutsrefs2012-05-08_183335.txt'
  )

#import combined data frame
infile<-paste(impdir,infilenames[2],sep='/')
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



################################################
###  Output Intermediate Calculations   #######
##############################################
# cols.out<-c(colnames(dfc)[1:8],colnames(dfc)[14:17])
# df.out<- dfc[,cols.out]
# 
# outfile<-file(description=paste(myOutDir,'GastricIntermediateCalcs.txt',sep=''),
#               open='w', encoding='UTF-8', raw=FALSE)
# #write data
# write.table(x=df.out, file=outfile,append=FALSE,quote=FALSE,sep='\t',eol='\n',na='NA',
#             dec='.', row.names=FALSE, col.names=TRUE, qmethod='escape')
# #close connection
# close(outfile)

################################################
###  Analysis                           #######
##############################################

# #Small test
# dfsm <- dfc[c(1:20),]
# dfsm$aaS[1:20] <- apply(dfc[1:20,], MARGIN=1, FUN=function(x){splat(doRow)(x)})
# 
# ### Small test RUN w/ VERBOSE ###
# #get list of character vectors
# retl <- apply(dfc[c(1:20),], MARGIN=1, FUN=function(x){splat(doRowVerbose)(x)})
# #Convert returned list to data frame
# retdf <- do.call('rbind',lapply(retl, "["))


#Start Sink
sink(paste(myOutDir,'AnalysisSink.txt',sep=''))
### Full Run ###
dfbig <- dfc[,-which(colnames(dfc)=='seq')]
dfbig$aaS <- apply(dfc, MARGIN=1, FUN=function(x){splat(doRow)(x)})
#End Sink
sink()

### Remove NAs ###
dfbig.i<-dfbig[!is.na(dfbig$aaS),]
### Make Unique Peptides Set ###
dupes <- duplicated(dfbig.i$aaS)
dfbig.u <- dfbig.i[!dupes,]
################################################
###  Output Data                        #######
##############################################


#Make AAStringSet
aaOutputSet <- AAStringSet(dfbig.i$aaS)
names(aaOutputSet)<-paste(dfbig.i$transcript,"leftflank",dfbig.i$leftflank)
#Make Uniques AAStringSet
aaUniqueSet <- AAStringSet(dfbig.u$aaS)
names(aaUniqueSet)<-paste(dfbig.u$transcript,"leftflank",dfbig.u$leftflank)

#write FASTA format output
outname <- 'mutantPeps_FULL.fa'
outfile <- paste(myOutDir,outname,sep='')
write.XStringSet(aaOutputSet, filepath=outfile, append=FALSE, format="fasta")

#write FASTA format output
outname <- 'mutantPeps_Unique.fa'
outfile <- paste(myOutDir,outname,sep='')
write.XStringSet(aaUniqueSet, filepath=outfile, append=FALSE, format="fasta")


### Output processed data frames ###
cols.out<-c(colnames(dfbig.i)[1:8],colnames(dfbig.i)[13:19])
df.out<- dfbig.i[,cols.out]
#open connection
outfile<-file(description=paste(myOutDir,'MutsToPeps_FULL.txt',sep=''),
                 open='w', encoding='UTF-8', raw=FALSE)
#write data
write.table(x=df.out, file=outfile,append=FALSE,quote=FALSE,sep='\t',eol='\n',na='NA',
            dec='.', row.names=FALSE, col.names=TRUE, qmethod='escape')
#close connection
close(outfile)


#select columns
cols.out<-c(colnames(dfbig.u)[1:8],colnames(dfbig.u)[13:19])
df.out<- dfbig.u[,cols.out]
#open connection
outfile<-file(description=paste(myOutDir,'MutsToPeps_Unique.txt',sep=''),
              open='w', encoding='UTF-8', raw=FALSE)
#write data
write.table(x=df.out, file=outfile,append=FALSE,quote=FALSE,sep='\t',eol='\n',na='NA',
            dec='.', row.names=FALSE, col.names=TRUE, qmethod='escape')
#close connection
close(outfile)


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













