################################################################################
# Title:  Wrapper script for Standard input to Mutation Data & Peptides        #
# Date:   2012-08-03                                                           #
# Author: Colin A. Gross                                                       #
# Desc:   Assuming that the standard input for the mutations to peptides alg   #
#   has already been obtained through other means, this just calls a list of   #
#   input files throught the mutations to peptides functions.                  #
################################################################################




################################################################################
#        SETUP                                                                 #
################################################################################
library("plyr")
#Input should already be std input columns
#input processing script file
#inscriptfile <- 'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/scripts/t2_PGDX_to_stdinput.R'

#standard input mutations to peptides processing script file
stdscriptfile <- 'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/scripts/t3_Std_Mutations_To_Peptides_call.R'

#input Directory
impdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata'

infilenames <- c(
  '582T_pgdx/PGDX582T_SomMuts_3-22-13_mutsrefs.txt',
  '579T_pgdx/PGDX579T_SomMuts_3-21-13_mutsrefs.txt',
  '580T_pgdx/PGDX580T_SomMuts_3-21-13_mutsrefs.txt',
  '581T_pgdx/PGDX581T_SomMuts_3-22-13_mutsrefs.txt')
  

#Set file prefixes for output
outprefixes <- c(
  '582T_pgdx_refseq',
  '579T_pgdx_refseq',
  '580T_pgdx_refseq',
  '581T_pgdx_refseq'
  )


#output directory to be used by script
outdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata/'

#Load correct UCSC Library for the input. Do NOT load both.
#library("BSgenome.Hsapiens.UCSC.hg18")
library("BSgenome.Hsapiens.UCSC.hg19")

starttime <- Sys.time()
print(starttime)

for(i in 1:length(infilenames)){
  cat("Run ",i," started... ")
  
  #infile to be used by script
  infile<-paste(impdir,infilenames[i],sep='/')
  
  #prefix to prepend to all output files
  outprefix<-paste(outprefixes[i],'_',sep='')
  
  #run script to process input file (dfc should be added to environment)
  #source(inscriptfile, local=TRUE, echo=FALSE, verbose=FALSE)
  
  #simply read infile as data frame into expected data frame name dfc
  dfc <-read.table(infile, header=TRUE, sep="\t",
                   comment.char="#",encoding="UTF-8",stringsAsFactors=FALSE)
  colnames(dfc)<-tolower(colnames(dfc)) #needs lower case col names
  dfc <- plyr::rename(dfc, c(name="transcript")) #name needs to be changed to transcript
  
  #run script to process standard input data frame (dfc)
  source(stdscriptfile, local=TRUE, echo=FALSE, verbose=FALSE)
  
  #Remove objects resulting from script run.
  objs <- ls(all=TRUE)
  exempt <- c('starttime','i','impdir','outprefixes','outdir','infile',
              'outprefix','inscriptfile','infilenames','stdscriptfile')
  objs <- objs[!(objs %in% exempt)] #exclude exempt variables from removal
  rm(list=objs)
  
  cat("Run ",i," complete:",Sys.time(),"\n")
}

runduration <- Sys.time() - starttime
print("Full duration: ")
print(runduration)


####################
### Cleanup     ###
##################

#Check which database was loaded, and detach it.
dbl <- c("package:BSgenome.Hsapiens.UCSC.hg18", 
         "package:BSgenome.Hsapiens.UCSC.hg19")

loaded <- unlist(sapply(dbl, FUN=function(v) grep(v,search(),value=TRUE)) )
#detach loaded packages
sapply(loaded, FUN=detach, character.only=TRUE)

### Detach Packages (reverse order from load) ###
detach("package:reshape")
detach("package:stringr")
detach("package:plyr")
detach("package:GenomicFeatures")
detach("package:BSgenome")
detach("package:Biostrings")

### REMOVE ALL OBJECTS ###
rm(list=ls(all=TRUE))



