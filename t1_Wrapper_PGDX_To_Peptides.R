################################################################################
# Title:  Wrapper script for PGDX to Mutation Data & Peptides                  #
# Date:   2012-07-10                                                           #
# Author: Colin A. Gross                                                       #
# Desc:   SnpEff output to text from VCF. Output mutated peptides              #
################################################################################




################################################################################
#        SETUP                                                                 #
################################################################################

#input processing script file
inscriptfile <- 'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/scripts/t1_SnpEff_To_Peptides_wArgs.R'

#standard input mutations to peptides processing script file
stdscriptfile <- 'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/scripts/t3_Std_Mutations_To_Peptides_call.R'

#input Directory
impdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata'

infilenames <- c(
  '2369/mutsrefs.txt'

  )

#Set file prefixes for output
outprefixes <- c(
  '2369'
  )

#output directory to be used by script
outdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata/'

starttime <- Sys.time()
print(starttime)

for(i in 1:length(infilenames)){
  cat("Run ",i," started... ")
  
  #infile to be used by script
  infile<-paste(impdir,infilenames[i],sep='/')
  
  #prefix to prepend to all output files
  outprefix<-paste(outprefixes[i],'_',sep='')
  
  #run script to process input file (dfc should be added to environment)
  source(inscriptfile, local=TRUE, echo=FALSE, verbose=FALSE)
  
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


