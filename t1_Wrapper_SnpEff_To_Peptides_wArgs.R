################################################################################
# Title:  Wrapper script for SNPEff to Mutation Data & Peptides                #
# Date:   2012-05-06                                                           #
# Author: Colin A. Gross                                                       #
# Desc:   SnpEff output to text from VCF. Output mutated peptides              #
################################################################################




################################################################################
#        SETUP                                                                 #
################################################################################

#script file
scriptfile <- 'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/scripts/t1_SnpEff_To_Peptides_wArgs.R'

#input Directory
impdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata'

infilenames <- c(
  '2219_mutsrefs_2012-06-04.txt',
  '2221_mutsrefs_2012-06-04.txt',
  '2246_mutsrefs_2012-06-04.txt',
  '2359_mutsrefs_2012-06-04.txt',
  '2556_mutsrefs_2012-06-04.txt',
  '3466_mutsrefs_2012-06-04.txt',
  'gastric2_mutsrefs_2012-06-04.txt'
  )

outprefixes <- c('2219','2221','2246','2359','2556','3466','gastric')


#output directory to be used by script
outdir<-'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata/'

for(i in 1:length(infilenames)){
  cat("Run ",i," started... ")
  
  #infile to be used by script
  infile<-paste(impdir,infilenames[i],sep='/')
  
  #prefix to prepend to all output files
  outprefix<-paste(outprefixes[i],'_',sep='')
  
  source(scriptfile, local=TRUE, echo=FALSE, verbose=FALSE)
  
  #Remove objects resulting from script run.
  objs <- ls(all=TRUE)
  exempt <- c('i','impdir','outprefixes','outdir','infile','outprefix','scriptfile','infilenames')
  objs <- objs[!(objs %in% exempt)] #exclude exempt variables from removal
  rm(list=objs)
  
  cat("Run ",i," complete.\n")
}





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


