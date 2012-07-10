################################################################################
# Title:  PGDX input --> standard input                                        #
# Date:   2012-07-09                                                           #
# Author: Colin A. Gross                                                       #
# Desc:   take PGDX input frame from args and leave standard input frame for   #
#          use by calling (parent) environment                                 #
################################################################################


################################################################################
#        SETUP                                                                 #
################################################################################
require(stringr)

################################################################################
#        PARSE ARGUMENTS OR USE DEFAULTS                                       #
################################################################################
str_rev <- function(x){
  paste(rev(substring(x,1:nchar(x),1:nchar(x)) ),collapse="")
}

argpattern <- '(.*?=\\w*?)\\s'
#get the arguements as a single string
rawargline <- (paste(commandArgs(trailingOnly=TRUE),collapse=' '))
#clean raw input
clnargline <- str_replace_all(rawargline,'"','') #replace quotes with nothing
clnargline <- str_replace_all(clnargline,'\\\\','/') #replace \ with /
clnargline <- str_pad(clnargline,nchar(clnargline)+1,side='left',pad=' ')
#paste in a dummy arg to address the case of no arguments
clnargline <- paste(clnargline, 'dum=dum',sep=' ')
#reverse cleaned input
revargline <- str_rev(clnargline)
#extract lav=rav matches (reversed var=val)
revarglist <- str_extract_all(revargline,argpattern)[[1]]
#trim and re-revers
charargs <- unlist(lapply(revarglist,function(x) str_trim(str_rev(x)) ))
#split var and values
listargs<-str_split(charargs,'=')
#make data frame of passed in var names and values
argdf<-data.frame(var=mapply(unlist,listargs)[1,],
                  val=mapply(unlist,listargs)[2,])

#Check if required variables are in the parsed command line arguments
if('infile' %in% argdf$var){ 
  infile <- argdf$val[argdf$var=='infile']
}
if('outdir' %in% argdf$var){
  outdirr <- argdf$val[argdf$var=='outdir']
}
if('outprefix' %in% argdf$var){
  outprefix <- argdf$val[argdf$var=='outprefix']
}

#Check if required variables exist
if( !exists('infile')){ stop('No infile specified.') } #Die
if( !exists('outdir')){ stop('No outdir specified.') } #Die
if( !exists('outprefix')){ 
  outprefix <- ''
}


################################################
###  Import Data                        #######
##############################################

#infile <- 'S:/TIL-LAB/Staff/Colin/Projects/MutationsToPeptides/procdata/2369/mutsrefs.txt'
dfc <-read.table(infile, header=TRUE, sep="\t",
                 comment.char="#",encoding="UTF-8",stringsAsFactors=FALSE)


################################################
###  Check Column Names                 #######
##############################################

expectedcols <- c('Transcript','leftflank.1','rightflank.1','Ref_allele',
                  'Var_allele','chrom','strand','txStart','txEnd','cdsStart',
                  'cdsEnd','exonCount','exonStarts','exonEnds','proteinID',
                  'alignID','geneSymbol')

if(length(expectedcols) - sum(expectedcols %in% colnames(dfc)) != 0){
  cat("Missing input columns")
}else{
  cat("All input columns present")
}


################################################
###  Process Input                      #######
##############################################

#Calculate back to left flank and right flank positions
dfc$leftflank <- dfc$leftflank.1 - 1
dfc$rightflank <- dfc$rightflank.1 + 1

#Drop old columns
dfc$leftflank.1 <- NULL
dfc$rightflank.1 <- NULL

#Convert all column names to lower case
colnames(dfc) <- tolower(colnames(dfc))


################################################
###  Cleanup Variables                  #######
##############################################

remove('infile','expectedcols')
#leave dfc for calling parent environment.

