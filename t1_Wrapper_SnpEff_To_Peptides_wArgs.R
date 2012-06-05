################################################################################
# Title:  Wrapper script for SNPEff to Mutation Data & Peptides                #
# Date:   2012-05-06                                                           #
# Author: Colin A. Gross                                                       #
# Desc:   SnpEff output to text from VCF. Output mutated peptides              #
################################################################################




################################################################################
#        SETUP                                                                 #
################################################################################










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

#import combined data frame
infile<-paste(impdir,infilenames[7],sep='/')



