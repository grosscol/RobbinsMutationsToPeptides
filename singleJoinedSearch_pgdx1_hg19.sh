#!/bin/bash

###############################################################################
#   Title: Single Sample Mutations Search cDNA in UCSC hg19
#   Author: Colin A. Gross
#   Date: 2012-04-27
#   Desc:
#
# This script is designed to take a list of transcript names and get the cDNA
# sequences from the UCSC database mirror hosted on biobase.nih.gov.
# Expected columns are:
#  Chrom Position  Reference  Change....
#
###############################################################################

############################
#   Function Definitions   #
############################
function throwError
{
  echo ERROR: $1
  exit 1
}


############################
#   Variables              #
############################
infile=     #name of arg 1, the input data file
sqlcmnd='mysql --[CONNECTION DETAILS]'
datetime=`date +%F_%H%M%S`
outdir='./'
outfile='../procdata/'"mutsrefs_$datetime"'.txt'


#####################################
######   Main                 ######
###################################

#################
### CHECK INPUT ##
###################
#Check length of arguements

if [ $# -ne 2 ]
then 
  echo "usage: $0 inputfile outputdir"
  exit 1   # status 1 signals an error
fi

#Check that input file (arg 1) exists,
#  is a regular file, not size 0, user has read permission

if [ -e "$1" ] && [ -f "$1" ] && [ -s "$1" ] && [ -r "$1" ]
then
  # assign arg 1 to the varriable infile
  infile=$1
else
  echo "Input file DOES NOT exist OR is size zero OR IS NOT readable."
  exit 1
fi

echo "Infile is $infile"

if [ -e "$2" ] && [ -d "$2" ] && [ -s "$2" ] && [ -w "$2" ]
then
  # assign arg 1 to the varriable infile
  outdir=$2
else
  echo "Input file DOES NOT exist OR is not a directory OR IS NOT writeable."
  exit 1
fi

outfile="$outdir/mutsrefs.txt"

echo "Outdir is $outdir"
echo "Outfile is $outfile"

#Check that it is tab delimited with at least four columns (head contains three tabs)

if [ `head -n 1 $infile | awk -F \t {'print NF}'` -gt 0 ]
then
  echo "Header has at least four columns."
else
  echo "Header is not tab delimited or has too few columns"
  exit 1
fi

#########################################
### DO SEARCH & WRITE mRNA DATA TO FILE ##
###########################################
#Calculate some useful variables
firstline=`grep -v -n \# $infile | head -n 1 | awk -F : '{print $1}'`
nlines=`awk 'END {print NR}' $infile`
let nloops="($nlines - $firstline) / 4000 + 1"


echo "Number of lines in infile: $nlines"
echo "First data line: $firstline"
echo "Number of loops $nloops"

#USE HERE DOC FOR MYSQL QUERY
#Find number of rows that begin with #

#Get list of unique chrom values
chromos=`cut --fields=1 $infile | tail -n +$firstline | sort -u`


#create or overwrite output file with header from dummy querry
$sqlcmnd << HEREDOC | head -n 1 | cat > $outfile
SELECT knownGene.*, kgXref.geneSymbol, list.*  FROM knownGene
LEFT JOIN kgXref ON kgID = name
JOIN (Select 'demo' as chrom, 0 AS pos, 'A' as ref, 'T' as var ) AS list
LIMIT 1;
HEREDOC

echo "Doing Querys"


for ((i=1; i<=$nloops; i++))
do
  let nstart="$firstline + ($i-1) * 4000"
  echo "Doing Query $i"
  echo "nstart = $nstart"
  
#Execute query but add switch to supress column headers being included.
# $sqlcmnd --skip-column-names
$sqlcmnd --skip-column-names 1>> $outfile << MYHEREDOC
SELECT knownGene.*, kgXref.geneSymbol, list.* FROM knownGene
LEFT JOIN kgXref ON kgID = name
JOIN (SELECT 'demo' as 'chrom', 0 as 'pos', 'A' as 'ref', 'T' as 'var' 
`tail -n +$nstart $infile | head -n 4000 | cut --fields=1,2,3,4 | \
awk '{print "UNION ALL SELECT \x27" $1 "\x27," $2",\x27" $3 "\x27,\x27" $4 "\x27"}'`
) AS list ON list.chrom = knownGene.chrom AND list.pos < knownGene.txend AND list.pos > knownGene.txstart 
MYHEREDOC

done

echo "Wrote mutation and mRNA data to $outfile"

#Done
echo "Done."

exit