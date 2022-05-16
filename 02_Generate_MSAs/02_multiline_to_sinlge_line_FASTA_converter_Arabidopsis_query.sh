#!/bin/bash
#
#SBATCH -J slFASTA
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL


## Script used to convert consensus sequences (fasta format) from multiline to single line.



module load bioinfo-tools



#read file name (with path)
READ1=${1:?msg}

#read file name (without path)
InFile1=${2:?msg}		

#output folder (root name)
RR=${3:?msg}

#sample name
sampleName=${4:?msg}




#remember initial path
SRCDIR_INI=$(pwd)                                           	


echo
echo $sampleName
echo $READ1
echo


# uncompress fasta
gunzip -c $READ1 > $SNIC_TMP/${sampleName}.fa


# convert to single line fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $SNIC_TMP/${sampleName}.fa > $SNIC_TMP/${sampleName}-bcftoolsConsensus_FINAL_singleLine_AUX.fa


##note: must edit output file and remove first line (empty)
cat $SNIC_TMP/${sampleName}-bcftoolsConsensus_FINAL_singleLine_AUX.fa | tail -n +2 > $SNIC_TMP/${sampleName}-bcftoolsConsensus_FINAL_singleLine.fa


# compress fasta
gzip -c $SNIC_TMP/${sampleName}-bcftoolsConsensus_FINAL_singleLine.fa > $RR/${sampleName}-bcftoolsConsensus_FINAL_singleLine.fa.gz


## clean auxiliary files
cd $SNIC_TMP
rm -f $SNIC_TMP/${sampleName}*


