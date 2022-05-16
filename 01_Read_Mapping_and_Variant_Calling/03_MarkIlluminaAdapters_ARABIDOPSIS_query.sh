#!/bin/bash -l
#SBATCH -J MarkIlluminaAdapters
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



## Script used to identify and mark position of Illumina adapters



# load modules
module load bioinfo-tools
module load picard/2.10.3





#read input file
BAM1=${1:?msg}                            

# output folder
RR=${2:?msg}			

#sample name
SNAME=${3:?msg}

echo
echo "Sample:"
echo $SNAME
echo

# current path
SRCDIR=$(pwd)                                                                    











########## Mark adapter sequences
####
#### MORE INFO:
#### https://gatkforums.broadinstitute.org/gatk/discussion/6483
#### https://gatkforums.broadinstitute.org/gatk/discussion/6484#latest#top
## Extra notes:
## picard's MarkIlluminaAdapters produces two files. 
## (1) The metrics file, XXXXX_metrics.txt bins the number of tagged adapter bases versus the number of reads. 
## (2) The XXXXX_markilluminaadapters.bam file is identical to the input BAM, XXXXX_revertsam.bam, except reads with adapter sequences will be marked with a tag in XT:i:# format, where # denotes the 5' starting position of the adapter sequence. 
##     At least six bases are required to mark a sequence. Reads without adapter sequence remain untagged.




cd $SNIC_TMP
mkdir __MIA_aux



########## Mark adapter sequences using MarkIlluminaAdapters
java -Xmx5G -jar $PICARD_HOME/picard.jar MarkIlluminaAdapters \
                                            INPUT=$BAM1 \
                                            OUTPUT=$RR/${SNAME}_markilluminaadapters.bam \
                                            METRICS=$RR/${SNAME}_markilluminaadapters_metrics.txt \
                                            TMP_DIR=$SNIC_TMP/__MIA_aux
   
   
   
   
# notes: 
# INPUT (File)				Required.
# OUTPUT (File)				If output is not specified, just the metrics are generated Default value: null.
# METRICS (File)			Histogram showing counts of bases_clipped in how many reads Required.
# TMP_DIR					    Working folder used to process large files






## clean scratch disk
rm -rf $SNIC_TMP/__MIA_aux





