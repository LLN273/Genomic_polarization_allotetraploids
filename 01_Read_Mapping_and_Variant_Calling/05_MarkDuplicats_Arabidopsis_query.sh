#!/bin/bash
#
#SBATCH -J duplicates
#SBATCH -p core 
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



## Script used to identify read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artifactual processes. 
## This step eliminates PCR duplicates (arising from library preparation) across all lanes in addition to optical duplicates (which are by definition only per-lane).


# Load software
module load bioinfo-tools
module load picard/2.10.3
module load samtools/1.6







########### input data


#read input folder name
AA=${1:?msg}                            

# output folder
RR=${2:?msg}			

#sample name
SNAME=${3:?msg}

echo
echo "Sample:" $SNAME
echo

echo "Output folder:" $RR
echo


# current path
SRCDIR_INI=$(pwd)   











########## Merge BAM files associated to the same sample

if [[ $SNAME == "Aarenosaarenosa1_SRR208278X" ]]; then
	
      java -Xmx5G -jar $PICARD_HOME/picard.jar MergeSamFiles \
	                                        INPUT=$AA/Aarenosaarenosa1_SRR2082782.LANE-8_MergeBamAlignment.bam \
	                                        INPUT=$AA/Aarenosaarenosa1_SRR2082785.LANE-7_MergeBamAlignment.bam \
						SORT_ORDER=coordinate \
	                                        O=$SNIC_TMP/${SNAME}.bam

elif [[ $SNAME == "AkamchaticaKWS_DRR054581" ]]; then

      java -Xmx5G -jar $PICARD_HOME/picard.jar MergeSamFiles \
	                                        INPUT=$AA/AkamchaticaKWS_DRR054581.LANE-2_MergeBamAlignment.bam \
	                                        INPUT=$AA/AkamchaticaKWS_DRR054581.LANE-7_MergeBamAlignment.bam \
						SORT_ORDER=coordinate \
	                                        O=$SNIC_TMP/${SNAME}.bam

elif [[ $SNAME == "Alyratapetraea3_SRR2040797" ]]; then

      java -Xmx5G -jar $PICARD_HOME/picard.jar MergeSamFiles \
	                                        INPUT=$AA/Alyratapetraea3_SRR2040797.LANE-2_MergeBamAlignment.bam \
	                                        INPUT=$AA/Alyratapetraea3_SRR2040797.LANE-3_MergeBamAlignment.bam \
						INPUT=$AA/Alyratapetraea3_SRR2040797.LANE-4_MergeBamAlignment.bam \
						SORT_ORDER=coordinate \
	                                        O=$SNIC_TMP/${SNAME}.bam

elif [[ $SNAME == "Alyratapetraea15_SRR2040828" ]]; then

      java -Xmx5G -jar $PICARD_HOME/picard.jar MergeSamFiles \
	                                        INPUT=$AA/Alyratapetraea15_SRR2040828_RUN1.LANE-2_MergeBamAlignment.bam \
	                                        INPUT=$AA/Alyratapetraea15_SRR2040828_RUN1.LANE-3_MergeBamAlignment.bam \
						INPUT=$AA/Alyratapetraea15_SRR2040828_RUN1.LANE-4_MergeBamAlignment.bam \
						INPUT=$AA/Alyratapetraea15_SRR2040828_RUN2.LANE-8_MergeBamAlignment.bam \
						SORT_ORDER=coordinate \
	                                        O=$SNIC_TMP/${SNAME}.bam 

else
	
   echo
   ls $AA/${SNAME}.LANE-*_MergeBamAlignment.bam

   rsync -ah $AA/${SNAME}.LANE-*_MergeBamAlignment.bam $SNIC_TMP/${SNAME}.bam
   rsync -ah $AA/${SNAME}.LANE-*_MergeBamAlignment.bam $SNIC_TMP/${SNAME}.bai

fi





########## Mark duplicates
echo
cd $SNIC_TMP
mkdir __RD_aux1

java -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
                                         INPUT=$SNIC_TMP/${SNAME}.bam \
                                         OUTPUT=$RR/${SNAME}_marked_duplicates.bam \
                                         M=$SNIC_TMP/${SNAME}_marked_dup_metrics.txt \
                                         TMP_DIR=$SNIC_TMP/__RD_aux1 
   
   
   


## Notes picard.jar MarkDuplicates
## INPUT (String)							One or more input SAM or BAM files to analyze. Must be coordinate sorted. Default value: null. This option may be specified 0 or more times.
## OUTPUT (File)							The output file to write marked records to Required.
## METRICS_FILE (File)						File to write duplication metrics to Required.
## REMOVE_DUPLICATES (Boolean)				If true do not write duplicates to the output file instead of writing them with appropriate flags set. Default value: false. This option can be set to 'null' to 
##											clear the default value. Possible values: {true, false}
## REMOVE_SEQUENCING_DUPLICATES (Boolean)	If true remove 'optical' duplicates and other duplicates that appear to have arisen from the sequencing process instead of the library preparation process, 
##											even if REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored. Default value: false. This option can 
##											be set to 'null' to clear the default value. Possible values: {true, false}




   


########## index BAM file
echo
java -Xmx5G -jar $PICARD_HOME/picard.jar BuildBamIndex \
                                         INPUT=$RR/${SNAME}_marked_duplicates.bam
		
		
		
		
		
								 
	
   
########### glean statistics about BAM file
echo
samtools flagstat $RR/${SNAME}_marked_duplicates.bam  > $RR/${SNAME}_marked_duplicates.flagstat


										 

## clean scratch disk
rm -rf $SNIC_TMP/__RD_aux1
rm -f $SNIC_TMP/*.bam
rm -f $SNIC_TMP/*.bai
rm -f $SNIC_TMP/*.txt








