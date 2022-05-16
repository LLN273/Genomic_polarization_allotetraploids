#!/bin/bash -l
#SBATCH -J fqTOuBAM
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


# Script converts fastq files to unmapped BAM (uBAM) format.
# Adds read-group tags.



# load modules
module load bioinfo-tools
module load picard/2.10.3
#module load samtools/1.6




#read input files
READ1=${1:?msg}                            
READ2=${2:?msg}   

# output folder
RR=${3:?msg}			

#library name
LIBNAME=${4:?msg}

#header type
HEADERT=${5:?msg}




####### sample name

# Aarenosaarenosa1_SRR2082782.LANE-8 and Aarenosaarenosa1_SRR2082785.LANE-7 are associated to the same sample >> Aarenosaarenosa1_SRR208278X
# Same goes for Alyratapetraea15_SRR2040828_RUN1 and Alyratapetraea15_SRR2040828_RUN2, which are both associated to the same sample
if [[ "$LIBNAME" = "Aarenosaarenosa1_SRR2082782.LANE-8" ]] ; then
   SNAME="Aarenosaarenosa1_SRR208278X"
elif [[ "$LIBNAME" = "Aarenosaarenosa1_SRR2082785.LANE-7" ]] ; then
   SNAME="Aarenosaarenosa1_SRR208278X"
elif [[ "$LIBNAME" = "Alyratapetraea15_SRR2040828_RUN1.LANE-2" ]] ; then
   SNAME="Alyratapetraea15_SRR2040828"
elif [[ "$LIBNAME" = "Alyratapetraea15_SRR2040828_RUN1.LANE-3" ]] ; then
   SNAME="Alyratapetraea15_SRR2040828"
elif [[ "$LIBNAME" = "Alyratapetraea15_SRR2040828_RUN1.LANE-4" ]] ; then
   SNAME="Alyratapetraea15_SRR2040828"
elif [[ "$LIBNAME" = "Alyratapetraea15_SRR2040828_RUN2.LANE-8" ]] ; then
   SNAME="Alyratapetraea15_SRR2040828"
else
   SNAME=$(echo $LIBNAME | cut -f1 -d '.')
fi



echo
echo "Library:" $LIBNAME
echo "Sample:" $SNAME
echo "Header type:" $HEADERT
echo "Input files:"
echo $READ1
echo $READ2



# current path
SRCDIR=$(pwd)                                                                    




#### Create uBAM for each lane
####
#### MORE INFO:
#### https://gatkforums.broadinstitute.org/gatk/discussion/6483
#### https://gatkforums.broadinstitute.org/gatk/discussion/6484#latest#top
#### https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
#### https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam







########## library number
########## Sample Aarenosaarenosa1 has two libraries (SRR2082782 & SRR2082785); all other samples have just one library

if [[ "$LIBNAME" = "Aarenosaarenosa1_SRR2082785.LANE-7" ]] ; then
   LIB_No="2"
else
   LIB_No="1"
fi







########## glean read group information   

if [[ "$HEADERT" = "0" ]] ; then

   # There is NO info on header
   # type 0: @SRR1945833.1 SRR492245.1/1
	   
   RG_ID=UNKNOWN.${SNAME}.1   			# read group identifier (aka read group name)
   SM=$SNAME				                # sample name
   LB=${SNAME}_LIB${LIB_No}			    # Library identifier >> use a different lib number (eg LIB2) if two DNA libraries are prepared for the same sample 
   PU=$RG_ID					              # platform unit
   PL=ILLUMINA					            # platform used to generate the data
   SC=NA					                  # sequencing center



elif [ "$HEADERT" = "1" ]; then
   #### lane is in field 4 (separator ':') 
   #### type 1 >>> @SRR2040797.1 BS-DSFCONTROL03:172:D17EJACXX:2:1101:1419:2084/1
	  
   # example of line format: @NS500198:23:H0VUMAGXX:3:11401:11643:1016 1:N:0:5
	  
   flowcell_ID=$(zcat $READ1 | head -1 | cut -d ' ' -f2 | cut -d ':' -f3)	    # flowcell id 			(H0VUMAGXX)
   flowcell_lane=$(zcat $READ1 | head -1 | cut -d ' ' -f2 | cut -d ':' -f4)	  # flowcell lane			(3)
   RG_ID=${flowcell_ID}.${SNAME}.${flowcell_lane}				# Because some barcodes have troublesome characters (eg '<' or ':'), we will use the sample name as the barcode
   SM=$SNAME									                          # sample name
   LB=${SNAME}_LIB${LIB_No}							                # Library identifier (during DNA preparation) >> use a different lib number (eg LIB2) if two DNA libraries are prepared for the same sample 
   PU=${flowcell_ID}.${SNAME}.${flowcell_lane}					# Because some barcodes have troublesome characters (eg '<' or ':'), we will use the sample name as the barcode
   PL=ILLUMINA									                        # platform used to generate the data
   SC=NA									                              # sequencing center 


elif [ "$HEADERT" = "2" ]; then
   #### lane is in field 2 (separator ':')
   #### type 2 >>> @SRR2040777.1 HWI-ST1176_0190:6:1101:4596:1991/1

   # example of line format: @ERR179410.1 FCC0BFYACXX:2:1101:1247:2097/1
	   
   RG_ID=$(zcat $READ1 | head -1 | cut -d ' ' -f2 | cut -d ':' -f1-2 | tr ':' '.')  		# read group identifier (aka read group name)
   SM=$SNAME											# sample name
   LB=${SNAME}_LIB${LIB_No}				# Library identifier >> use a different lib number (eg LIB2) if two DNA libraries are prepared for the same sample 
   PU=$RG_ID										  # platform unit
   PL=ILLUMINA										# platform used to generate the data
   SC=NA									        # sequencing center
   																									
fi
   
echo
echo "flowcell_lane:" $flowcell_lane
echo "read group identifier:" $RG_ID
echo "Library identifier:" $LB
echo "platform unit:" $PU
echo 




########## Convert FASTQ to unmapped BAM (uBAM) and add read group information

java -Xmx10G -jar $PICARD_HOME/picard.jar FastqToSam \
                                             FASTQ=$READ1 \
                                             FASTQ2=$READ2 \
                                             OUTPUT=$RR/${LIBNAME}_fastqtosam.bam \
                                             READ_GROUP_NAME=${RG_ID} \
                                             SAMPLE_NAME=${SM} \
                                             LIBRARY_NAME=${LB} \
                                             PLATFORM_UNIT=${PU} \
                                             PLATFORM=${PL} \
                                             SEQUENCING_CENTER=${SC} \
					                                   SORT_ORDER=queryname
   
   





   
# notes (picard.jar FastqToSam): 
# FASTQ (File)					Input fastq file (optionally gzipped) for single end data, or first read in paired end data. Required.
# FASTQ2 (File)					Input fastq file (optionally gzipped) for the second read of paired end data. Default value: null.
# OUTPUT (File)					Output SAM/BAM file. Required.
# READ_GROUP_NAME (String)		Read group name Default value: A. This option can be set to 'null' to clear the default value.
# SAMPLE_NAME (String)			Sample name to insert into the read group header Required.
# LIBRARY_NAME (String)			The library name to place into the LB attribute in the read group header Default value: null.
# PLATFORM_UNIT (String)		The platform unit (often run_barcode.lane) to insert into the read group header Default value: null.
# PLATFORM (String)				The platform type (e.g. illumina, solid) to insert into the read group header Default value: null.
# SEQUENCING_CENTER (String)	The sequencing center from which the data originated Default value: null.
# SORT_ORDER (SortOrder)		The sort order for the output sam/bam file. Default value: queryname. 
#								This option can be set to 'null' to clear the default value. Possible values: {unsorted, queryname, coordinate, duplicate, unknown}









## clean scratch disk
rm -f $SNIC_TMP/*.bam


