#!/bin/bash
#
#SBATCH -J map_bwa
#SBATCH -p core 
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


## Script used to convert bam to fastq, then pipe into BWA. Finally, merge original uBAM and BAM produced by BWA to generate a clean BAM



# Load software
module load bioinfo-tools
module load picard/2.10.3
module load bwa/0.7.17
module load samtools/1.6






########### input data


#read input file
INBAM=${1:?msg}                            

# output folder
RR=${2:?msg}			

#sample name
SNAME=${3:?msg}

# reference genome
refGenome=${4:?msg}			

# reference genome folder
refGenomefolder=${5:?msg}

echo
echo "Sample:" $SNAME
echo "Input BAM:" $INBAM
echo "Reference genome:" $refGenome
echo "Reference genome folder:" $refGenomefolder
echo "Output folder:" $RR
echo

# current path
SRCDIR_INI=$(pwd)   






########### copy reference genome to scratch disk

cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%.fa*}* .




#### Convert bam to fastq, then pipe into BWA. Finally, merge original uBAM and BAM produced by BWA to generate a clean BAM

   
cd $SNIC_TMP
mkdir __BWA_aux1
mkdir __BWA_aux3
   
   
java -Xmx40G -jar $PICARD_HOME/picard.jar SamToFastq \
                                          INPUT=$INBAM \
                                          FASTQ=/dev/stdout \
                                          CLIPPING_ATTRIBUTE=XT \
   					  CLIPPING_ACTION=2 \
					  INTERLEAVE=true \
					  NON_PF=true \
                                          TMP_DIR=$SNIC_TMP/__BWA_aux1 | \
bwa mem -t $SLURM_NTASKS \
            -M \
   	    -p \
   	    $SNIC_TMP/$refGenome \
            /dev/stdin | \
java -Xmx40G -jar $PICARD_HOME/picard.jar  MergeBamAlignment \
                                              ALIGNED_BAM=/dev/stdin \
                                              UNMAPPED_BAM=$INBAM \
                                              OUTPUT=$RR/${SNAME}_MergeBamAlignment.bam \
                                              REFERENCE_SEQUENCE=$SNIC_TMP/$refGenome \
    					      CREATE_INDEX=true \
     					      ADD_MATE_CIGAR=true \
                                              CLIP_ADAPTERS=false \
 					      CLIP_OVERLAPPING_READS=true \
                                              INCLUDE_SECONDARY_ALIGNMENTS=true \
 					      MAX_INSERTIONS_OR_DELETIONS=-1 \
                                              PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
 					      ATTRIBUTES_TO_RETAIN=XS \
                                              TMP_DIR=$SNIC_TMP/__BWA_aux3







########### glean statistics about BAM file
samtools flagstat $RR/${SNAME}_MergeBamAlignment.bam > $RR/${SNAME}_MergeBamAlignment.flagstat

   # total number of reads reported in flagstat file include:
   #	2* paired reads aligned concordantly exactly 1 time
   #	2* paired aligned concordantly >1 times
   #	2* paired reads aligned discordantly 1 time
   #	and, from pairs aligned 0 times concordantly or discordantly:
   #		1* reads aligned exactly 1 time
   #		1* aligned >1 times






########## clean temporary folder
rm -rf $SNIC_TMP/__BWA_aux1
rm -rf $SNIC_TMP/__BWA_aux3

   









##### Notes: picard.jar SamToFastq
# INPUT (File)					Input SAM/BAM file to extract reads from Required.
# FASTQ (File)					Output FASTQ file (single-end fastq or, if paired, first end of the pair FASTQ). Required. Cannot be used in conjuction with option(s) OUTPUT_PER_RG (OPRG)
# CLIPPING_ATTRIBUTE (String)	The attribute that stores the position at which the SAM record should be clipped Default value: null.
#								NOTE: In the current run, the clipping attributed, XT, has been set when we marked the position of Illumina adapters
# CLIPPING_ACTION (String)		The action that should be taken with clipped reads: 
#									'X' means the reads and qualities should be trimmed at the clipped position; 
#									'N' means the bases should be changed to Ns in the clipped region; 
#									and any integer means that the base qualities should be set to that value in the clipped region. Default value: null.  [NOTE: this is the option we are using here]
# INTERLEAVE (Boolean)			Will generate an interleaved fastq if paired, each line will have /1 or /2 to describe which end it came from Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
#INCLUDE_NON_PRIMARY_ALIGNMENTS (Boolean)	If true, include non-primary alignments in the output. Support of non-primary alignments in SamToFastq is not comprehensive, so there may be exceptions if this is set to true and there are paired reads 
#											with non-primary alignments. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}



##### Notes: BWA
# -t INT        number of threads [1]
# -R STR        Complete read group header line. It is particularly important to ensure that the @RG information here is correct as this information is used by later tools. 
#				The SM field must be set to the name of the sample being processed, and LB field to the library. 
#				'\t' can be used in STR and will be converted to a TAB in the output SAM. 
#				The read group ID will be attached to every read in the output.  
#				Example: '@RG\tID:foo\tSM:bar' [null]
# -M			Mark shorter split hits as secondary (for Picard compatibility).
#   			The BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments for different part of a query sequence. This is a crucial feature
#				for long sequences. However, some tools such as Picard's markDuplicates does not work with split alignments. One may consider to use option -M to flag shorter
#				split hits as secondary.
# -p        	smart pairing (ignoring in2.fq)
# -a			Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.



##### Notes: MergeBamAlignment
# ALIGNED_BAM (File)						SAM or BAM file(s) with alignment data. Default value: null. This option may be specified 0 or more times. Cannot be used in conjuction with option(s) READ1_ALIGNED_BAM (R1_ALIGNED) READ2_ALIGNED_BAM (R2_ALIGNED)
# UNMAPPED_BAM (File)						Original SAM or BAM file of unmapped reads, which must be in queryname order. Required.
# OUTPUT (File)								Merged SAM or BAM file to write to. Required.
# REFERENCE_SEQUENCE (File)					Path to the fasta file for the reference sequence. Required.
# ADD_MATE_CIGAR (Boolean)					Adds the mate CIGAR tag (MC) if true, does not if false. Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# CLIP_ADAPTERS (Boolean)					Whether to clip adapters where identified. Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# CLIP_OVERLAPPING_READS (Boolean)			For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate. Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# INCLUDE_SECONDARY_ALIGNMENTS (Boolean)	If false, do not write secondary alignments to output. Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# MAX_INSERTIONS_OR_DELETIONS (Integer)		The maximum number of insertions or deletions permitted for an alignment to be included. Alignments with more than this many insertions or deletions will be ignored. Set to -1 to allow any number of insertions or deletions. Default value: 1. This option can 
#											be set to 'null' to clear the default value.
# PRIMARY_ALIGNMENT_STRATEGY (PrimaryAlignmentStrategy)	Strategy for selecting primary alignment when the aligner has provided more than one alignment for a pair or fragment, and none are marked as primary, more than one is marked as primary, or the primary alignment is filtered out for some reason.
#												- BestMapq expects that multiple alignments will be correlated with HI tag, and prefers the pair of alignments with the largest MAPQ, in the absence of a primary selected by the aligner. 
#												- EarliestFragment prefers the alignment which maps the earliest base in the read. Note that EarliestFragment may not be used for paired reads. 
#												- BestEndMapq is appropriate for cases in which the aligner is not pair-aware, and does not output the HI tag. It simply picks the alignment for each end with the highest MAPQ, and makes those alignments primary, regardless of whether the two alignments 
#												make sense together.
#												- MostDistant is also for a non-pair-aware aligner, and picks the alignment pair with the largest insert size. If all alignments would be chimeric, it picks the alignments for each end with the best MAPQ. 
#												For all algorithms, ties are resolved arbitrarily. Default value: BestMapq. This option can be set to 'null' to clear the default value. Possible values: {BestMapq, EarliestFragment, BestEndMapq, MostDistant}
# ATTRIBUTES_TO_RETAIN (String)				Reserved alignment attributes (tags starting with X, Y, or Z) that should be brought over from the alignment data when merging. Default value: null. This option may be specified 0 or more times.








## clean scratch disk
rm -f $SNIC_TMP/*.bam
rm -f $SNIC_TMP/*.sam
rm -f $SNIC_TMP/*.fa*
rm -f $SNIC_TMP/${refGenome%.fa*}*


