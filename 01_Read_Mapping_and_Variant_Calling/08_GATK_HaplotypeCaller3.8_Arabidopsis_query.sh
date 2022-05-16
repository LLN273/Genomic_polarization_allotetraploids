#!/bin/bash
#
#SBATCH -J GATK_HC
#SBATCH -p core 
#SBATCH -n 10
#SBATCH -t 48:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



## Script used to perform variant calling using GATK; first step: HaplotypeCaller


module load bioinfo-tools
module load GATK/3.8-0
module load samtools/1.6



############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)   

#read BAM file (with path)
READ1=${1:?msg}

#read BAM file (without path)
BAM1=${2:?msg}		

#output file (root name)
RR=${3:?msg}

#sample name
SNAME=${4:?msg}

# reference genome
refGenome=${5:?msg}			

# reference genome folder
refGenomefolder=${6:?msg}

# ploidy
SPLOIDY=${7:?msg}

# Intervals list (BED format) exome data 
ANNF=${8:?msg}



echo
echo "Input samples (BAM):" $READ1
echo "Sample name:" $SNAME
echo "Reference genome:" $refGenome
echo "Reference genome folder:" $refGenomefolder
echo "Output folder:" $RR
echo "Sample ploidy:" $SPLOIDY
echo "Intervals list:" $ANNF
echo






########### copy reference genome to scratch disk

cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%.fa*}* .






################# Create intervals list (target genes)

### Convert position 0 to position 1 (GATK does not accept position 0)
echo
awk 'BEGIN {OFS="\t"} $2 == "0" {$2 = "1"}; 1'  $ANNF > $SNIC_TMP/__aux0.txt

### Convert intervals list to chr1:100-200 format
echo
cat $SNIC_TMP/__aux0.txt | cut -f1 > $SNIC_TMP/__aux1.txt
cat $SNIC_TMP/__aux0.txt | cut -f2 > $SNIC_TMP/__aux2.txt
cat $SNIC_TMP/__aux0.txt | cut -f3 > $SNIC_TMP/__aux3.txt
paste -d':' $SNIC_TMP/__aux1.txt $SNIC_TMP/__aux2.txt > $SNIC_TMP/__auxA.txt
paste -d'-' $SNIC_TMP/__auxA.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__targets.interval_list






##### Run GATK HaplotypeCaller in GVCF mode

	
java -Xmx50G -Djava.io.tmpdir=$SNIC_TMP -jar $GATK_HOME/GenomeAnalysisTK.jar \
	    -T HaplotypeCaller \
	    --emitRefConfidence GVCF \
	    -ploidy $SPLOIDY \
	    --intervals $SNIC_TMP/__targets.interval_list \
	    --num_cpu_threads_per_data_thread $SLURM_NTASKS \
	    -R $SNIC_TMP/${refGenome} \
	    -I $READ1 \
	    -o $RR/${SNAME}-individualHaplotypeCaller.g.vcf

	




      
# Notes 

# For HELP, type:
# java -jar $GATK_HOME/GenomeAnalysisTK.jar --help
# See also: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

# -T,--analysis_type <analysis_type>                	Name of the tool to run
#						 	HaplotypeCaller: Call SNPs and indels via local re-assembly of haplotypes
# --emitRefConfidence					Mode for emitting reference confidence scores
#							Records whether the trimming intervals are going to be used to emit reference confidence, 
#							{@code true}, or regular HC output {@code false}. 
#							Mode for emitting reference confidence scores:
#							NONE: 			Regular calling without emitting reference confidence calls.
#							BP_RESOLUTION:	Reference model emitted site by site.
#							GVCF:		    Reference model emitted with condensed non-variant blocks, i.e. the GVCF format. 
# -genotyping_mode -gt_mode 	 			Specifies how to determine the alternate alleles to use for genotyping  [DISCOVERY]
# --sample_ploidy -ploidy 				Ploidy per sample [2]	
# --num_cpu_threads_per_data_thread	-nct		Number of CPU threads to allocate per data thread

# --min_mapping_quality_score -mmq 	20 		Minimum read mapping quality required to consider a read for analysis with the HaplotypeCaller [20]
# -rf UnmappedRead					Filters out unmapped reads
# -rf DuplicateRead					Filters out duplicate reads
# -rf MappingQualityUnavailable				Filters out reads with no mapping quality information [APPLIED BY DEFAULT]
# -rf NotPrimaryAlignment				Filters out reads tagged as secondary alignments 
# --min_base_quality_score   -mbq 	 		Minimum base quality required to consider a base for calling [10]
# -standard_min_confidence_threshold_for_calling / -stand_call_conf	[10]			The minimum phred-scaled confidence threshold at which variants should be called. Only variant sites with QUAL equal or greater than this threshold will be called. Note that
#																					since version 3.7, we no longer differentiate high confidence from low confidence calls at the calling step. The default call confidence threshold is set low intentionally to
# 																					achieve high sensitivity, which will allow false positive calls as a side effect. Be sure to perform some kind of filtering after calling to reduce the amount of false positives
#																					in your final callset. Note that when HaplotypeCaller is used in GVCF mode (using either -ERC GVCF or -ERC BP_RESOLUTION) the call threshold is automatically set to zero. Call
#																					confidence thresholding will then be performed in the subsequent GenotypeGVCFs command.
#
# --unsafe / -U	ALLOW_SEQ_DICT_INCOMPATIBILITY			Enable unsafe operations: nothing will be checked at runtime.  'ALLOW_SEQ_DICT_INCOMPATIBILITY' allows for contig legth differences between dict index file and lengths listed in SAM file

# -I,--input_file <input_file>                       	Input file containing sequence data (BAM or CRAM)







#### clean scratch disk
rm -f $SNIC_TMP/${SNAME}*
rm -f $SNIC_TMP/${refGenome%.fa*}*


