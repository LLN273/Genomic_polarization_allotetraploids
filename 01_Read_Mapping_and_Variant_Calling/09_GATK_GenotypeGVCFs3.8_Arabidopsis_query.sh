#!/bin/bash
#
#SBATCH -J GATK_SG
#SBATCH -p core 
#SBATCH -n 3
#SBATCH -t 04:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



## Script used to perform variant calling using GATK; second step: GenotypeGVCFs


module load bioinfo-tools
module load GATK/3.8-0
#module load samtools/1.6


############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)   

# gVCF file (wih path)
IN_FILE=${1:?msg}

# gVCF input file (without path)
IN_FILE_name=${2:?msg}
	
#output file
RR=${3:?msg}

# Sample name
SNAME=${4:?msg}

# reference genome
refGenome=${5:?msg}

# reference genome folder
refGenomefolder=${6:?msg}










echo
echo "Sample name:" $SNAME
echo "Input file:" $IN_FILE
echo "Output folder:" $RR
echo "Reference genome:" $refGenome
echo "Reference genome folder:" $refGenomefolder
echo




########### copy reference genome to scratch disk

cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%.fa*}* .





##### Run HaplotypeCaller variant calling in GVCF mode INDIVIDUALLY FOR EACH SAMPLE (not in joint genotyping mode)
##### IMPORTANT NOTICE: must includeNonVariantSites when preparing data for phylogenetic analysis

java -Xmx40G -Djava.io.tmpdir=$SNIC_TMP -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    --num_threads 3 \
    --includeNonVariantSites \
    --variant $IN_FILE \
    -R $SNIC_TMP/$refGenome \
    -o $RR/${SNAME}.vcf











#     --sample_ploidy 2 \
# Special note on ploidy (from https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php):
# This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy for non-diploid organisms.
#
# From slurm file:
# GenotypeGVCFs - Notice that the -ploidy parameter is ignored in GenotypeGVCFs tool as this is automatically determined by the input variant files


      
# Notes GATK

# For HELP, type:
# java -jar $GATK_HOME/GenomeAnalysisTK.jar --help
# See also: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php

# GATK: Read filters
# These Read Filters are automatically applied to the data by the Engine before processing by GenotypeGVCFs.
#	MalformedReadFilter
#	BadCigarFilter
#	UnmappedReadFilter
#	NotPrimaryAlignmentFilter
#	FailsVendorQualityCheckFilter
#	DuplicateReadFilter





# -T,--analysis_type <analysis_type>                	Name of the tool to run
#							GenotypeGVCFs: Perform joint genotyping on gVCF files produced by HaplotypeCaller 

#							GenotypeGVCFs merges gVCF records that were produced as part of the Best Practices workflow for variant 
#							discovery (see Best Practices documentation for more details) using the '-ERC GVCF' or '-ERC 
#							BP_RESOLUTION' mode of the HaplotypeCaller, or result from combining such gVCF files using 
#							CombineGVCFs. This tool performs the multi-sample joint aggregation step and merges the records 
#							together in a sophisticated manner: at each position of the input gVCFs, this tool will combine all 
#							spanning records, produce correct genotype likelihoods, re-genotype the newly merged record, and then 
#							re-annotate it.
# --num_threads -nt					Number of data threads to allocate to this analysis [1]
# --includeNonVariantSites -allSites			Include loci found to be non-variant after genotyping
# --variant -V						One or more input gVCF files
# --reference_sequence -R				Reference sequence file

# -rf UnmappedRead					Filter out unmapped reads
# -rf DuplicateRead					Filter out duplicate reads
# -rf NotPrimaryAlignment				Filter out read records that are secondary alignments. This filter recognizes the SAM flag that 
#							identifies secondary alignments. It is intended to ensure that only records that are likely to be 
#							mapped in the right place, and therefore to be informative, will be used in analysis. To be clear, it 
#							does NOT filter out read records that are supplementary alignments.








#### clean scratch disk
rm -r $SNIC_TMP/*.fa*





