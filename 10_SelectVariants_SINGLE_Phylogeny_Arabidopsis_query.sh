#!/bin/bash
#
#SBATCH -J SelectVar
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


# Script used to split VCF file into two separate files, one for SNP+invariant sites and the other containing only INDEL variants


module load bioinfo-tools
module load GATK/3.8-0







############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)


#read VCF file (with path)
READ1=${1:?msg}

#read VCF file (without path)
VCF1=${2:?msg}		

#output file (root name)
RR=${3:?msg}

# interval list (BED format) [target genes]
ANNF=${4:?msg}

# reference genome
refGenome=${5:?msg}			

# reference genome folder
refGenomefolder=${6:?msg}




### Output name (root)
OUTF=${VCF1%.vcf}



echo
echo "VCF file:" $READ1
echo "Interval list:" $ANNF
echo "Reference genome:" $refGenome
echo "Reference genome folder:" $refGenomefolder
echo "Output folder:" $RR
echo
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








################ Select variants	



####### Extract SNPs and NON-VARIANT sites from call set
java -Xmx5G -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $SNIC_TMP/${refGenome} \
    -V $READ1 \
    -selectType SNP \
    -selectType NO_VARIATION \
    --intervals $SNIC_TMP/__targets.interval_list \
    -o $RR/${OUTF}_snps.vcf



echo
echo

####### Extract INDELS from call set
java -Xmx5G -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $SNIC_TMP/${refGenome} \
    -V $READ1 \
    -selectType INDEL \
    --intervals $SNIC_TMP/__targets.interval_list \
    -o $RR/${OUTF}_indels.vcf



echo
echo



#### NOTE ON VARIANT TYPES
# INDEL
# SNP
# MIXED		Reference = 'T', Sample = 'A,TCC'
# MNP (aka complex substitution) 	Reference = 'GC', Sample = 'TTA'
# SYMBOLIC	* or <NONREF>
# NO_VARIATION





# For HELP, type:
# java -jar $GATK_HOME/GenomeAnalysisTK.jar --help
# See also: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php
#
# -T SelectVariants							elect a subset of variants from a larger callset
#--selectTypeToInclude (-selectType)		Select only a certain type of variants from the input file
# --selectexpressions  (-select)			One or more criteria to use when selecting the data
# --num_threads -nt							Number of data threads to allocate to this analysis [1]
# --reference_sequence -R					Reference sequence file
# --variant	-V								Input VCF file






#### clean scratch disk
rm -f $SNIC_TMP/${OUTF}*
rm -f $SNIC_TMP/${refGenome%.fa*}*



