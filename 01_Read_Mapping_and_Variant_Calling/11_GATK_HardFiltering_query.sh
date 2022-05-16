#!/bin/bash
#
#SBATCH -J HardFilter
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail


module load bioinfo-tools
module load GATK/3.8-0



# Script used to perform hard filtering on SNPs, invariant sites, and indel variants



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

# mapper flag
FLAG_MAPPER=${7:?msg}

# Variant type
VART=${8:?msg}

# Ploidy
PLOIDY=${9:?msg}

### Output name (root)
OUTF=${VCF1%.vcf}


echo
echo "VCF file:" $READ1
echo "Variant type:" $VART
echo "Ploidy:" $PLOIDY
echo "Mapper:" $FLAG_MAPPER
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






####### Filer variants



### SNPs

if [[ $FLAG_MAPPER = "BWA" && $VART = "snps" && $PLOIDY = "2" ]] ; then

    echo

    java -Xmx5G -jar $GATK_HOME/GenomeAnalysisTK.jar \
                -T VariantFiltration \
                -R $SNIC_TMP/${refGenome} \
                -V $READ1 \
				-filterName "QD_filter" \
				-filter "QD < 2.0" \
				-filterName "FS_filter" \
				-filter "FS > 10.0" \
				-filterName "MQ_filter" \
				-filter "MQ < 40.0" \
				-filterName "SOR_filter" \
				-filter "SOR > 5.0" \
				-filterName "MQRankSum_filter_min" \
				-filter "MQRankSum < -2.5" \
				-filterName "MQRankSum_filter_max" \
				-filter "MQRankSum > 2.5" \
				-filterName "ReadPosRankSum_filter_min" \
				-filter "ReadPosRankSum < -2.5" \
				-filterName "ReadPosRankSum_filter_max" \
				-filter "ReadPosRankSum > 2.5" \
				-filterName "Depth_filter_min" \
				-filter "DP < 3.0" \
	            --intervals $SNIC_TMP/__targets.interval_list \
				--downsampling_type NONE \
				--logging_level ERROR \
                -o $RR/${OUTF}-filtered_tranche_2.vcf



fi 


if [[ $FLAG_MAPPER = "BWA" && $VART = "snps" && $PLOIDY = "4" ]] ; then

    echo

    java -Xmx5G -jar $GATK_HOME/GenomeAnalysisTK.jar \
                -T VariantFiltration \
                -R $SNIC_TMP/${refGenome} \
                -V $READ1 \
				-filterName "QD_filter" \
				-filter "QD < 2.0" \
				-filterName "FS_filter" \
				-filter "FS > 10.0" \
				-filterName "MQ_filter" \
				-filter "MQ < 40.0" \
					-filterName "SOR_filter" \
					-filter "SOR > 5.0" \
					-filterName "MQRankSum_filter_min" \
					-filter "MQRankSum < -2.5" \
					-filterName "MQRankSum_filter_max" \
					-filter "MQRankSum > 2.5" \
					-filterName "ReadPosRankSum_filter_min" \
					-filter "ReadPosRankSum < -2.5" \
					-filterName "ReadPosRankSum_filter_max" \
					-filter "ReadPosRankSum > 2.5" \
					-filterName "Depth_filter_min" \
					-filter "DP < 5.0" \
	            --intervals $SNIC_TMP/__targets.interval_list \
				--downsampling_type NONE \
				--logging_level ERROR \
                -o $RR/${OUTF}-filtered_tranche_2.vcf
	

fi 






### indels

if [[ $VART = "indels" && $PLOIDY = "2" ]] ; then

    echo

    java -Xmx5G -jar $GATK_HOME/GenomeAnalysisTK.jar \
                -T VariantFiltration \
                -R $SNIC_TMP/${refGenome} \
                -V $READ1 \
				-filterName "QD_filter" \
				-filter "QD < 2.0" \
				-filterName "FS_filter" \
				-filter "FS > 10.0" \
					-filterName "SOR_filter" \
					-filter "SOR > 5.0" \
					-filterName "MQRankSum_filter_min" \
					-filter "MQRankSum < -2.5" \
					-filterName "MQRankSum_filter_max" \
					-filter "MQRankSum > 2.5" \
					-filterName "ReadPosRankSum_filter_min" \
					-filter "ReadPosRankSum < -2.5" \
					-filterName "ReadPosRankSum_filter_max" \
					-filter "ReadPosRankSum > 2.5" \
					-filterName "Depth_filter_min" \
					-filter "DP < 3.0" \
	            --intervals $SNIC_TMP/__targets.interval_list \
				--downsampling_type NONE \
				--logging_level ERROR \
                -o $RR/${OUTF}-filtered_tranche_2.vcf
	

fi 


if [[ $VART = "indels" && $PLOIDY = "4" ]] ; then

    echo

    java -Xmx5G -jar $GATK_HOME/GenomeAnalysisTK.jar \
                -T VariantFiltration \
                -R $SNIC_TMP/${refGenome} \
                -V $READ1 \
				-filterName "QD_filter" \
				-filter "QD < 2.0" \
				-filterName "FS_filter" \
				-filter "FS > 10.0" \
					-filterName "SOR_filter" \
					-filter "SOR > 5.0" \
					-filterName "MQRankSum_filter_min" \
					-filter "MQRankSum < -2.5" \
					-filterName "MQRankSum_filter_max" \
					-filter "MQRankSum > 2.5" \
					-filterName "ReadPosRankSum_filter_min" \
					-filter "ReadPosRankSum < -2.5" \
					-filterName "ReadPosRankSum_filter_max" \
					-filter "ReadPosRankSum > 2.5" \
					-filterName "Depth_filter_min" \
					-filter "DP < 5.0" \
	            --intervals $SNIC_TMP/__targets.interval_list \
				--downsampling_type NONE \
				--logging_level ERROR \
                -o $RR/${OUTF}-filtered_tranche_2.vcf
	

fi 







#### clean scratch disk
rm -f $SNIC_TMP/${OUTF}*
rm -f $SNIC_TMP/${refGenome%.fa*}*
rm -f $SNIC_TMP/__*





