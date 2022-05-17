#!/bin/bash
#
#SBATCH -J bcfConsensus
#SBATCH -p core 
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL


## Script used to create a consensus sequence for each sample


module load bioinfo-tools
module load bcftools/1.12
module load BEDOPS/2.4.39
module load BEDTools/2.29.2



#input folder
AA=${1:?msg}

#read file name (without path)
InFile1=${2:?msg}		

#output folder (root name)
RR=${3:?msg}

#sample name
SNAME=${4:?msg}

# reference genome
refGenome=${5:?msg}			

# reference genome folder
refGenomefolder=${6:?msg}

# interval list (BED format) [target genes]  >> Not in use
ANNF=${7:?msg}


#remember initial path
SRCDIR_INI=$(pwd)                                           	


echo
echo "Sample:" $SNAME
echo "VCF file (SNPs):" $AA/$InFile1
echo "VCF file (indels):" $AA/${SNAME}_indels.vcf.gz
echo "Reference genome:" $refGenomefolder/$refGenome $SNAME		
echo "Interval list:" $ANNF
echo "output folder:" $RR
echo






##### Move reference genome to scratch disk 
cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%.fa*}* .

gunzip -c $AA/$InFile1 > $SNIC_TMP/${SNAME}.vcf








#### 1. Variants
#### Note: all (./.) sites have already been set to FAIL
echo
echo '1. Preparing variants file' ${SNAME}_variants.vcf
cat $SNIC_TMP/${SNAME}.vcf | grep '^#' > $SNIC_TMP/${SNAME}_variants.vcf
cat $SNIC_TMP/${SNAME}.vcf | grep -v '^#' | grep PASS > $SNIC_TMP/${SNAME}_variants_aux.vcf
awk 'length($4) == 1 { print $0 }' $SNIC_TMP/${SNAME}_variants_aux.vcf > $SNIC_TMP/${SNAME}_variants_aux2.vcf	# reference allele must be a single nucleotide 
awk '$4 != "N" { print $0 }' $SNIC_TMP/${SNAME}_variants_aux2.vcf > $SNIC_TMP/${SNAME}_variants_aux3.vcf	# excludes masked sites (N) 
awk '$5 !~ /*/ { print }' $SNIC_TMP/${SNAME}_variants_aux3.vcf >> $SNIC_TMP/${SNAME}_variants.vcf		# excludes alternate sites that include '*' ('*' but also 'A,*' sites etc -- but this is a tiny number) >> necessary otherwise these sites are REMOVED from consensus sequence !!







##### 2. Adjust vcf file format and index (required by bcftools)
echo
echo '2. Indexing variants file'
bgzip -c $SNIC_TMP/${SNAME}_variants.vcf > $SNIC_TMP/${SNAME}_variants.vcf.gz
bcftools index $SNIC_TMP/${SNAME}_variants.vcf.gz
bcftools view $SNIC_TMP/${SNAME}_variants.vcf.gz -Oz -o $SNIC_TMP/${SNAME}_variants_ADJ.vcf.gz
bcftools index $SNIC_TMP/${SNAME}_variants_ADJ.vcf.gz



##### 3. Create consensus sequence 

echo
echo '3. Running bcftools consensus'
bcftools consensus $SNIC_TMP/${SNAME}_variants_ADJ.vcf.gz \
                    --fasta-ref $SNIC_TMP/$refGenome \
                    --iupac-codes \
                    --sample $SNAME \
                    -o $SNIC_TMP/${SNAME}-bcftoolsConsensus.fa




## Notes:
#  -f, --fasta-ref FILE			reference sequence in fasta format 
#  -e, --exclude EXPRESSION		exclude sites for which EXPRESSION is true.     
#  -M, --missing CHAR			instead of skipping the missing genotypes, output the character CHAR (e.g. "?") 
#  -s, --sample NAME			apply variants of the given sample 
#  -o, --output FILE                    write output to a file 
#  -I, --iupac-codes			output variants in the form of IUPAC ambiguity codes 
#  -H, --haplotype 1|2|R|A|I|LR|LA|SR|SA|1pIu|2pIu	choose which allele from the FORMAT/GT field to use (the codes are case-insensitive):
#		1	the first allele, regardless of phasing 
#		2	the second allele, regardless of phasing 
#		R	the REF allele (in heterozygous genotypes) 
#		A	the ALT allele (in heterozygous genotypes) 
#		I	IUPAC code for all genotypes 
#		LR, LA	the longer allele. If both have the same length, use the REF allele (LR), or the ALT allele (LA) 
#		SR, SA	the shorter allele. If both have the same length, use the REF allele (SR), or the ALT allele (SA) 
#		1pIu, 2pIu	first/second allele for phased genotypes and IUPAC code for unphased genotypes. This option requires *-s*, unless exactly one sample is present in the VCF








##### 4a. Create snp-mask file (SNP)
##### This file contains all ./. sites and failed variants
##### These sites will be set to N in the consensus sequence 
echo
echo '4a. Creating snp-mask file' ${SNAME}_snpmask.vcf
cat $SNIC_TMP/${SNAME}.vcf | grep '^#' > $SNIC_TMP/${SNAME}_snpmask.vcf
cat $SNIC_TMP/${SNAME}.vcf | grep -v '^#' | grep -v PASS > $SNIC_TMP/${SNAME}_snap-mask_aux.txt
awk 'length($4) == 1 { print $0 }' $SNIC_TMP/${SNAME}_snap-mask_aux.txt > $SNIC_TMP/${SNAME}_snap-mask_aux2.txt	# reference allele must be a single nucleotide 
awk '$4 != "N" { print $0 }' $SNIC_TMP/${SNAME}_snap-mask_aux2.txt > $SNIC_TMP/${SNAME}_snap-mask_aux3.txt	# excludes masked sites (N) 
awk '$5 !~ /*/ { print }' $SNIC_TMP/${SNAME}_snap-mask_aux3.txt >> $SNIC_TMP/${SNAME}_snpmask.vcf		# excludes alternate sites that include '*' ('*' but also 'A,*' sites etc -- but this is a tiny number) >> necessary otherwise these sites are REMOVED from consensus sequence !!



##### 4b. Create snp-mask file (deletions only)
echo
echo '4b. Create snp-mask file (deletions only)'
VCF_indels=${SNAME}_indels.vcf.gz
zcat $AA/$VCF_indels | grep '^#' > $SNIC_TMP/VCF_indels1.vcf
zcat $AA/$VCF_indels | grep -v '^#' | grep PASS | grep -vE "([.]/[.])" > $SNIC_TMP/__aux_pass1.txt
awk '$5 !~ /*/ { print }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt		    # excludes alternate sites that include '*' ('*' but also 'A,*' sites etc -- but this is a tiny number) >> necessary otherwise these sites are REMOVED from consensus sequence !!
awk '(length($4) > 1 && length($4) > length($5)) { print $0 }' $SNIC_TMP/__aux_pass2.txt >> $SNIC_TMP/VCF_indels1.vcf    # include deletions only
vcf2bed --deletions < $SNIC_TMP/VCF_indels1.vcf | cut -f1-3 > $SNIC_TMP/VCF_indels1.bed		# convert to bed format






#### 5. Use bedtools to mask ./. sites and failed variants

echo
echo '5a. Running bedtools maskfasta'
bedtools maskfasta \
	     -fi $SNIC_TMP/${SNAME}-bcftoolsConsensus.fa \
	     -bed $SNIC_TMP/${SNAME}_snpmask.vcf \
	     -fo $SNIC_TMP/${SNAME}-bcftoolsConsensus_aux.fa

echo
echo '5b. Running bedtools maskfasta, deletions'
bedtools maskfasta \
	     -fi $SNIC_TMP/${SNAME}-bcftoolsConsensus_aux.fa \
	     -bed $SNIC_TMP/VCF_indels1.bed \
	     -fo $SNIC_TMP/${SNAME}-bcftoolsConsensus_FINAL.fa


#Notes:
#  -fi 		<input FASTA>
# -bed 		<BED/GFF/VCF>
# -fo 		<output FASTA>
# -soft 	Soft-mask (that is, convert to lower-case bases) the FASTA sequence. By default, hard-masking (that is, conversion to Ns) is performed.
# -mc 		Replace masking character. That is, instead of masking with Ns, use another character.





#### 6. Stats
echo 'Compiling stats'

# Number of nucleotides in reference genome
aux_nucR="$(cat $SNIC_TMP/$refGenome | grep -v '^>' | grep -io [A-Z] | wc -l )"
echo  "Number of nucleotides in reference genome:" $aux_nucR

# Number of nucleotides in consensus sequence (before masking)
aux_nucC="$(cat $SNIC_TMP/${SNAME}-bcftoolsConsensus.fa | grep -v '^>' | grep -io [A-Z] | wc -l )"
echo  "Number of nucleotides in consensus sequence (before masking failed variants):" $aux_nucC

# Number of nucleotides in consensus sequence (FINAL)
aux_nucF="$(cat $SNIC_TMP/${SNAME}-bcftoolsConsensus_FINAL.fa | grep -v '^>' | grep -io [A-Z] | wc -l )"
echo  "Number of nucleotides in consensus sequence (FINAL):" $aux_nucF

# Number of PASS sites (variant or invariant)
aux_var="$(cat $SNIC_TMP/${SNAME}_variants.vcf | grep -v '^#' | wc -l )"
echo  "Number of PASS sites (variant or invariant):" $aux_var

# Total number of ./. sites and failed variants (excludes sites marked as N in reference genome)
aux_N2="$(cat $SNIC_TMP/${SNAME}_snpmask.vcf | grep -v '^#' | cut -f4 | grep -v N | wc -l )"
echo "Total number of ./. sites and failed variants (excludes sites marked as N in reference genome):" $aux_N2

# Number of masked sites in reference genome
aux_nucRN="$(cat $SNIC_TMP/$refGenome | grep -v '^>' | grep -io [N] | wc -l )"
echo  "Number of masked sites in reference genome:" $aux_nucRN

# Number of masked sites in consensus sequence (FINAL)
aux_nucFN="$(cat $SNIC_TMP/${SNAME}-bcftoolsConsensus_FINAL.fa | grep -v '^>' | grep -io [N] | wc -l )"
echo  "Number of masked sites in consensus sequence (FINAL):" $aux_nucFN

# Difference in number of masked nucleotides between consensus and reference
N_DIFF="$(awk '{print $1-$2}' <<<"$aux_nucFN $aux_nucRN")"
echo  "Difference in number of masked nucleotides between consensus and reference:" $N_DIFF




##### Copy files to results folder
echo 'Saving files to output folder'
cd $RR
rsync -ah $SNIC_TMP/${SNAME}_variants.vcf .
rsync -ah $SNIC_TMP/${SNAME}_snpmask.vcf .
#rsync -ah $SNIC_TMP/${SNAME}_snpmask.bed .
gzip -c $SNIC_TMP/${SNAME}-bcftoolsConsensus.fa > $RR/${SNAME}-bcftoolsConsensus.fa.gz
gzip -c $SNIC_TMP/${SNAME}-bcftoolsConsensus_FINAL.fa > $RR/${SNAME}-bcftoolsConsensus_FINAL.fa.gz
rsync -ah  $SNIC_TMP/VCF_indels1.bed $RR/${SNAME}_indels.bed


##### Remove auxiliary files
cd $SNIC_TMP
rm -f $SNIC_TMP/${SNAME}*
rm -f $SNIC_TMP/${refGenome%.fa*}*
rm -f $SNIC_TMP/_*
rm -f $SNIC_TMP/VCF*


echo 'Done!'


