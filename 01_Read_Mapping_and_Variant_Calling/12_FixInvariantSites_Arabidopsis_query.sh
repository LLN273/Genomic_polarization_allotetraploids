#!/bin/bash
#
#SBATCH -J FixInvariants
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 3:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail


### Script used to fix invariant sites in VCF file
### i) if invariant site has DP>0 in GT filed, then PASS and 0/0
### ii) if invariant site has DP=0 in GT field, then FAIL and ./.




module load bioinfo-tools
module load python3/3.6.0





############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)


# input folder
AA=${1:?msg}

#read VCF file (without path)
VCF1=${2:?msg}		

# output folder
RR=${3:?msg}




echo
echo "VCF file:" $VCF1
echo "Input folder:" $AA
echo "Output folder:" $RR
echo
echo






# Fix invariant sites
cd $SNIC_TMP
gunzip -c $AA/$VCF1 > $SNIC_TMP/${VCF1%.gz}
rsync -ah $SRCDIR_INI/12s_fix_invariant_sites.py .

python3 12s_fix_invariant_sites.py ${VCF1%.gz} 




# copy results to output folder
cd $AA
gzip -c $SNIC_TMP/${VCF1%.vcf.gz}-CLEAN.vcf > $RR/${VCF1%.vcf.gz}-CLEAN.vcf.gz


# clean scratch disk
rm -f $SNIC_TMP/12s_fix_invariant_sites.py
rm -f $SNIC_TMP/${VCF1%.gz}
rm -f $SNIC_TMP/*-CLEAN.vcf


