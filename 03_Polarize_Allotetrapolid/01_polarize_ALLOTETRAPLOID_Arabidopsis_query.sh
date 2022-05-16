#!/bin/bash
#
#SBATCH -J polarize
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=END
ulimit -c unlimited



## Script used to polarize allotetraploid sequence in the MSA.
## Polarization done separately for each locus.


# load modules
module load bioinfo-tools
module load python3/3.8.7



#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# input folder
AA=${1:?msg}

# output folder 
RR=${2:?msg}

# Number of species
NSPEC=${3:?msg}

# Allotetraploid ID
POLYP=${4:?msg}

# Ref sequence ID (polarizing sequence)
REFspecies=${5:?msg}

## gene list
GL=${6:?msg}

## genomic feature
GO_FEATURE=${7:?msg}

# exclude species
GO_EXCLUDE=${8:?msg}



echo
echo $GO_FEATURE
echo $AA
echo $RR
echo $GL
echo "Number of species:" $NSPEC
echo "Polyploid:" $POLYP
echo "Reference sequence (polarizing sequence):" $REFspecies
echo "Exclude:" $GO_EXCLUDE
echo





######## read gene list

i=1

while read -r LINE                                                                      
do
   GENE_LIST[$i]=$LINE
   let "i+=1"
done < ${GL}

N_GENES=${#GENE_LIST[@]}





cd $RR

for i in `seq 1 1 $N_GENES`; do

   # prepare MSA
   GO_MSA=${GENE_LIST[$i]}.${GO_FEATURE}_trimN.fa
   AUX_fa=${GO_MSA%.fa}_SLIM-${POLYP}.fa
   sed "/^>${GO_EXCLUDE}/{N;d}" $AA/$GO_MSA > $RR/$AUX_fa

   # adjust number of species
   NSPEC_aux=`expr $NSPEC - 1`



   echo $i $GO_MSA

   # polarize sequences (first check whether source MSA is empty)
   if [ -s $RR/$AUX_fa ]
   then
      python3 01s_polarizeTETRA.py $AUX_fa $NSPEC_aux $REFspecies $POLYP
   else
      touch $RR/${AUX_fa%.fa}-ALT.fasta    	#if MSA file is empty
   fi



done


