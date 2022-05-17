#!/bin/bash
#
#SBATCH -J refseq
#SBATCH -p devcore 
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2017-7-149
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=END
ulimit -c unlimited
set -eo pipefail


## Script used to randomly select ID of reference sequence (polarizing sequence) during the first iteration (N=1).
## Run independently for each replicate



# load modules
module load bioinfo-tools


#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# input folder
AA=${1:?msg}

# output folder 
RR=${2:?msg}

# Number of species
NSPEC=${3:?msg}

# subfolder name (it is also the root name for each library)
GO_FOLDER=${4:?msg}

## number of replicates
REP=${5:?msg}

## number of loci
NLOCUS=${6:?msg}

echo
echo $GO_FOLDER
echo $AA
echo "Number of species:" $NSPEC
echo "Number of replicates:" $REP
echo "Number of loci:" $NLOCUS
echo



for r_aux in $( eval echo {001..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

   ### create output file

   OUTF=${RR}/REPLICATE_${r_aux}_RefSequenceID.txt
   rm -f $OUTF

   ### Select random reference genome among listed species (excluding tetraploid and outgroup)
   cat ${AA}/REPLICATE_${r_aux}_SAMPLE_LIST.txt | head -n-1 | tail -n +2 > $SNIC_TMP/_REFspecies_aux
   
   aux_2="2"
   NSPEC_redux="$(awk '{print $1-$2}' <<<"$NSPEC $aux_2")"
   REF_i="$(echo $((1 + $RANDOM % $NSPEC_redux)))"
   
   REFspecies="$(cat $SNIC_TMP/_REFspecies_aux | head -n $REF_i | tail -n 1 )"

   echo $REFspecies > $OUTF

   echo "Replicate:" $r_aux
   echo "Selected reference species:" $REFspecies
   echo

done


# clean auxiliary file
rm -f $SNIC_TMP/_REFspecies_aux

