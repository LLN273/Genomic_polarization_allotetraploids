#!/bin/bash
#
#SBATCH -J tetra
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2020-15-40
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=END
ulimit -c unlimited
set -eo pipefail


## Script used to randomly select ID of the two species that will be used to create the allotetraploip individual.
## Run independently for each replicate.


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



echo
echo $GO_FOLDER
echo $AA
echo "Number of species:" $NSPEC
echo "Number of replicates:" $REP
echo


mkdir -p $RR

for r_aux in $( eval echo {001..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

   ### Select ancestral species for allotetraploid (can be any species pair with the exception of the root species [0])
   ANC_01="$(echo $((1 + $RANDOM % $NSPEC)))"
   ANC_02=$ANC_01
   while [ "${ANC_02}" -eq "${ANC_01}" ]
   do
     ANC_02="$(echo $((1 + $RANDOM % $NSPEC)))"
   done

   echo "Replicate:" $r_aux
   echo "Parent 1:" $ANC_01
   echo "Parent 2:" $ANC_02
   echo


   ### create new sample list (ancestral species removed; allotetraploid ID added to sample list)

   NEW_SL=$RR/REPLICATE_${r_aux}_SAMPLE_LIST.txt
   rm -f $SNIC_TMP/_NEW_SL

   for k in `seq 0 1 $NSPEC`; do
      echo $k >> $SNIC_TMP/_NEW_SL
   done

   cat $SNIC_TMP/_NEW_SL | grep -vx ${ANC_01} | grep -vx ${ANC_02} > $NEW_SL
   echo T-${ANC_01}-${ANC_02} >> $NEW_SL


done


# clean auxiliary files
rm -f $SNIC_TMP/_NEW_SL

