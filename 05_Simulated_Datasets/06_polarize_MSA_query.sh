#!/bin/bash
#
#SBATCH -J polarize
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2017-7-149
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail


## Script used to process simulated MSAs:
## (1) creates diploid species (collapses sequences for the two haploid individuals belonging to the same species);
## (2) creates allotetraploid sequence (collapses sequences for two diploid species);
## (3) polarizes allotetraploid sequence.

### Input MSA files must contain two haploid sequences per species (only one sequence for outgroup [0_0_0])
### It is assumed that input file contains NO indels (i.e. sequences are aligned)
###


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

# subfolder name (it is also the root name for each library)
GO_FOLDER=${4:?msg}

## replicate number
r_aux=${5:?msg}

## number of loci
NLOCUS=${6:?msg}

## folder with tetraploid and reference sequence info
BB=${7:?msg}

## Iteration 
NALT=${8:?msg}


echo
echo $GO_FOLDER
echo $AA
echo $BB
echo $RR
echo "Replicate number:" $r_aux
echo "Number of species:" $NSPEC
echo "Number of loci:" $NLOCUS
echo "Iteration:" $NALT



# read reference sequence ID
REFspecies="$(cat ${BB}/REPLICATE_${r_aux}_RefSequenceID.txt )"

echo 'Reference sequence ID:' $REFspecies

# read tetraploid ID (parental species ID)
anc_1="$(cat ${BB}/REPLICATE_${r_aux}_SAMPLE_LIST.txt | tail -1 | cut -f2 -d'-' )"
anc_2="$(cat ${BB}/REPLICATE_${r_aux}_SAMPLE_LIST.txt | tail -1 | cut -f3 -d'-' )"

echo 'Parental species IDs:' $anc_1 $anc_2
echo

mkdir -p $RR/$r_aux
cd $RR/$r_aux

rsync -ah $SRCDIR_INI/06s_polarizeMSA_TETRA.py .


for l_aux in $( eval echo {001..${NLOCUS}} ); do		# run for each locus; must include leading zeros (LOCUS_0001, LOCUS_0002, etc)

   echo $l_aux
   
   # clean input MSA
   cat $AA/$r_aux/data_${l_aux}_TRUE.fasta | tr -d "[:blank:]" | grep -v '^$' > $RR/$r_aux/data_${l_aux}_TRUE_CLEAN.fasta

   # collapse sequences and polarize allotetraploid
   python3 06s_polarizeMSA_TETRA.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies $anc_1 $anc_2

done








# Runtime: < 10 min

