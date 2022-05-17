#!/bin/bash
#


## Script used to process simulated MSAs:
## (1) creates diploid species (collapses sequences for the two haploid individuals belonging to the same species);
## (2) creates allotetraploid sequence (collapses sequences for two diploid species);
## (3) polarizes allotetraploid sequence.

### Input MSA files must contain two haploid sequences per species (only one sequence for outgroup [0_0_0])
### It is assumed that input file contains NO indels (i.e. sequences are aligned)
###



#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# Input folder (root)
AA=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/01_simphy_ILS

# Input/output subfolder (different simulation conditions)
INfolder[1]=A_modILS_recSpeciation_16t_CLEAN
INfolder[2]=B_modILS_deepSpeciation_16t_CLEAN
INfolder[3]=C_HIGH-ILS_recSpeciation_16t_CLEAN
INfolder[4]=D_HIGH-ILS_deepSpeciation_16t_CLEAN


# Folder with tetraploid and reference sequence info (root name)
BB=$SRCDIR_INI/ID_tetraploids

# output folder (root)
RR=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/03_polarized_sequences

# iteration
NALT=1

# Number of species (including root species)
NSPEC[1]=17
NSPEC[2]=17
NSPEC[3]=17
NSPEC[4]=17

## number of replicates
REP=100		

# number of loci
NLOCUS=1200



for k in `seq 1 1 4`; do 			# different simulation conditions

   GO_FOLDER=${INfolder[$k]}
   GO_AA=${AA}/${GO_FOLDER}
   GO_RR=${RR}/${GO_FOLDER}
   GO_NSPEC=${NSPEC[$k]}
   
   for r_aux in $( eval echo {001..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

      echo $r_aux $GO_NSPEC $REP $NLOCUS
      echo $GO_FOLDER
      echo $GO_AA
      echo $GO_RR


      sbatch ./06_polarize_MSA_query.sh $GO_AA \
                                        $GO_RR \
                                        $GO_NSPEC \
                                        $GO_FOLDER \
                                        $r_aux \
                                        $NLOCUS \
                                        $BB/$GO_FOLDER \
                                        $NALT

      sleep 0.1
      echo

   done     
done





