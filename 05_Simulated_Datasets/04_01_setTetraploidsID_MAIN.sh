#!/bin/bash
#

## Script used to randomly select ID of the two species that will be used to create the allotetraploip individual.
## Run independently for each replicate.


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


# Number of species (excluding root species)
NSPEC[1]=16
NSPEC[2]=16
NSPEC[3]=16
NSPEC[4]=16


## number of replicates
REP=100		


for k in `seq 1 1 4`; do 			# different simulation conditions

   GO_AA=${AA}/${INfolder[$k]}
   GO_RR=$SRCDIR_INI/ID_tetraploids/${INfolder[$k]}
   GO_FOLDER=${INfolder[$k]}
   GO_NSPEC=${NSPEC[$k]}
   
   echo $GO_NSPEC $REP
   echo $GO_FOLDER
   echo $GO_AA
   echo $GO_RR


   ./01_setTetraploidsID_query.sh $GO_AA \
                                        $GO_RR \
                                        $GO_NSPEC \
                                        $GO_FOLDER \
                                        $REP

   sleep 0.1
   echo
     
done





