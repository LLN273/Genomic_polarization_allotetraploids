#!/bin/bash
#


## Script used to randomly select ID of reference sequence (polarizing sequence) during the first iteration (N=1).
## Run independently for each replicate.



#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# Input folder (root)
AA=$SRCDIR_INI/ID_tetraploids

# Input/output subfolder (different simulation conditions)
INfolder[1]=A_modILS_recSpeciation_16t_CLEAN
INfolder[2]=B_modILS_deepSpeciation_16t_CLEAN
INfolder[3]=C_HIGH-ILS_recSpeciation_16t_CLEAN
INfolder[4]=D_HIGH-ILS_deepSpeciation_16t_CLEAN


# Number of species
NSPEC[1]=16
NSPEC[2]=16
NSPEC[3]=16
NSPEC[4]=16


## number of replicates
REP=100		

# number of loci/genes
NLOCUS=1200



for k in `seq 1 1 4`; do 			# different simulation conditions

   GO_FOLDER=${INfolder[$k]}
   GO_AA=${AA}/${GO_FOLDER}
   GO_RR=$GO_AA
   GO_NSPEC=${NSPEC[$k]}
   
   echo $GO_NSPEC $REP $NLOCUS
   echo $GO_FOLDER
   echo $GO_AA
   echo $GO_RR

   ./05_setRefSequenceID_query.sh $GO_AA \
                                        $GO_RR \
                                        $GO_NSPEC \
                                        $GO_FOLDER \
                                        $REP \
                                        $NLOCUS

   sleep 0.1
   echo
     
done

