#!/bin/bash


## Script used to enforce absence of short internal branches (quasi-polytomies).
## This was done by equalizing branch lengths.



####################################### paths and folder names

##### remember initial path
SRCDIR_INI=$(pwd)                                           	 


##### input folder 
AA=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/01_simphy_ILS

##### Input/output subfolder
INfolder[1]=X_NoILS_16t


##### Output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/20_No-ILS-No-polytomy_Simulated_Dataset


##### Number of species
NSPEC[1]=16


##### number of replicates
REP=100		


##### number of loci/genes
NLOCUS=1



echo

for k in `seq 1 1 1`; do 			# different simulation conditions

      echo
      GO_NSPEC=${NSPEC[$k]}
      GO_FOLDER=${INfolder[$k]}
      GO_AA=${AA}/${GO_FOLDER}
      GO_RR=${RR}/${GO_FOLDER}_CLEAN
      mkdir -p $GO_RR

      cd $GO_RR
      rsync -ah $GO_AA/*.command .
      rsync -ah $GO_AA/*.db .
      rsync -ah $GO_AA/*.params .


      for r_aux in $( eval echo {001..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_002, etc)

         mkdir -p $GO_RR/$r_aux
         cd $GO_RR/$r_aux
         rsync -ah $SRCDIR_INI/02_equalize_branch_length.R .

         echo $r_aux $GO_NSPEC $NLOCUS
         echo $GO_AA/$r_aux


         cd $SRCDIR_INI
         $SRCDIR_INI/02_equalize_branch_length_query.sh $GO_AA/$r_aux \
                                                        $GO_RR/$r_aux \
                                                        $r_aux \
                                                        $GO_NSPEC \
                                                        $NLOCUS
                                          

         sleep 0.1
         echo

      done
done









