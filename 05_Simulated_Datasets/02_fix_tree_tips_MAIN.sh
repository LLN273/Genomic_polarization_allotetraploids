#!/bin/bash


## Script used to fix tree tips so that the two (haploid) individuals belonging to the same species are monophyletic 
## (as these will be recoded as one single individual at a later stage)


####################################### paths and folder names

##### remember initial path
SRCDIR_INI=$(pwd)                                           	 


##### input folder (root)

AA=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/01_simphy_ILS


##### Input/output subfolder (different sumulation conditions)
INfolder[1]=A_modILS_recSpeciation_16t
INfolder[2]=B_modILS_deepSpeciation_16t
INfolder[3]=C_HIGH-ILS_recSpeciation_16t
INfolder[4]=D_HIGH-ILS_deepSpeciation_16t


# Number of species (including outgroup and tetraploid; excludes allotetraploid's parents)
NSPEC[1]=16
NSPEC[2]=16
NSPEC[3]=16
NSPEC[4]=16


## number of replicates
REP=100		


# number of loci/genes
NLOCUS=1200


echo

for k in `seq 1 1 4`; do 			# different simulation conditions

      echo
      GO_NSPEC=${NSPEC[$k]}
      GO_FOLDER=${INfolder[$k]}
      GO_AA=${AA}/${GO_FOLDER}
      GO_RR=${AA}/${GO_FOLDER}_CLEAN
      mkdir -p $GO_RR

      cd $GO_RR
      rsync -ah $GO_AA/*.command .
      rsync -ah $GO_AA/*.db .
      rsync -ah $GO_AA/*.params .


      for r_aux in $( eval echo {001..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

         mkdir -p $GO_RR/$r_aux
         cd $GO_RR/$r_aux
         rsync -ah $SRCDIR_INI/02_fix_tree_tips.R .

         echo $r_aux $GO_NSPEC $NLOCUS
         echo $GO_AA/$r_aux
         echo $GO_RR/$r_aux

         cd $SRCDIR_INI
         sbatch $SRCDIR_INI/02_fix_tree_tips_query.sh $GO_AA/$r_aux \
                                                      $GO_RR/$r_aux \
                                                      $r_aux \
                                                      $GO_NSPEC \
                                                      $NLOCUS
                                          

         sleep 0.1
         echo

      done
done


