#!/bin/bash
#
ulimit -c unlimited



## Script used to enforce absence of short internal branches (quasi-polytomies).
## This was done by equalizing branch lengths.
## Also, fix tree tips so that the two (haploid) individuals belonging to the same species are monophyletic  
## (as these will be recoded as one single diploid individual at a later stage)



# load modules
module load bioinfo-tools
module load R/4.0.0




#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

#input folder
GO_AA=${1:?msg}	

#output folder (root name)
GO_RR=${2:?msg}

#Replicate
r_aux=${3:?msg}

# Number of species
GO_NSPEC=${4:?msg}

# number of loci/genes	
NLOCUS=${5:?msg}



echo
echo 'Replicate:' $r_aux
echo 'Number of species:' $GO_NSPEC
echo 'Number of loci/genes:' $NLOCUS
echo 'Input folder:' $GO_AA
echo 'Output folder:' $GO_RR
echo



cd $GO_RR

rsync -ah $GO_AA/l_trees.trees .
rsync -ah $GO_AA/s_tree.trees .


for i in $( eval echo {0001..${NLOCUS}} ); do		# locus/genes; must include leading zeros (LOCUS_0001, LOCUS_02, etc)

   ## input phylogeny 
   gene_tree="$(cat $GO_AA/g_trees${i}.trees)"

   ## output file
   OUTFile="g_trees${i}.trees"

   echo "Replicate, gene, No. species:" ${r_aux} $i $GO_NSPEC

   ## remove tip, if needed
   Rscript --no-save ./02_equalize_branch_length.R $gene_tree \
                                                   $GO_RR \
                                                   $OUTFile \
                                                   $GO_NSPEC


   ## if tip removed above, add new tip close to remaining species tip
   if test -f "$GO_RR/id_kept.txt"; then
      unset TIP_LIST
      x=1

      while read -r LINE                                                                      
      do
         TIP_LIST[$x]=$LINE
         let "x+=1"
      done < $GO_RR/id_kept.txt

      N_TIPS=${#TIP_LIST[@]}

      for n in `seq 1 1 $N_TIPS`; do
         tip_ref=${TIP_LIST[$n]}
         tip_sis=${tip_ref%_0}_1
         tip_group="(${tip_ref}:0.15000000,${tip_sis}:0.15000000)"
         if grep --quiet ",${tip_ref}" "$OUTFile"; then 
            sed -e "s/,${tip_ref}/,${tip_group}/" $OUTFile > ./_aux_outfile1.txt
            mv ./_aux_outfile1.txt $OUTFile
         elif grep --quiet "(${tip_ref}" "$OUTFile"; then
            sed "s/(${tip_ref}/(${tip_group}/" $OUTFile > ./_aux_outfile1.txt
            mv ./_aux_outfile1.txt $OUTFile
         fi
      done
   fi

   #scan tree branch lengths, convert from scientific notation to numerical (otherwise INDELIble crashes during the next step)
   # approach: add tabs (\t) before and after any of "(", ")", "," or ":". This separates all fields by tabs. Then do conversion to numerical (for each separate field, first checking how they are formated). Finally, delete tabs.
   sed -e 's/(/\t(\t/g' $OUTFile > ./_aux_outfile_num1.txt
   sed -e 's/)/\t)\t/g' ./_aux_outfile_num1.txt > ./_aux_outfile_num2.txt
   sed -e 's/,/\t,\t/g' ./_aux_outfile_num2.txt > ./_aux_outfile_num3.txt
   sed -e 's/:/\t:\t/g' ./_aux_outfile_num3.txt > ./_aux_outfile_num4.txt
   awk 'BEGIN{CONVFMT="%.9f"; FS=OFS="\t"}{for(i=1; i<=NF; i++)if($i~/[1-9.]+[eE][+-][0-9]+?/)$i+=0;}1' ./_aux_outfile_num4.txt > ./_aux_outfile_num5.txt
   sed -e 's/\t//g' ./_aux_outfile_num5.txt > ./_aux_outfile_num6.txt
   mv ./_aux_outfile_num6.txt $OUTFile


   rm -f $GO_RR/id_kept.txt
   rm -f $GO_RR/_*

done


