#!/bin/bash
#
#SBATCH -J preASTRAL
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited

 

## Script used to create ASTRAL input file.
## It collects the consensus trees obtained for each locus, using IQ-TEST2, that passed the filtering criteria (previous step).



# load modules
module load bioinfo-tools





######################################## paths and folder names

######### remember initial path
SRCDIR_INI=$(pwd)                                           	 

######### input folder 

declare -A AA
declare -A RR

AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P05a_phylogeny_IQ-TREE2_Trim_IUPAC_TETRA/BWA_unmasked_Alyrata	
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P05a_phylogeny_IQ-TREE2_Trim_IUPAC_TETRA/BWA_unmasked_Ahalleri-ensembl

######### output folder
RR[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P07a_phylogeny_ASTRAL_Trim_IUPAC_TETRA/BWA_unmasked_Alyrata_100
RR[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P07a_phylogeny_ASTRAL_Trim_IUPAC_TETRA/BWA_unmasked_Ahalleri-ensembl

######### Gene lists folder (loci that passed filtering criteria)
GLF=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P05_phylogeny_IQ-TREE2/INFORMATIVE_SITES_TETRA_100


######### Polyploid ID
POLYP="AKAMCHATICAKWS_DRR054581"	# allo
#POLYP="ASUECICAASO5_SRR2084157"	# allo


######### model name
MODEL_outfilename=MFP_ModelFinder


######### Ref genome (mapping)
ATYPE[1]=Alyrata_REF
ATYPE[2]=Ahalleri-ensembl_REF


######### subfolders containing input data (exons or introns) 
SUBFOLDER[1]=SINGLE_TRANSCRIPTS
SUBFOLDER[2]=SINGLE_INTRONS


######### FEATURE: EXON or INTRON
FT[1]=exons
FT[2]=introns





for n in `seq 1 1 2`; do		# genomic feature (exon intron)

   GO_FEATURE=${FT[$n]}

   for k in `seq 1 1 1`; do 			# cycle through ref genomes 

      GO_AA=${AA[$k]}/${SUBFOLDER[$n]}
      GO_RR=${RR[$k]}/${SUBFOLDER[$n]} 
      GO_ATYPE=${ATYPE[$k]}
      mkdir -p $GO_RR

      GO_POLYP=${POLYP}

      # NO SLIM
      AUX_prefix_out=""

      # Gene list
      GL=$GLF/00_Arabidopsis_gene_list_${GO_FEATURE}_${GO_ATYPE}_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}_TRIM_MAX-00000.txt
         

      ######## read gene list

      unset GENE_LIST
      x=1

      while read -r LINE                                                                      
      do
         GENE_LIST[$x]=$LINE
         let "x+=1"
      done < $GL

      N_GENES=${#GENE_LIST[@]}

   
      # input file suffix (IQ-TREE2 output file)
      InFileSUF=${GO_FEATURE}_trimN_SLIM-${GO_POLYP}-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.iqtree

      #output file
      OUTfile=GENE_TREES.${GO_FEATURE}_trimN_SLIM-${GO_POLYP}-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.newick


      rm -f ${GO_RR}/$OUTfile



      echo $InFileSUF
      echo $OUTfile
      echo $GL
      echo $GO_AA
      echo $GO_RR
      echo "Total number of genes:" $N_GENES
      echo

      cd $GO_AA

      for i in `seq 1 1 $N_GENES`; do                 # cycle through genes

         GO_GENE=${GENE_LIST[$i]}
         
         # get IQ-TEST2 consensus tree for each locus           
         cat ${GO_GENE}.${InFileSUF} | grep -A 2 'Consensus tree in newick format:' | tail -1 >> ${GO_RR}/${OUTfile}
      
      done
      
   done
done


echo "Done!"

