#!/bin/bash
#
#SBATCH -J ASTRAL
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail
 

## Script used to estimate species tree using ASTRAL, based on gene trees inferred using IQ-TREE2 



# load modules
module load bioinfo-tools
module load newick_utils/20160413

ASTRAL=/home/luisleal/MYAPPS/astral.5.7.3/Astral/astral.5.7.3.jar



###################################### paths and folder names

###### remember initial path
SRCDIR_INI=$(pwd)                                           	 



###### input folder 

declare -A AA
declare -A RR

AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P07a_phylogeny_ASTRAL_Trim_IUPAC_TETRA/BWA_unmasked_Alyrata_100
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P07a_phylogeny_ASTRAL_Trim_IUPAC_TETRA/BWA_unmasked_Ahalleri-ensembl


###### Allopolyploid ID
POLYP="AKAMCHATICAKWS_DRR054581"	# allo
#POLYP="ASUECICAASO5_SRR2084157"	# allo


###### model name
MODEL_outfilename=MFP_ModelFinder


###### subfolders containing input data (exons or introns) 
SUBFOLDER[1]=SINGLE_TRANSCRIPTS
SUBFOLDER[2]=SINGLE_INTRONS


###### FEATURE: EXON INTRON
FT[1]=exons
FT[2]=introns




for n in `seq 1 1 2`; do		# genomic feature (exon intron)

   GO_FEATURE=${FT[$n]}

   for k in `seq 1 1 1`; do 			# cycle through ref genomes 

      GO_AA=${AA[$k]}/${SUBFOLDER[$n]}
      GO_RR=${GO_AA} 

      GO_POLYP=${POLYP}

      # NO SLIM
      AUX_prefix_out=""

      #input file suffix
      InFileSUF=GENE_TREES.${GO_FEATURE}_trimN_SLIM-${GO_POLYP}-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.newick

      #output file (root)
      OUTfile=ASTRAL.${GO_FEATURE}_trimN_SLIM-${GO_POLYP}-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}

      rm -f $GO_RR/${OUTfile}-BS20.newick
      rm -f $GO_RR/${OUTfile}-BS20.quartet.t8.newick

      echo
      echo $InFileSUF
      echo $OUTfile
      echo $GO_AA
      echo $GO_RR
      echo

      cd $GO_RR

      #####  Contracting very low support branches
      nw_ed  $GO_AA/$InFileSUF 'i & b<=30' o > $GO_AA/${InFileSUF%.newick}-BS30.newick

      ##### run ASTRAL
      java -jar $ASTRAL \
              -i $GO_AA/${InFileSUF%.newick}-BS30.newick \
              -o $GO_RR/${OUTfile}-BS30.newick \
              -x

      ##### Alternative: provide quartet support (-t 8): Outputs q1, q2, q3.
      java -jar $ASTRAL \
              -i $GO_AA/${InFileSUF%.newick}-BS30.newick \
              -o $GO_RR/${OUTfile}-BS30.quartet.t8.newick \
              -t 8 \
              -x


      # Note: -x	run the exact version of the ASTRAL algorithm (the entire tree space is tested; this is only possible if the number of species is small[<=18])






   done
done


