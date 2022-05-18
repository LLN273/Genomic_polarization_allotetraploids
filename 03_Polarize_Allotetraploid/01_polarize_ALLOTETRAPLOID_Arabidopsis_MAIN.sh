#!/bin/bash
#

## Script used to polarize allotetraploid sequence in the MSA.
## Polarization done separately for each locus.




#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# Input folder
AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P03a_MS_alignment_IUPAC/BWA_unmasked_Alyrata
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P03a_MS_alignment_IUPAC/BWA_unmasked_Ahalleri-ensembl


# output folder (output folder name must be changed if POLYP or REFSEQ is changed)
RR[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P04a_Polarize_IUPAC_ALT1_TETRA/BWA_unmasked_Alyrata
RR[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P04a_Polarize_IUPAC_ALT1_TETRA/BWA_unmasked_Ahalleri-ensembl



# Polyploid ID  [NOTE: output folder name (RR) must be changed if POLYP is changed]
POLYP="AKAMCHATICAKWS_DRR054581"    # allo 4X
#POLYP="ASUECICAASO5_SRR2084157"     # allo 4X


# Exclude (if MSA contains more than one allopolyploid, exclude species not selected as POLYP)
EXCLUDE="ASUECICAASO5_SRR2084157"
#EXCLUDE="AKAMCHATICAKWS_DRR054581"


# initial reference sequence ID [can be any sequence in the MSA] [NOTE: output folder name (RR) must be changed if REFSEQ is changed]
REFSEQ="ALYRATAPETRAEA3_SRR2040797"      # Aly 2X
#REFSEQ="ATHALIANA7102_SRR1945833"      # Ath 2X



# Number of species (including outgroup species)
NSPEC=12



#file containing gene list
GL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list_Alyrata.txt
GL[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list_Ahalleri-ensembl.txt


# subfolders containing input data (exons or introns) 
SUBFOLDER[1]=SINGLE_TRANSCRIPTS
SUBFOLDER[2]=SINGLE_INTRONS


# genomic feature
GF[1]=exons
GF[2]=introns




for n in `seq 1 1 2`; do		# genomic feature (exons or introns)

   GO_FEATURE=${GF[$n]}

   for k in `seq 1 1 1`; do 			# ref genome (mapping)

      GO_AA=${AA[$k]}/${SUBFOLDER[$n]}
      GO_RR=${RR[$k]}/${SUBFOLDER[$n]}
      GO_GL=${GL[$k]}
      GO_NSPEC=${NSPEC}
      mkdir -p $GO_RR

      cd $GO_RR
      rsync -ah $SRCDIR_INI/01s_polarizeTETRA.py .

      GO_POLYP=${POLYP}
      GO_EXCLUDE=${EXCLUDE}
      GO_REF=${REFSEQ}
   
      echo $GO_NSPEC $GO_FEATURE
      echo $GO_POLYP
      echo $GO_EXCLUDE
      echo $GO_REF
      echo $GO_GL
      echo $GO_AA
      echo $GO_RR

      cd $SRCDIR_INI

      sbatch ./01_polarize_ALLOTETRAPLOID_Arabidopsis_query.sh $GO_AA \
                                        $GO_RR \
                                        $GO_NSPEC \
                                        $GO_POLYP \
                                        $GO_REF \
                                        $GO_GL \
                                        $GO_FEATURE \
                                        $GO_EXCLUDE

      sleep 0.1
      echo


   done     
done
