#!/bin/bash
#
#SBATCH -J infogenes
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
 


## Script used to get list of informative loci (based on IQ-TREE2 output log file)
## Filter-out locus if any of the following is satisfied:
## - number of parsimony-informative sites < 100
## - MSA contains sequences that contain more than 50% gaps/ambiguities
## - MSA contains sequences with identical nucleotide composition
## - IQ-TREE2 bootstrap analysis did not converge



# load modules
module load bioinfo-tools





#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

declare -A AA
declare -A RR

################ input folder
AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P05a_phylogeny_IQ-TREE2_Trim_IUPAC_TETRA/BWA_unmasked_Alyrata
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P05a_phylogeny_IQ-TREE2_Trim_IUPAC_TETRA/BWA_unmasked_Ahalleri-ensembl

##### output folder
RR=${SRCDIR_INI}/INFORMATIVE_SITES_TETRA_100


# Allopolyploid ID
POLYP="AKAMCHATICAKWS_DRR054581"	# allo
#POLYP="ASUECICAASO5_SRR2084157"	# allo

# Number of species in MSA
NSPEC=11   

## Reference genome (mapping)
ATYPE[1]=Alyrata_REF
ATYPE[2]=Ahalleri-ensembl_REF

## subfolders containing input data (exons or introns) 
SUBFOLDER[1]=SINGLE_TRANSCRIPTS
SUBFOLDER[2]=SINGLE_INTRONS

#### FEATURE: EXON or INTRON
FT[1]=exons
FT[2]=introns

# Gene list
GL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list_Alyrata.txt
GL[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list_Ahalleri-ensembl.txt

# model name
MODEL_outfilename=MFP_ModelFinder





######## read gene list

declare -A GENE_LIST

i=1

while read -r LINE                                                                      
do
   GENE_LIST[1,$i]=$LINE
   GENE_LIST_1[$i]=$LINE
   let "i+=1"
done < ${GL[1]}

N_GENES[1]=${#GENE_LIST_1[@]}



i=1

while read -r LINE                                                                      
do
   GENE_LIST[2,$i]=$LINE
   GENE_LIST_2[$i]=$LINE
   let "i+=1"
done < ${GL[2]}

N_GENES[2]=${#GENE_LIST_2[@]}










for n in `seq 1 1 2`; do		# genomic feature (exon intron)

   GO_FEATURE=${FT[$n]}

   for j in `seq 1 1 1`; do  		# cycle through ref genomes 

      GO_AA=${AA[$j]}/${SUBFOLDER[$n]}
      GO_RR=${RR}
      GO_ATYPE=${ATYPE[$j]}
      mkdir -p $GO_RR

      GO_POLYP=${POLYP}

      # NO SLIM
      AUX_prefix_out=""

      #output files
      OUTfile=INFORMATIVE_SITES_03.${GO_FEATURE}_${GO_ATYPE}_UFBoot_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}.txt
      OUT_FILE2=00_Arabidopsis_gene_list_${GO_FEATURE}_${GO_ATYPE}_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}_TRIM_MAX-00000.txt
      rm -f ${GO_RR}/$OUTfile
      rm -f ${GO_RR}/$OUT_FILE2

      echo
      echo 
      echo $OUTfile
      echo $GO_AA
      echo $GO_RR

      for i in `seq 1 1 ${N_GENES[$j]}`; do                 # cycle through genes

            GO_GENE=${GENE_LIST[$j,$i]}

            #input file (IQ-TREE2 output log file)
            InFileSUF=${GO_GENE}.${GO_FEATURE}_trimN_SLIM-${GO_POLYP}-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.log

            cd $GO_AA

            # get number of informative sites
            NPS_AUX="$(cat ${InFileSUF} | grep 'parsimony-informative' | cut -f1 -d' ')"
            NPS[i]=$NPS_AUX

            if [[ "${NPS[$i]}" -ge 100 ]]
            then                  
               NPSF[i]=OK
            else
               NPSF[i]=FAIL
            fi

            # get number of sequences in alignment
            NSEQ[i]="$(cat ${InFileSUF} | grep 'Alignment has' | cut -f3 -d' ')"

            if [[ "${NSEQ[$i]}" = ${NSPEC} ]]
            then                  
               NSEQF[i]=OK
            else
               NSEQF[i]=FAIL
            fi

            # check whether there are sequences that contain more than 50% gaps/ambiguity
            GA_AUX="$(grep 'sequences contain more than 50% gaps/ambiguity' ${InFileSUF})"
         
            if [[ "${GA_AUX}" = "" ]]
            then                  
               GA[i]=OK
            else
               GA[i]=FAIL
            fi
 
            # check whether there are sequences with identical nucleotide composition
            INC_AUX="$(grep 'is identical to' ${InFileSUF})"

            if [[ "${INC_AUX}" = "" ]]
            then                  
               INC[i]=OK
            else
               INC[i]=FAIL
            fi

            # check whether bs analysis converged
            BS_AUX="$(grep 'WARNING: bootstrap analysis did not converge. You should rerun with higher number of iterations (-nm option)' ${InFileSUF})"

            if [[ "${BS_AUX}" = "" ]]
            then                  
               BS[i]=OK
            else
               BS[i]=FAIL
            fi

      done

      # save to file
      rm -f $GO_RR/_AUX_${OUTfile}
      for i in `seq 1 1 ${N_GENES[$j]}`; do
        echo ${GENE_LIST[$j,$i]} ${NPS[$i]} ${NPSF[$i]} ${GA[$i]} ${INC[$i]} ${BS[$i]} ${NSEQ[$i]} ${NSEQF[$i]} >> $GO_RR/_AUX_${OUTfile}
      done

      # sort output file
      cat $GO_RR/_AUX_${OUTfile} | sort -k2 -nr > $GO_RR/${OUTfile}
      rm -f $GO_RR/_AUX_${OUTfile}
      
      # list of loci that pass all filtering criteria
      cat $GO_RR/${OUTfile} | grep -v FAIL | cut -f1 -d' ' > $GO_RR/$OUT_FILE2

         

   done
done




