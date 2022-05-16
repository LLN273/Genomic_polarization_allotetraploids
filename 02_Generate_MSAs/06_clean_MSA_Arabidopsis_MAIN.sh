#!/bin/bash
#
#SBATCH -J filterN
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


## Scrip used for removing regions with too many masked positions (N) in the MSAs (site removed if # masked sites > 2)


# load modules
module load bioinfo-tools
module load trimAl/1.4.1




##### Paths and folders

#remember initial path
SRCDIR_INI=$(pwd) 

# Path to folder containing input data
AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P03a_MS_alignment_IUPAC/BWA_unmasked_Alyrata
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P03a_MS_alignment_IUPAC/BWA_unmasked_Ahalleri-ensembl


#file containing gene list
GL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list_Alyrata.txt
GL[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list_Ahalleri-ensembl.txt	


# gap threshold (accept up to 2 masked sites per column)
gapTH=0.833 		# maximum number of masked sites allowed is 2 [2/12=0.166]  >>>>>>> 1-0.166=0.833


# Output file containing stats (root)
OUTF=Glean_06_filter-masked_Arabidopsis




######## read gene list

declare -A GENE_LIST

i=1

while read -r LINE                                                                      
do
   GENE_LIST[1,$i]=$LINE
   GENE_LIST_1[$i]=$LINE
   let "i+=1"
done < ${GL[1]}


### Number of genes
N_GENES[1]=${#GENE_LIST_1[@]}



i=1

while read -r LINE                                                                      
do
   GENE_LIST[2,$i]=$LINE
   GENE_LIST_2[$i]=$LINE
   let "i+=1"
done < ${GL[2]}


### Number of genes
N_GENES[2]=${#GENE_LIST_2[@]}





# subfolders containing input data (exon or intron) 
SUBFOLDER[1]=SINGLE_TRANSCRIPTS
SUBFOLDER[2]=SINGLE_INTRONS

# genomic feature
GF[1]=exons
GF[2]=introns

# Reference genome
REFGEN[1]=Alyrata
REFGEN[2]=Ahalleri-ensembl


cd $SRCDIR_INI



for n in `seq 1 1 2`; do		# genomic feature (exons or introns)

   GO_FEATURE=${GF[$n]}

   for k in `seq 1 1 2`; do 						# cycle through ref genomes
   
      GO_AA=${AA[$k]}/${SUBFOLDER[$n]}
      GO_RR=${GO_AA}

      echo $GO_AA
      echo $GO_RR
      echo

      OUTF_GO=${OUTF}_${GO_FEATURE}_${REFGEN[$k]}.txt
      rm -f $SRCDIR_INI/$OUTF_GO
      echo "gene" "Masked_before" "Masked_after" "Percentage_masked_before" "Percentage_masked_after"  > $SRCDIR_INI/$OUTF_GO
			
      for i in `seq 1 1 ${N_GENES[$k]}`; do                 # cycle through genes

         GO_GENE=${GENE_LIST[$k,$i]}
	       REF_FASTA=${GO_GENE}.${GO_FEATURE}.fa
               
 	       # convert sequence to big caps
         tr a-z A-Z < $GO_AA/$REF_FASTA > $SNIC_TMP/__${GO_GENE}_AUX_1.fa

         # in sequence headers, convert '-' to '_'
         sed '/^>/ s/-/_/g' $SNIC_TMP/__${GO_GENE}_AUX_1.fa > $SNIC_TMP/__${GO_GENE}_AUX_2.fa

         # convert gaps (-) to Z
         sed 's/-/Z/g' $SNIC_TMP/__${GO_GENE}_AUX_2.fa > $SNIC_TMP/__${GO_GENE}_AUX_3.fa
	 
         # convert masked nucleotides (N) to gaps (-)
         sed 's/N/-/g' $SNIC_TMP/__${GO_GENE}_AUX_3.fa > $SNIC_TMP/__${GO_GENE}_AUX_4.fa

         # remove sites with too many missing sites (use trimAl, with missing sites coded as gaps)
         trimal -in $SNIC_TMP/__${GO_GENE}_AUX_4.fa \
                -out $SNIC_TMP/__${GO_GENE}_AUX_5.fa \
                -gapthreshold $gapTH \
                -cons 0

        # Notes:
	      # -gapthreshold <n>    1 - (fraction of sequences with a gap allowed). Range: [0 - 1]
        # -gapthreshold 0.904 >>>>> maximum number of masked sites allowed is 2 [2/21=0.096]  >>>>>>> 1-0.096=0.904
        # -gapthreshold 0.93 >>>>> maximum number of masked sites allowed is 2 [2/30=0.0667]  >>>>>>> 1-0.07=0.93
        # -cons <n>                Minimum percentage of the positions in the original alignment to conserve. Range: [0 - 100]

         # convert gaps (-) to masked nucleotides (N)
         sed 's/-/N/g' $SNIC_TMP/__${GO_GENE}_AUX_5.fa > $SNIC_TMP/__${GO_GENE}_AUX_6.fa

         # convert Z sites to gaps (-)
         sed 's/Z/-/g' $SNIC_TMP/__${GO_GENE}_AUX_6.fa > $SNIC_TMP/__${GO_GENE}_AUX_7.fa
	   	   
	       # convert sequence to big caps
         tr a-z A-Z < $SNIC_TMP/__${GO_GENE}_AUX_7.fa > $SNIC_TMP/__${GO_GENE}_AUX_8.fa

         # convert to single line fasta
         awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $SNIC_TMP/__${GO_GENE}_AUX_8.fa > $SNIC_TMP/__${GO_GENE}_AUX_9.fa
         cat $SNIC_TMP/__${GO_GENE}_AUX_9.fa | tail -n +2 > $GO_RR/${REF_FASTA%.fa}_trimN.fa


         # STATS

         # Number of masked sites before trimming
         N_before="$(cat $GO_AA/$REF_FASTA | grep -v '^>' | grep -io [N] | wc -l )"

         # Number of masked sites after trimming
         N_after="$(cat $GO_RR/${REF_FASTA%.fa}_trimN.fa | grep -v '^>' | grep -io [N] | wc -l )"

 	       # Total number of sites before trimming
         T_before="$(cat $GO_AA/$REF_FASTA | grep -v '^>' | grep -io [A-T] | wc -l )"

         # Total number of sites after trimming
         T_after="$(cat $GO_RR/${REF_FASTA%.fa}_trimN.fa | grep -v '^>' | grep -io [A-T] | wc -l )"
		   
         # Percentage masked, before trimming
         aux_100="100"   
         if [ -z "$N_before" ]; then
            PERCT_before="0"
         else
            PERCT_before="$(awk '{print $1*$2/$3}' <<<"$aux_100 $N_before $T_before")"		   
         fi

         # Percentage masked, after trimming   
         if [ -z "$N_after" ]; then
            PERCT_after="0"
         else
            PERCT_after="$(awk '{print $1*$2/$3}' <<<"$aux_100 $N_after $T_after")"
         fi	

         echo $GO_GENE $N_before $N_after $PERCT_before $PERCT_after >> $SRCDIR_INI/$OUTF_GO
		   
	       # remove auxiliary files
         rm -f $SNIC_TMP/__*.fa	   


      done
   done
done



echo
echo done!
echo


