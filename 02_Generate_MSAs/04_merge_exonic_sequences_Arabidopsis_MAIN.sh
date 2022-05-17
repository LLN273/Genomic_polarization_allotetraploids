#!/bin/bash
#
#SBATCH -J mergeExome
#SBATCH -p devcore 
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=END
ulimit -c unlimited


## Scrip used to concatenate exonic (or intragenic) sequences associated to the same gene.
## Performed independently for each sample.


# load modules
module load bioinfo-tools




##### Paths and folders

#remember initial path
SRCDIR_INI=$(pwd) 

# Path to folder containing input data
AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P02a_fetch_gene_sequences_IUPAC_SG/BWA_unmasked_Alyrata		        # Based on Alyrata reference genome
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P02a_fetch_gene_sequences_IUPAC_SG/BWA_unmasked_Ahalleri-ensembl	# based on Ahalleri-ensembl reference genome


# Sample list
SPL=/crex1/proj/snic2017-7-149/private/Luis/P1_DemSim/08_VCFgatk3.8/00_samples_list_Arabidopsis.txt





######## read sample list

declare -A SAMPLE_LIST

i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[1,$i]=$LINE
   SAMPLE_LIST_1[$i]=$LINE
   let "i+=1"
done < ${SPL}


### Number of samples
N_SAMPLES=${#SAMPLE_LIST_1[@]}




# Genomic feature
Gfeature[1]=exons
Gfeature[2]=introns


cd $SRCDIR_INI


for n in `seq 1 1 2`; do 						# loop through genomic features (exons, introns)

   rm -f $SRCDIR_INI/slurm_merge_${Gfeature[$n]}_sequences_Alyrata.txt
   rm -f $SRCDIR_INI/slurm_merge_${Gfeature[$n]}_sequences_Ahalleri-ensembl.txt

   for k in `seq 1 1 2`; do 						# cycle through ref genomes
      for i in `seq 1 1 ${N_SAMPLES}`; do                 		# cycle through samples

         GO_AA=${AA[$k]}
         GO_RR=$GO_AA     
         DNAFILE_A=${SAMPLE_LIST[1,$i]}.${Gfeature[$n]}.fa
	 GO_SAMPLE=${SAMPLE_LIST[1,$i]}  
	     
         echo $GO_SAMPLE   
         echo $DNAFILE_A   
         echo $GO_AA
         echo $GO_RR
                                           


	 # Order fasta headers alphabetically 
         perl -pe 's/[\r\n]+/;/g; s/>/\n>/g' $GO_AA/$DNAFILE_A | sort -t"[" -k2,2V | sed 's/;/\n/g' | sed '/^$/d' > $GO_RR/${DNAFILE_A%.fa}_aux1.fa

		 
         # reformat fasta headers (keep only gene name)
         cat $GO_RR/${DNAFILE_A%.fa}_aux1.fa | awk -F"." '/>/{$0=">"$1}1' > $GO_RR/${DNAFILE_A%.fa}_aux2.fa


         #### concatenate exon sequences associated to the same gene
         awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' $GO_RR/${DNAFILE_A%.fa}_aux2.fa > $GO_RR/_AUX_${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa
         cat $GO_RR/_AUX_${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa | tail -n +2 > $GO_RR/${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa

         rm -f $GO_RR/_AUX_${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa
         rm -f $GO_RR/${DNAFILE_A%.fa}_aux*.fa

         # Check number of nucleotides and number of headers in fasta file
         aux_1="$(cat $GO_AA/$DNAFILE_A | grep -v '^>' | grep -io [A-Z] | wc -l  )"
         aux_1_concat="$(cat $GO_RR/${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa | grep -v '^>' | grep -io [A-Z] | wc -l  )"
         aux_headers="$(cat $GO_RR/${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa | grep '^>' | wc -l  )"
         
         if [[ $k = "1" ]] ; then
            echo $GO_RR/${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa $aux_1 $aux_1_concat $aux_headers >> $SRCDIR_INI/slurm_merge_${Gfeature[$n]}_sequences_Alyrata.txt
         elif [[ $k = "2" ]] ; then
            echo $GO_RR/${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa $aux_1 $aux_1_concat $aux_headers >> $SRCDIR_INI/slurm_merge_${Gfeature[$n]}_sequences_Ahalleri-ensembl.txt
         fi

         echo
      done
   done
done


echo "done!"
echo



