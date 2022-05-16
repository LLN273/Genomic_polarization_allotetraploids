#!/bin/bash
#
#SBATCH -J MSA
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


# Script used to generate a MSA for each locus (based on exonic sequences only).



# load modules
module load bioinfo-tools



#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

########### folder containing input fasta files
AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P02a_fetch_gene_sequences_IUPAC_SG/BWA_unmasked_Alyrata		        # Based on Alyrata reference genome   
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P02a_fetch_gene_sequences_IUPAC_SG/BWA_unmasked_Ahalleri-ensembl	# Based on Ahalleri-ensembl reference genome
    

###########  output folder
RR[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P03a_MS_alignment_IUPAC/BWA_unmasked_Alyrata/SINGLE_TRANSCRIPTS
RR[2]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P03a_MS_alignment_IUPAC/BWA_unmasked_Ahalleri-ensembl/SINGLE_TRANSCRIPTS


# Sample list
SPL=/crex1/proj/snic2017-7-149/private/Luis/P1_DemSim/08_VCFgatk3.8/00_samples_list_Arabidopsis.txt

			               			
#file containing gene list (or transcript list)
GL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list_Alyrata.txt		
GL[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list_Ahalleri-ensembl.txt
			               			



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






######## read gene list

declare -A GENE_LIST

i=1

while read -r LINE                                                                      
do
   GENE_LIST[1,$i]=$LINE
   GENE_LIST_1[$i]=$LINE
   let "i+=1"
done < ${GL[1]}


### Number of samples
N_GENES[1]=${#GENE_LIST_1[@]}




i=1

while read -r LINE                                                                      
do
   GENE_LIST[2,$i]=$LINE
   GENE_LIST_2[$i]=$LINE
   let "i+=1"
done < ${GL[2]}


### Number of samples
N_GENES[2]=${#GENE_LIST_2[@]}





##reads each filename at a time, stores names in array
declare -A READS  

for i in `seq 1 1 ${N_SAMPLES}`; do       #loop through samples     
      READS[1,$i]=${SAMPLE_LIST[1,$i]}.exons.concatenated.fa
done








########### Generate a MSA for each locus

cd $RR

echo 'Create fasta file for each gene'
for k in `seq 1 1 2`; do 						# loop through ref genomes 

   GO_AA=${AA[$k]}
   GO_RR=${RR[$k]}   
   mkdir -p $GO_RR  
   
   for n in `seq 1 1 ${N_GENES[$k]}`; do 				# loop through gene list
   
      GO_GENE=${GENE_LIST[$k,$n]}
      GO_OUTFILE=${GO_GENE}.exons.fa
      rm -f $GO_RR/$GO_OUTFILE
   
         for i in `seq 1 1 ${N_SAMPLES}`; do  			# loop through samples
   
	          DNAFILE_exons=${READS[1,$i]} 
            GO_SAMPLE=${SAMPLE_LIST[1,$i]}
         
	          echo $GO_GENE
	          echo $GO_SAMPLE
	          echo $DNAFILE_exons      
            echo $GO_AA
            echo $GO_RR
   
            # paste sequences associated to gene into a single fasta file
             header_AUX=${GO_SAMPLE}
	           body_AUX[1]="$(grep -A 1 $GO_GENE $GO_AA/$DNAFILE_exons | grep -v "^>")"
	  
	           echo ">"${header_AUX} >> $GO_RR/$GO_OUTFILE 
	           echo ${body_AUX[1]} >> $GO_RR/$GO_OUTFILE
             echo

         done
   done
done




