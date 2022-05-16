#!/bin/bash -l
#SBATCH -J Ortho_genes
#SBATCH -p devcore 
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2021-22-291
##SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited

## Script used to get list of genes in single copy orthogroups (obtained using OrthoFinder/diamond)


# load modules
module load bioinfo-tools


###### START UPPMAX RUNS 

echo
echo "Starting Uppmax jobs ..."
echo





##### Paths and folders

###### remember initial path
SRCDIR_INI=$(pwd)   

###### Path to OrthoFinder results (Orthogroups_SingleCopyOrthologues.txt and Orthogroups.tsv) (Note: orthogroup codes are NOT shared between methods)

## Alyrata
echo
echo "diamond" 
AA1=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P1_DemSim/20_OrthoFinder/diamond_transcriptMix_Alyrata/Results_Dec19/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
BB1=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P1_DemSim/20_OrthoFinder/diamond_transcriptMix_Alyrata/Results_Dec19/Orthogroups/Orthogroups.tsv
gene_extension=v2.1
echo $AA1
echo $BB1

## Ahalleri_ensembl
#echo
#echo "diamond" 
#AA1=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P1_DemSim/20_OrthoFinder/diamond_transcriptMix_Ahalleri_ensembl/Results_Apr19/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
#BB1=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P1_DemSim/20_OrthoFinder/diamond_transcriptMix_Ahalleri_ensembl/Results_Apr19/Orthogroups/Orthogroups.tsv
#gene_extension=""
#echo $AA1
#echo $BB1


###### Output files
## Alyrata
SCD=singleCopyOrthologues_diamond_Alyrata.txt
## Ahalleri_ensemb
#SCD=singleCopyOrthologues_diamond_Ahalleri_ensemb.txt



###### Path to output folder 
RR=$SRCDIR_INI









##### read list of orthogroup codes associated to sigle copy orthogroups

i=1

while read -r LINE                                                                      
do
   DIAMOND_SCO_LIST[$i]=$LINE
   let "i+=1"
done < ${AA1}

# number of sigle copy orthogroups
N_DSCO=${#DIAMOND_SCO_LIST[@]}
echo 
echo "Diamond | Number of single-copy orthogroups:" $N_DSCO




##### fetch list of genes associated to each orthogroup

cd $RR
rm -f $SCD


for i in `seq 1 1 $N_DSCO`; do 

   header_diamond="$(head -1 ${BB1})"
   aux_d1="$(grep ${DIAMOND_SCO_LIST[$i]} ${BB1} | cut -f2 | cut -f1 -d' ' | cut -f1 -d'.')"
   control_aux="$(grep ${DIAMOND_SCO_LIST[$i]} ${BB1} | cut -f 2 | wc -l)"
   if [ "$control_aux" != "1" ];then
      echo "ERROR: found" $control_aux "entries in diamond's Orthogroups.tsv file when searching orthogroup with code" ${LINE} ; exit 1
   fi

   if [ "$SCD" == "singleCopyOrthologues_diamond_Ahalleri_ensemb.txt" ];then
      echo ${aux_d1} >> $RR/$SCD
   else
      echo ${aux_d1}.${gene_extension} >> $RR/$SCD
   fi

done





#### Keep unique entries only

# number of entries per gene
cat $RR/$SCD | sort | uniq > $SRCDIR_INI/__aux1.txt
# remove spaces at beginning of line
sed 's/^[[:space:]]*//g' $SRCDIR_INI/__aux1.txt > $SRCDIR_INI/__aux2.txt
# keep genes with a single entry
awk '$1 == "1"' $SRCDIR_INI/__aux2.txt > $SRCDIR_INI/__aux3.txt
cat $SRCDIR_INI/__aux3.txt | cut -f2 -d' ' > $SRCDIR_INI/__aux4.txt
mv $SRCDIR_INI/__aux4.txt $RR/$SCD





