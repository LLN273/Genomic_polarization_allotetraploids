#!/bin/bash -l
#SBATCH -J OrthoFinder
#SBATCH -p devcore 
##SBATCH -p node 
#SBATCH -n 20
#SBATCH -t 01:00:00
#SBATCH -A snic2021-22-291
##SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL



## Script used to identify single copy orthogroups with ORTHOFINDER.


# load modules
module load bioinfo-tools
module load mcl/14-137


# local python installation
export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
pyenv global 2.7.13





##### Paths and folders

# Path to file with list of proteomes
# Note: Provide amino acid sequences, in FASTA format, for the species you want to analyze. If you
# have the option, it is best to use a version containing a single representative/longest transcript-variant for each gene.
AA=/crex1/proj/snic2017-7-149/private/Luis/P1_DemSim/12_OrthoFinder_Arabidopsis/00_proteome_list_primaryTranscript_Alyrata.txt
#AA=/crex1/proj/snic2017-7-149/private/Luis/P1_DemSim/12_OrthoFinder_Arabidopsis/00_proteome_list_primaryTranscript_halleri_ensembl.txt


# Path to output folder (root output folder)
RR=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P1_DemSim/20_OrthoFinder

# Proteome folder (to be created in results folder)
PTF=proteome_PrimaryTranscripsMix_Alyrata
#PTF=proteome_PrimaryTranscripsMix_Ahalleri_ensembl

# OrthoFinder results folder (to be created in results folder)
RRO=diamond_transcriptMix_Alyrata
#RRO=diamond_transcriptMix_Ahalleri_ensembl

# Path to OrthoFinder
ORTHOFINDER=/home/luisleal/MYAPPS/OrthoFinder-2.3.1/orthofinder
#ORTHOFINDER=/home/luisleal/MYAPPS/OrthoFinder/orthofinder/orthofinder.py


# path to FastMe
# Note: must rename one of binaries in this folder to 'fastme' (Uppmax:    cp fastme-2.1.5-linux64 fastme)
export PATH=$PATH:/home/luisleal/MYAPPS/fastme-2.1.5/binaries

# path to diamond (version installed on Uppamx is old)
export PATH=$PATH:/home/luisleal/MYAPPS/diamond_v0.9.22

# path to MMseq2 folder
export PATH=$PATH:/home/luisleal/MYAPPS/mmseqs2/bin

# path to DLCpar
export PATH=$PATH:/home/luisleal/MYAPPS/dlcpar-1.0/bin

#remember initial path
SRCDIR_INI=$(pwd)   








######## read proteome list

i=1

while read -r LINE                                                                      
do
   PROTEOME_LIST[$i]=$LINE
   let "i+=1"
done < ${AA}



### Number of proteomes
N_PROT=${#PROTEOME_LIST[@]}





# check whether OUTPUT folder exists
cd $RR
if [ -d "$RRO" ]                            #check whether output folder already exists
   then
      { echo "Error: $RRO folder already exists."; echo ; exit 1; }
fi





##### copy proteomes to results folder

cd $RR

#check whether PROTEOME folder already exists
if [ -d "$PTF" ]; then
   echo; echo "Warning: $PTF folder already exists."; echo 
else
   mkdir $PTF
   cd $PTF
   
   for i in `seq 1 1 $N_PROT`; do 
      rsync -ah ${PROTEOME_LIST[$i]} .
   done
fi






# run OrthoFinder

cd $RR


# when using Diamond (fastest method)
$ORTHOFINDER -f $RR/$PTF \
	     -t $SLURM_NTASKS \
             -a $SLURM_NTASKS \
	     -S diamond \
             -o $RR/$RRO



# when using MMseqs2 (compiled version not working; use source python script) [MMseqs2 is slower than Diamond, about twice as slow]
#$ORTHOFINDER -f $RR/$PTF \
#	      -t $SLURM_NTASKS \
#             -a $SLURM_NTASKS \
#	      -S mmseqs

##             -o $RR/$RRO   # -o option not valid when using python script   >> results saved in genomes folder when using python script



# when using blast (very, very slow)
#$ORTHOFINDER -f $RR/$PTF \
#             -t $SLURM_NTASKS \
#             -a $SLURM_NTASKS \
#	      -S blast \
#             -o $RR/$RRO 



# run OrthoFinder (test dataset)

#$ORTHOFINDER -f /home/luisleal/MYAPPS/OrthoFinder-2.3.1/ExampleData \
#	      -S blast \
#             -o $RR 


# OPTIONS:
# -t <int>          Number of parallel sequence search threads [Default = 40]
# -a <int>          Number of parallel analysis threads [Default = 1]
# -M <txt>          Method for gene tree inference. Options 'dendroblast' & 'msa'
#                   [Default = dendroblast]
# -S <txt>          Sequence search program [Default = diamond]
#                   Options: blast, mmseqs, blast_gz, diamond
# -A <txt>          MSA program, requires '-M msa' [Default = mafft]
#                   Options: muscle, mafft
# -T <txt>          Tree inference method, requires '-M msa' [Default = fasttree]
#                   Options: iqtree, raxml-ng, fasttree, raxml
# -s <file>         User-specified rooted species tree
# -I <int>          MCL inflation parameter [Default = 1.5]
# -x <file>         Info for outputting results in OrthoXML format
# -p <dir>          Write the temporary pickle files to <dir>
# -1                Only perform one-way sequence search 
# -n <txt>          Name to append to the results directory
# -o <txt>          Non-default results directory






