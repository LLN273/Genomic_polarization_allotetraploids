#!/bin/bash -l
#SBATCH -J splitbylane
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 14:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
#set -eo pipefail


## Takes a fastq file as input and produces a separate fastq for each sequencing lane.
## It is assumed that input fastq files have already been split by sequencing-machine/run/flow-cell, if required.



# load modules
module load bioinfo-tools





#read name of files containing paired reads (with path)
READ1=${1:?msg}                            
READ2=${2:?msg}

#read name of files containing paired reads (no path)
rawReads1=${3:?msg}                         
rawReads2=${4:?msg}

# output folder
RR=${5:?msg}			

#sample name
SNAME=${6:?msg}

#header type
HEADERT=${7:?msg}

# current path
SRCDIR=$(pwd)                                                                    

echo
echo "Sample:" $SNAME
echo "Header type:" $HEADERT
echo "Libraries:"
echo $READ1
echo $READ2






if [ "$HEADERT" = "0" ];then
   # There is NO lane info on header  >> rename libraries and copy to output folder
   # type 0: @SRR1945833.1 SRR492245.1/1
   rsync -ah $READ1 $RR/${SNAME}.LANE-1_R1.paired.fastq.gz
   rsync -ah $READ2 $RR/${SNAME}.LANE-1_R2.paired.fastq.gz

else

   ### unzip input files to scratch disk
   gunzip -c $READ1 > $SNIC_TMP/rawReads1
   gunzip -c $READ2 > $SNIC_TMP/rawReads2


   
   ###  Split Reads For Different Flowcell Lanes

   cd $SNIC_TMP

   for i in `seq 1 1 2`; do 

      rm -f lane.*.fastq

      if [ "$HEADERT" = "1" ]; then
         #### lane is in field 4 (separator ':') 
         #### type 1 >>> @SRR2040797.1 BS-DSFCONTROL03:172:D17EJACXX:2:1101:1419:2084/1
   
         # create fastq file for each lane
         awk 'BEGIN {FS = ":"} {lane=$4 ; print > "LANE-"lane".fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "LANE-"lane".fastq"}}' < $SNIC_TMP/rawReads${i}

      elif [ "$HEADERT" = "2" ]; then
         #### lane is in field 2 (separator ':')
         #### type 2 >>> @SRR2040777.1 HWI-ST1176_0190:6:1101:4596:1991/1

         # create fastq file for each lane
         awk 'BEGIN {FS = ":"} {lane=$2 ; print > "LANE-"lane".fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "LANE-"lane".fastq"}}' < $SNIC_TMP/rawReads${i}

      fi

      # get list fastq files created
      ls LANE-*.fastq > lane.list

      # copy files to results folder
      while read -r LINE                                                                      
      do
         gzip -c $LINE > $RR/${SNAME}.${LINE%.fastq}_R${i}.paired.fastq.gz
      done < lane.list


   

   done

fi



rm -f $SNIC_TMP/rawReads*

