#!/bin/bash
#
#SBATCH -J INDELIble
#SBATCH -p devcore 
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail


## Script used to create fasta sequences (simulated MSA) with INDELIble, based on simulated gene phylogenies.
## For wrapper's usage, see simphy manual: https://github.com/adamallo/SimPhy/wiki/Manual#521-input-files-newick-tree-format


# load modules
module load bioinfo-tools
module load perl/5.26.2

# path to INDELIble
INDELIBLE="/crex1/proj/snic2017-7-149/private/Luis/z_APPS/INDELible/INDELibleV1.03/src/"
export PATH="$INDELIBLE:$PATH"

# simphy's INDELIble wrapper
IDW=/crex1/proj/snic2017-7-149/private/Luis/z_APPS/SimPhy/SimPhy_1.0.2/scripts/INDELIble_wrapper.pl


#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# Input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/20_No-ILS-No-polytomy_Simulated_Dataset/X_NoILS_16t_CLEAN


# Configuration file
CONFIG=/crex1/proj/snic2017-7-149/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/10_No-ILS-No-polytomy_Simulated_Dataset/INDELible_simple.txt


# Run INDELIble (perl wrapper)
# Run as: perl INDELIble_wrapper.pl directory input_config seed numberofcores

perl $IDW \
     $AA \
     $CONFIG \
     22 \
     1
     
