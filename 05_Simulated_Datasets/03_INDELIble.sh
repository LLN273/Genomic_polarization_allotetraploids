#!/bin/bash
#
#SBATCH -J INDELIble
#SBATCH -p core 
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -A snic2020-15-40
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=END
ulimit -c unlimited
set -eo pipefail


## Script used to create fasta sequences with INDELIble, based on simulated gene phylogenies.
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
AA=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/01_simphy_ILS/A_modILS_recSpeciation_16t_CLEAN
#AA=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/01_simphy_ILS/B_modILS_deepSpeciation_16t_CLEAN
#AA=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/01_simphy_ILS/C_HIGH-ILS_recSpeciation_16t_CLEAN
#AA=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/01_simphy_ILS/D_HIGH-ILS_deepSpeciation_16t_CLEAN


# Configuration file
CONFIG=/crex1/proj/snic2017-7-149/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/01_simphy_ILS/INDELible_complex.txt


# Run INDELIble (perl wrapper)
# Run as: perl INDELIble_wrapper.pl directory input_config seed numberofcores

perl $IDW \
     $AA \
     $CONFIG \
     22 \
     $SLURM_NTASKS
     


