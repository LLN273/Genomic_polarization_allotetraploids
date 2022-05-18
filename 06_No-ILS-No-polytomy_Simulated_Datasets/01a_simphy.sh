#!/bin/bash
#
#SBATCH -J simphy
#SBATCH -p devcore 
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2017-7-149
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=END
ulimit -c unlimited



## Script used to generate simulated gene phylogenies with SimPhy [No ILS scenario]
## For each reference tree (created randomly), a single gene tree was created. This process was repeated for 100 replicates.

## simphy manual: https://github.com/adamallo/SimPhy/wiki/Manual#521-input-files-newick-tree-format


# load modules
module load bioinfo-tools

# path to simphy
simphy=/crex1/proj/snic2017-7-149/private/Luis/z_APPS/SimPhy/SimPhy_1.0.2/bin/simphy_lnx64



#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/01_simphy_ILS






##### Conditions

# General conditions:
# 100 replicates (-rs)
# 1 locus (-rl)
# 1 gene tree per locus (-rg)
# 2 individuals per species (-si)
# Ne=400K (-sp)


# no ILS, deep speciation, 16 taxa
b=0.0000001
OUTF=X_NoILS_16t
ntaxa=16






## Run simphy (No ILS scenario)

mkdir -p $RR/$OUTF

$simphy -rs 100 \
        -rl f:1 \
        -rg 1 \
        -sb f:${b} \
        -sd f:0 \
        -st ln:14.70055,0.25 \
        -sl f:${ntaxa} \
        -so f:1 \
        -si f:2 \
        -sp f:400000 \
        -su ln:-17.27461,0.6931472 \
        -hh f:1 \
        -hs ln:1.5,1 \
        -hl ln:1.551533,0.6931472 \
        -hg ln:1.5,1 \
        -cs 55277 \
        -v 3 \
        -o $RR/$OUTF \
        -ot 0 \
        -op 1 \
        -od 1



#### Notes:

## Replicates
# -RS i: 	Number of species tree replicates (i.e., study replicates)	[50 replicates (simulation is repeated 50 times, each time with a different species tree)]
# -RL *: 	Number of locus trees per species tree.				[1200 locus computed for each species tree, that is, for each replicate]
# -RG i: 	Number of gene trees per locus tree. (Not for general usage).	[1]  (each locus contains one gene)

## Species tree parameters
# -S		Fixed species tree 			For example: (((A:10000,B:10000):30000,C:40000):1000,D:41000); 
# -SB *: 	Speciation rate (events/time unit).	[birth rate]
# -SD *: 	Extinction rate (events/time unit).	[0] (no extinction events (birth-only process))
# -ST *: 	Species tree height (time units).
# -SL *: 	Number of taxa.				[16 + root]
# -SO *: 	Ratio between ingroup height and the branch from the root to the ingroup. If this parameter is not set the outgroup is not simulated. [1]
# -SI *: 	Number of individuals per species.	[2]
# -SP *: 	Tree-wide effective population size.	[400k] (effective population size, fixed (same for all species))
# -SU *: 	Tree-wide substitution rate.

##  Substitution rate heterogeneity parameters
# -HH *: 	Gene-by-lineage-specific locus tree parameter (to use with the HG argument below) [1]
# -HS *: 	Species-specific branch rate heterogeneity modifiers.
# -HL *: 	Gene-family-specific rate heterogeneity modifiers.
# -HG *: 	Gene-by-lineage-specific rate heterogeneity modifiers.

## Global options
# -CS i: 	Random number generator seed.

## Verbosity
#    -V 0: Only warnings and errors.
#    -V 1: Global settings summary, simulation progress per replicate (number of simulated gene trees), warnings and errors.
#    -V 2: Global settings summary, simulation progress per gene tree, warnings and errors.
#    -V 3: Global settings summary, sampled settings per species , locus and gene trees, simulation progress per gene tree, warnings and errors.
#    -V 4: Global settings summary, sampled settings per species, locus and gene trees, simulation progress per gene tree, simulated trees, warnings and errors.
#    -V 5: Global settings summary, sampled settings per species, locus and gene trees, simulation progress detailing the internal algorithms (events), simulated trees, IO progress, warnings and errors.
#    -V 6: Global settings summary, sampled settings per species, locus and gene trees, simulation progress detailing the internal algorithms with the highest detail (event times, rejections, probabilities…), simulated trees, IO progress, warnings and errors.

## Output parameters
# -O c: 	Common output prefix-name (for folders and files).
# -OT b: 	Determines whether the species and locus tree branches are written in number of generations (0) or time units (1).	[0]
# -OP b: 	Activates logging of sampled options.
# -OD b: 	Activates the SQLite database output.

# Distribution parameters:
# 			Code 	Parameter 1 	Parameter 2 	Parameter 3
# Point (fixed value) 	F 	Value 		
# Uniform 		U 	Minimum 	Maximum 	
# Normal 		N 	Mean 		Scale 	
# Exponential 		E 	Rate 		
# Gamma 		G 	Shape 		Scale 	
# Lognormal 		LN 	Location 	Scale 	
# LogUniform 		LU 	Minimum 	Maximum 	
# Lognormal * constant 	SL 	Location 	Scale 		Multiplier
