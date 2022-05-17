#!/bin/bash
#
#SBATCH -J IQT
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 250:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



## Script used to perform phylogenetic inference for each locus using IQ-TREE2.


# load modules
module load bioinfo-tools
module load iqtree/2.0-rc2-omp-mpi





#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

#input folder
AA=${1:?msg}	

#input file suffix (eg "-ALT.fasta" if MSAs are named "XXXX-ALT.fasta", where "XXXX" is the locus name)
InFileSUF=${2:?msg}	

#output folder
RR=${3:?msg}

#number of bootstraps (eg 1000)
BS=${4:?msg}

# output file extension
OFE=${5:?msg}

# substitution model (use 'MFP' if using ModelFinder to determine best substitution model)
MODEL_IQTREE2=${6:?msg}

#gene list
GL=${7:?msg}	

# Outgroup sequence ID (eg "CRUBELLA1_SRR2082718")
OUTG=${8:?msg}


echo
echo "input file suffix:" $InFileSUF
echo "Input folder:" $AA
echo 'Gene list:' $GL
echo 
echo "Number of bootstraps replicates:" $BS
echo 'Substitution model:' $MODEL_IQTREE2
echo 'Outgroup:' $OUTG




######## read gene list

i=1

while read -r LINE                                                                      
do
   GENE_LIST[$i]=$LINE
   let "i+=1"
done < ${GL}

N_GENES=${#GENE_LIST[@]}
echo "Total number of genes:" $N_GENES





########### run IQ-TREE2 with ultrafast bootstrap approximation 

cd $RR

for i in `seq 1 1 $N_GENES`; do                 # cycle through genes

   GO_GENE=${GENE_LIST[$i]}                      # locus ID
   GO_MSA=${GO_GENE}${InFileSUF}                 # MSA filename
   AUX_prefix_out=${GO_MSA%.fa*}_${OFE}          # prefix for all output files

   rm -f $RR/${AUX_prefix_out}*

   echo
   echo $i $GO_MSA

   iqtree-omp --ufboot $BS \
              -alrt $BS \
              -s $AA/$GO_MSA \
              --seqtype DNA \
              -m ${MODEL_IQTREE2} \
              -nt $SLURM_NTASKS \
              --seed 4545 \
              --prefix ${AUX_prefix_out} \
              --redo \
              --safe \
              --allnni \
	            -o $OUTG \
              -bnni


done




##### Notes:
#
## GENERAL OPTIONS:
# -s 		FILE[,...,FILE]   	PHYLIP/FASTA/NEXUS/CLUSTAL/MSF alignment file(s)
# --seqtype 	STRING      		BIN, DNA, AA, NT2AA, CODON, MORPH (default: auto-detect)
# -o 		TAX[,...,TAX]     	Outgroup taxon (list) for writing .treefile
# --prefix 	STRING      		Prefix for all output files (default: aln/partition)
# --seed 	NUM           		Random seed number, normally used for debugging purpose
# --safe               			Safe likelihood kernel to avoid numerical underflow
# --runs 	NUM           		Number of indepedent runs (default: 1)
# --redo               			Ignore checkpoint and overwrite outputs (default: OFF)
# --threads-max NUM    			Max number of threads for -T AUTO (default: all cores)

# --ninit NUM          Number of initial parsimony trees (default: 100)
# --nstop NUM          Number of unsuccessful iterations to stop (default: 100)
# --nbest NUM          Number of best trees retained during search (defaut: 5)
# --allnni             Perform more thorough nearest neighbor interchange (NNI) search (default: OFF)

# NON-PARAMETRIC BOOTSTRAP/JACKKNIFE:   					>>>> SLOW
#  -b, --boot NUM       Replicates for bootstrap + ML tree + consensus tree
#  -j, --jack NUM       Replicates for jackknife + ML tree + consensus tree
#  --jack-prop NUM      Subsampling proportion for jackknife (default: 0.5)
#  --bcon NUM           Replicates for bootstrap + consensus tree
#  --bonly NUM          Replicates for bootstrap only
#  --tbe                Transfer bootstrap expectation

#ULTRAFAST BOOTSTRAP/JACKKNIFE:							>>>> FAST
#  -B, --ufboot NUM     Replicates for ultrafast bootstrap (>=1000)
#  -J, --ufjack NUM     Replicates for ultrafast jackknife (>=1000)
#  --jack-prop NUM      Subsampling proportion for jackknife (default: 0.5)
#  --sampling STRING    GENE|GENESITE resampling for partitions (default: SITE)
#  --boot-trees         Write bootstrap trees to .ufboot file (default: none)
#  --wbtl               Like --boot-trees but also writing branch lengths
#  --nmax NUM           Maximum number of iterations (default: 1000)
#  --nstep NUM          Iterations for UFBoot stopping rule (default: 100)
#  --bcor NUM           Minimum correlation coefficient (default: 0.99)
#  --beps NUM           RELL epsilon to break tie (default: 0.5)
#  --bnni               Optimize UFBoot trees by NNI on bootstrap alignment
#                       In order to reduce the impact of severe model violations while using UFBoot, starting with IQ-TREE version 1.6 we provide a new option -bnni to reduce 
#			the risk of overestimating branch supports with UFBoot due to severe model violations. With this option UFBoot will further optimize each bootstrap tree 
#			using a hill-climbing nearest neighbor interchange (NNI) search based directly on the corresponding bootstrap alignment.
#			If severe model violations are present in the data set at hand, users are advised to append -bnni to the regular UFBoot command(-B).

# SINGLE BRANCH TEST:
#  --alrt NUM           Replicates for SH approximate likelihood ratio test
#  --alrt 0             Parametric aLRT test (Anisimova and Gascuel 2006)
#  --abayes             approximate Bayes test (Anisimova et al. 2011)
#  --lbp NUM            Replicates for fast local bootstrap probabilities


## PARTITION MODEL:
#  -p FILE|DIR          NEXUS/RAxML partition file or directory with alignments
#                       Edge-linked proportional partition model
#  -q FILE|DIR          Like -p but edge-linked equal partition model 
#  -Q FILE|DIR          Like -p but edge-unlinked partition model
#  -S FILE|DIR          Like -p but separate tree inference
#  --subsample NUM      Randomly sub-sample partitions (negative for complement)
#  --subsample-seed NUM Random number seed for --subsample

## LIKELIHOOD/QUARTET MAPPING:
#  --lmap NUM           Number of quartets for likelihood mapping analysis
#  --lmclust FILE       NEXUS file containing clusters for likelihood mapping
#  --quartetlh          Print quartet log-likelihoods to .quartetlh file


# MODEL-FINDER:
#  -m TESTONLY          Standard model selection (like jModelTest, ProtTest)
#  -m TEST              Standard model selection followed by tree inference
#  -m MF                Extended model selection with FreeRate heterogeneity
#  -m MFP               Extended model selection followed by tree inference
#  -m ...+LM            Additionally test Lie Markov models
#  -m ...+LMRY          Additionally test Lie Markov models with RY symmetry
#  -m ...+LMWS          Additionally test Lie Markov models with WS symmetry
#  -m ...+LMMK          Additionally test Lie Markov models with MK symmetry
#  -m ...+LMSS          Additionally test strand-symmetric models
#  --mset STRING        Restrict search to models supported by other programs
#                       (raxml, phyml or mrbayes)
#  --mset STR,...       Comma-separated model list (e.g. -mset WAG,LG,JTT)
#  --msub STRING        Amino-acid model source
#                       (nuclear, mitochondrial, chloroplast or viral)
#  --mfreq STR,...      List of state frequencies
#  --mrate STR,...      List of rate heterogeneity among sites
#                       (e.g. -mrate E,I,G,I+G,R is used for -m MF)
#  --cmin NUM           Min categories for FreeRate model [+R] (default: 2)
#  --cmax NUM           Max categories for FreeRate model [+R] (default: 10)
#  --merit AIC|AICc|BIC  Akaike|Bayesian information criterion (default: BIC)
#  --mtree              Perform full tree search for every model
#  --mredo              Ignore .model.gz checkpoint file (default: OFF)
#  --madd STR,...       List of mixture models to consider
#  --mdef FILE          Model definition NEXUS file (see Manual)
#  --modelomatic        Find best codon/protein/DNA models (Whelan et al. 2015)

# The special MFP key word stands for ModelFinder Plus, which tells IQ-TREE to perform ModelFinder and the remaining analysis using the selected model. ModelFinder computes the log-likelihoods of an initial parsimony tree for many different models and the Akaike information criterion (AIC), corrected Akaike information criterion (AICc), and the Bayesian information criterion (BIC). Then ModelFinder chooses the model that minimizes the BIC score (you can also change to AIC or AICc by adding the option -AIC or -AICc, respectively).

#SUBSTITUTION MODEL:
#  -m STRING            Model name string (e.g. GTR+F+I+G)
#                 DNA:  HKY (default), JC, F81, K2P, K3P, K81uf, TN/TrN, TNef,
#                       TIM, TIMef, TVM, TVMef, SYM, GTR, or 6-digit model
#                       specification (e.g., 010010 = HKY)
#             Protein:  LG (default), Poisson, cpREV, mtREV, Dayhoff, mtMAM,
#                       JTT, WAG, mtART, mtZOA, VT, rtREV, DCMut, PMB, HIVb,
#                       HIVw, JTTDCMut, FLU, Blosum62, GTR20, mtMet, mtVer, mtInv, Q.LG
#			Q.pfam, Q.pfam_gb, Q.bird, Q.mammal, Q.insect, Q.plant, Q.yeast
#     Protein mixture:  C10,...,C60, EX2, EX3, EHO, UL2, UL3, EX_EHO, LG4M, LG4X
#              Binary:  JC2 (default), GTR2
#     Empirical codon:  KOSI07, SCHN05
#   Mechanistic codon:  GY (default), MG, MGK, GY0K, GY1KTS, GY1KTV, GY2K,
#                       MG1KTS, MG1KTV, MG2K
#Semi-empirical codon:  XX_YY where XX is empirical and YY is mechanistic model
#      Morphology/SNP:  MK (default), ORDERED, GTR
#      Lie Markov DNA:  1.1, 2.2b, 3.3a, 3.3b, 3.3c, 3.4, 4.4a, 4.4b, 4.5a,
#                       4.5b, 5.6a, 5.6b, 5.7a, 5.7b, 5.7c, 5.11a, 5.11b, 5.11c,
#                       5.16, 6.6, 6.7a, 6.7b, 6.8a, 6.8b, 6.17a, 6.17b, 8.8,
#                       8.10a, 8.10b, 8.16, 8.17, 8.18, 9.20a, 9.20b, 10.12,
#                       10.34, 12.12 (optionally prefixed by RY, WS or MK)
#      Non-reversible:  STRSYM (strand symmetric model, equiv. WS6.6),
#                       NONREV, UNREST (unrestricted model, equiv. 12.12)
#           Otherwise:  Name of file containing user-model parameters

#STATE FREQUENCY:
#  -m ...+F             Empirically counted frequencies from alignment
#  -m ...+FO            Optimized frequencies by maximum-likelihood
#  -m ...+FQ            Equal frequencies
#  -m ...+FRY           For DNA, freq(A+G)=1/2=freq(C+T)
#  -m ...+FWS           For DNA, freq(A+T)=1/2=freq(C+G)
#  -m ...+FMK           For DNA, freq(A+C)=1/2=freq(G+T)
#  -m ...+Fabcd         4-digit constraint on ACGT frequency
#                       (e.g. +F1221 means f_A=f_T, f_C=f_G)
#  -m ...+FU            Amino-acid frequencies given protein matrix
#  -m ...+F1x4          Equal NT frequencies over three codon positions
#  -m ...+F3x4          Unequal NT frequencies over three codon positions

#RATE HETEROGENEITY AMONG SITES:
#  -m ...+I             A proportion of invariable sites
#  -m ...+G[n]          Discrete Gamma model with n categories (default n=4)
#  -m ...*G[n]          Discrete Gamma model with unlinked model parameters
#  -m ...+I+G[n]        Invariable sites plus Gamma model with n categories
#  -m ...+R[n]          FreeRate model with n categories (default n=4)
#  -m ...*R[n]          FreeRate model with unlinked model parameters
#  -m ...+I+R[n]        Invariable sites plus FreeRate model with n categories
#  -m ...+Hn            Heterotachy model with n classes
#  -m ...*Hn            Heterotachy model with n classes and unlinked parameters
#  --alpha-min NUM      Min Gamma shape parameter for site rates (default: 0.02)
#  --gamma-median       Median approximation for +G site rates (default: mean)
#  --rate               Write empirical Bayesian site rates to .rate file
#  --mlrate             Write maximum likelihood site rates to .mlrate file



