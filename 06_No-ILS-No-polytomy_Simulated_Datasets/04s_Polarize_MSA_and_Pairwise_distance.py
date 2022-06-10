#  
# Luis Leal (2021)
#
# Written in Python 3



## Script used to process simulated MSAs:
## (1) creates diploid species (collapses sequences for the two haploid individuals belonging to the same species);
## (2) creates allotetraploid sequence (collapses sequences for two diploid species);
## (3) polarizes allotetraploid sequence.

### Also computes pairwise distance between each sequence in MSA and the polarized polyploid.
###

### Input MSA files must contain two haploid sequences per species (only one sequence for outgroup [0_0_0])
### It is assumed that input file contains NO indels (i.e. sequences are aligned)


### Input arguments
### 1. MSA for each locus simulated using INDELible (fasta format)
### 2. Reference sequence ID
### 3. ID of tetraploid ancestral species


#print('\n Parsing started ...')



######################################################### LOAD STANDARD MODULES

import sys
import re                                                 # 're' module provides regular expression matching operations
import ast                                                # required to convert reference tree from string to list
import os                                                 # required to create new folders
import time
from collections import Counter
from pdb import set_trace as bp                           # when using break points during debugging >> useage: bp()
import random						  # to select random element from array
import math







######################################################## OPEN INPUT FILES


try:
    fhand = open(sys.argv[1], 'r')                  			# input fasta file
except:
    print('\n Error: input file missing.')
    exit() 



try:
    RefGenID = sys.argv[2]                  				# Reference sequence ID
except:
    print('\n Error: provide the following arguments:')
    print('03s_polarizeMSA_fullCODED_ALT.py [input.fasta] [refSequence ID] [ancestral ID1] [ancestral ID2] \n')
    exit() 


try:
    AncestralID_1 = sys.argv[3]                  				# ID first ancestral sequence  
except:
    print('\n Error: provide the following arguments:')
    print('03s_polarizeMSA_fullCODED_ALT.py [input.fasta] [refSequence ID] [ancestral ID1] [ancestral ID2] \n')
    exit() 


try:
    AncestralID_2 = sys.argv[4]                  				# ID second ancestral sequence  
except:
    print('\n Error: provide the following arguments:')
    print('03s_polarizeMSA_fullCODED_ALT.py [input.fasta] [refSequence ID] [ancestral ID1] [ancestral ID2] \n')
    exit() 


















### AUXILIARY FUNCTION: code nucleotide; extended coding for heterozygotic sites according to IUPAC nomenclature:
### See also: http://www.iqtree.org/doc/Frequently-Asked-Questions#how-does-iq-tree-treat-gapmissingambiguous-characters

def extendedCODE(input_list):
   species_nuc_uniq = list(set(input_list))  # get unique nucleotides
   species_nuc_uniq.sort()		     # sort list alphabetically

   outvalue = ''
   if species_nuc_uniq == ['A'] : outvalue = 'A'
   if species_nuc_uniq == ['C'] : outvalue = 'C'
   if species_nuc_uniq == ['G'] : outvalue = 'G'
   if species_nuc_uniq == ['T'] : outvalue = 'T'
   if species_nuc_uniq == ['A','C'] : outvalue = 'M'
   if species_nuc_uniq == ['A','G'] : outvalue = 'R'
   if species_nuc_uniq == ['A','T'] : outvalue = 'W'
   if species_nuc_uniq == ['C','G'] : outvalue = 'S'
   if species_nuc_uniq == ['C','T'] : outvalue = 'Y'
   if species_nuc_uniq == ['G','T'] : outvalue = 'K'
   if species_nuc_uniq == ['A','C','G'] : outvalue = 'V'
   if species_nuc_uniq == ['A','C','T'] : outvalue = 'H'
   if species_nuc_uniq == ['A','G','T'] : outvalue = 'D'
   if species_nuc_uniq == ['C','G','T'] : outvalue = 'B'
   if species_nuc_uniq == ['A','C','G','T'] : outvalue = 'N'

   return outvalue


### end function




### AUXILIARY FUNCTION: expand IUPAC code

def unfoldCODE(species_nuc):

   outvalue = ''
   if species_nuc == 'A' : outvalue = ['A']
   if species_nuc == 'C' : outvalue = ['C']
   if species_nuc == 'G' : outvalue = ['G']
   if species_nuc == 'T' : outvalue = ['T']
   if species_nuc == 'M' : outvalue = ['A','C']
   if species_nuc == 'R' : outvalue = ['A','G']
   if species_nuc == 'W' : outvalue = ['A','T']
   if species_nuc == 'S' : outvalue = ['C','G']
   if species_nuc == 'Y' : outvalue = ['C','T']
   if species_nuc == 'K' : outvalue = ['G','T']
   if species_nuc == 'V' : outvalue = ['A','C','G']
   if species_nuc == 'H' : outvalue = ['A','C','T']
   if species_nuc == 'D' : outvalue = ['A','G','T']
   if species_nuc == 'B' : outvalue = ['C','G','T']
   if species_nuc == 'N' : outvalue = ['N']
   if outvalue == '' :
      print('\n Error: sequences contain non-IUPAC nucleotide character\n')
      exit()

   return outvalue

### end function










######################################################## Read MSA

INfasta_species = list()
INfasta_seq = list()
FLAG_aux = 0
for line in fhand:
    aux_NC = line.replace("\n", "")
    if FLAG_aux == 0 :
      INfasta_species.append(aux_NC)
      FLAG_aux = 1
    else :
      INfasta_seq.append(aux_NC)
      FLAG_aux = 0


INfasta_species_CLEAN = INfasta_species
INfasta_seq_CLEAN = INfasta_seq


# save species/haplotype and associated sequence to dictionary
INfasta_dict = dict()
for i in range(len(INfasta_species_CLEAN)) :
   INfasta_dict[INfasta_species_CLEAN[i]] = INfasta_seq_CLEAN[i]

#for i in range(len(INfasta_species_CLEAN)) :
#   print(INfasta_species_CLEAN[i], INfasta_dict[INfasta_species_CLEAN[i]])

#print(len(INfasta_dict))



# get species list from fasta file
species_list = list()
for tag in INfasta_species_CLEAN :
   tag = tag.replace(">", "")                   #remove > from sequence name
   species_list.append(tag[:(len(tag) - 4)])	#remove "_0_0" from "14_0_0" ["14" is the species name]

species_list = list(set(species_list))		# unique elements only


# Check number of sequences
#DIC_LEN=len(INfasta_dict)
#if DIC_LEN != (2 * int(Nspecies) -1) :
#   print('\n Error: check number of sequences in fasta file\n')
#   exit()

#if len(species_list) != int(Nspecies) :
#   print('\n Error: check number of sequences in fasta file\n')
#   exit()




# Check all sequences have the same length
FLAG_aux = 0
for i in range(len(INfasta_species_CLEAN)) :
   len_seq = len(INfasta_dict[INfasta_species_CLEAN[i]])
   #print(len_seq)
   if FLAG_aux == 0 :
      len_ref = len_seq
      FLAG_aux = 1
   else :
      if len_seq != len_ref :
         print('\n Error: sequences have different lenghts\n')
         exit()

          
# Check outgroup is not one of the ancestral sequences
#if (AncestralID_1 == "0") or (AncestralID_2 == "0") :
#    print('\n Error: ancestral sequence cannot be the outgroup\n')
#    exit()










###################################################### Create diploid/tetraploid sequences



## Create tetraploid sequences, IUPAC coded
ANC1_ID = '>' + AncestralID_1 + '_0_0'
ANC2_ID = '>' + AncestralID_1 + '_0_1'
ANC3_ID = '>' + AncestralID_2 + '_0_0'
ANC4_ID = '>' + AncestralID_2 + '_0_1'

ANC1seq = list()
ANC2seq = list()
ANC3seq = list()
ANC4seq = list()

for i in INfasta_dict[ANC1_ID] :
   ANC1seq.append(i)

for i in INfasta_dict[ANC2_ID] :
   ANC2seq.append(i)

for i in INfasta_dict[ANC3_ID] :
   ANC3seq.append(i)

for i in INfasta_dict[ANC4_ID] :
   ANC4seq.append(i)

TETRAseq = list()
TETRAseq_aux = list(zip(ANC1seq,ANC2seq,ANC3seq,ANC4seq))


# reduce each locus to its variants, IUPAC coded (mimicking MSA obtained from WGS data)
for i in TETRAseq_aux :
   species_nuc_uniq = extendedCODE(i)
   TETRAseq.append(species_nuc_uniq)
   #print(i,species_nuc_uniq)

TETRAseq = "".join(TETRAseq)



## Create diploid sequences (IUPAC coded)

DIPdic = INfasta_dict	# dictionary only with diploid sequences + outgroup (before linking haplotypes for same individual)
#del DIPdic[ANC1_ID]
#del DIPdic[ANC2_ID]
#del DIPdic[ANC3_ID]
#del DIPdic[ANC4_ID]

DEP_IDs = list()
for key, value in DIPdic.items():
    tag = key.replace(">", "")                   #remove > from sequence name
    DEP_IDs.append(tag[:(len(tag) - 4)])	#remove "_0_0" from "14_0_0" ["14" is the species name]

DEP_IDs = list(set(DEP_IDs))		# unique elements only
DEP_IDs.remove('0')			# remove outgroup


DICT_seq = dict()	# dictionary containing linked haplotypes, all species

for item in DEP_IDs :
   hap_A = '>' + item + '_0_0'
   hap_B = '>' + item + '_0_1'
   seq_A = list()
   seq_B = list()
   
   for i in INfasta_dict[hap_A] :
      seq_A.append(i)

   for i in INfasta_dict[hap_B] :
      seq_B.append(i)

   hap_AB = list()
   hap_AB_aux = list(zip(seq_A,seq_B))
   
   for i in hap_AB_aux :
      species_nuc_uniq = extendedCODE(i)
      hap_AB.append(species_nuc_uniq)
      #print(item, i,species_nuc_uniq)

   hap_AB = "".join(hap_AB)

   DICT_seq[item] = hap_AB



 
## add outgroup to dictionary
hap_A = '>' + '0' + '_0_0'
DICT_seq['0'] = INfasta_dict[hap_A]


## add tetraploid to dictionary
POLYPID = 'T-' + AncestralID_1 + '-' + AncestralID_2
DICT_seq[POLYPID] = TETRAseq


## tag reference sequence
REFG_ID =  RefGenID + '_Ref'
DICT_seq[REFG_ID] = DICT_seq.pop(RefGenID)

## tag parental species
PAR1_TAG = AncestralID_1 + '_PAR'
PAR2_TAG = AncestralID_2 + '_PAR'
DICT_seq[PAR1_TAG] = DICT_seq.pop(AncestralID_1)
DICT_seq[PAR2_TAG] = DICT_seq.pop(AncestralID_2)


#for key, value in DICT_seq.items():
#   print(key,value)

#for key, value in DICT_seq.items():
#   print(key,len(DICT_seq[key]))







###################################################### Unfold tetraploid and ref sequences

## Reference sequence
REFseq = list()
for i in DICT_seq[REFG_ID] :
   nuc_ste = unfoldCODE(i)
   REFseq.append(nuc_ste)
   #print(i,nuc_ste)


## Polyploid sequence
TETRAseq_u = list()
for i in DICT_seq[POLYPID] :
   nuc_ste = unfoldCODE(i)
   TETRAseq_u.append(nuc_ste)
   #print(i,nuc_ste)







###################################################### Create polarized sequence for ALLOTETRAPLOID


DICT_POLAR_ALT = dict()	        # dictionary containing polarized polyploid sequence (also all other sequences, but which stay as in the input MSA)


###  Add non-polyploid sequences to dictionary
for key, value in DICT_seq.items():
   aux_ID = '>' + key
   DICT_POLAR_ALT[aux_ID] = value
   #print(key)
   #print(key,DICT_POLAR_ALT[aux_ID])


# remove unpolarized polyploid from dictionary (polarized polyploid will be added in the next step)
POLY_ID = '>' + POLYPID
in_POLY = DICT_POLAR_ALT[POLY_ID]
del DICT_POLAR_ALT[POLY_ID]

#for key, value in DICT_POLAR_ALT.items():
   #print(key,len(value))
   #print(key,value)



### Allotetraploid

ALT_list = list()
species_seq = TETRAseq_u[:]
for i in range(len(REFseq)) :		# for each nucleotide in reference sequence
   ref_nuc_ls = REFseq[i]
   species_nuc = list(species_seq[i])
   #print(ref_nuc_ls, species_nuc)
   #
   if species_nuc[0] == 'N' :		# if polyploid is masked
      aux_nuc = 'N' 
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc) 
   elif ref_nuc_ls[0] == 'N' :		# if reference is masked
      aux_nuc = 'N' 
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   elif (species_nuc == ['A'] or species_nuc == ['C'] or species_nuc == ['G'] or species_nuc == ['T']) :		# polyploid if fixed
      ALT_list.append(species_nuc[0])
      #print(ref_nuc_ls, species_nuc, species_nuc[0])
   elif species_nuc == ref_nuc_ls :			# if reference alleles were found in the polyploid and there are no other alleles present
      aux_nuc = extendedCODE(species_nuc)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   elif all(x0 in species_nuc for x0 in ref_nuc_ls) :		# if all reference alleles were found also in polyploid, but there are other alleles present as well
      aux_nuc0 = [pp for pp in species_nuc if pp not in ref_nuc_ls]		# remove reference alleles from polyploid allele set
      aux_nuc = extendedCODE(aux_nuc0)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   elif any(x0 in ref_nuc_ls for x0 in species_nuc) :		# if at least one of the reference alleles is found also in the polyploid
      aux_nuc0 = [pp for pp in species_nuc if pp not in ref_nuc_ls]		# remove reference alleles from polyploid allele set, if present
      aux_nuc = extendedCODE(aux_nuc0)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   else :							# reference allele not found in polyploid
      aux_nuc = extendedCODE(species_nuc)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)



# update dictionary
DICT_POLAR_ALT[POLY_ID]=''.join(ALT_list)

#for key, value in DICT_POLAR_ALT.items():
#   print(key)
#   print(key,len(DICT_POLAR_ALT[key]))






## Count differences between polarized and non-polarized polyploid sequence


tetra_RAW = in_POLY
tetra_ALT = DICT_POLAR_ALT[POLY_ID]

#counter_ALT_RAW = 0

#for i in range(len(tetra_RAW)) :		
#   hap_RAW = tetra_RAW[i]
#   hap_ALT = tetra_ALT[i]
#   if hap_RAW != hap_ALT :
#      counter_ALT_RAW += 1
#      #print(hap_RAW,hap_ALT)



## Compute pairwise distance between each sequence and the polarized polyploid [distance = SQRT(id/l) ; where id:the pairwise sequence identity; l:sequence length]
species_list_aux = list()
dist_by_species = list()
for key, value in DICT_POLAR_ALT.items() :
   counter_ALT = 0
   aux_spec = key
   aux_seq = value
   aux_len = len(value)
   for i in range(len(aux_seq)) :		
      hap_SPEC = aux_seq[i]
      hap_ALT = tetra_ALT[i]
      if hap_SPEC != hap_ALT :
         nuc_SPEC = unfoldCODE(hap_SPEC)
         nuc_ALT = unfoldCODE(hap_ALT)
         #print(nuc_SPEC, nuc_ALT)
         if (set(nuc_SPEC).issubset(nuc_ALT) or set(nuc_ALT).issubset(nuc_SPEC)) :			# accounts for multiple alleles
            counter_ALT += 0
         else :
            counter_ALT += 1
   #
   counter_ALT = math.sqrt(counter_ALT/aux_len)
   species_list_aux.append(key)
   dist_by_species.append(counter_ALT)



##### Save polarized MSA to file

#open output file 
straux1=str(sys.argv[1])
straux1 = straux1[:(len(straux1)-3)]
outfile_ALT = straux1 + '-ALT.fasta'
outfile_stats1 = straux1 + '.stats'

outfile2 = open(outfile_ALT, 'w')   
for key, value in DICT_POLAR_ALT.items() :
    outfile2.write(key)
    outfile2.write('\n')
    outfile2.write(value)
    outfile2.write('\n')

outfile3 = open(outfile_stats1, 'w')   
outfile3.write('\n')
outfile3.write('Number of nucleotide differences between polyploid and all other sequences:')
outfile3.write('\n')
for i in range(len(species_list_aux)) :
   out_aux = str(dist_by_species[i]) + " " + species_list_aux[i]
   outfile3.write(out_aux)
   outfile3.write('\n')
   #print(out_aux)

outfile2.close()
outfile3.close()



 
      

