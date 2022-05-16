#  
# Luis Leal (2020)
#
# Written in Python 3
# UPPMAX: module load python3/3.6.0



### Script used to fix invariant sites in VCF file
### i) if invariant site has DP>0 in GT filed, then PASS and 0/0
### ii) if invariant site has DP=0 in GT field, then FAIL and ./.




print('\n Parsing started ...')




######################################################### LOAD STANDARD MODULES

import sys
import re                                                 # 're' module provides regular expression matching operations
import ast                                                # required to convert reference tree from string to list
import os                                                 # required to create new folders
import time
from pdb import set_trace as bp                           # when using break points during debugging >> useage: bp()






######################################################## OPEN INPUT FILES


try:
    fhand = open(sys.argv[1], 'r')                  			# vcf file name
except:
    print('\n Error: input file missing.')
    exit() 




######################################################## Read input file, save parameters of interest to file


### AUXILIARY FUNCTION: searches 'Fixed fields' string (VCF) and gets numeric value of a specific parameter

def getValue(rd_aux00, queryString):
    rm_mark0 = rd_aux00.find(queryString)
    if rm_mark0 != -1 :                                             #if parameter present in x[7]
        rd_aux0 = rd_aux00[(rm_mark0+len(queryString)):]
        rm_mark = rd_aux0.find(';')
        if rm_mark != -1 :                                          #when parameter is in the middle of x[7]
            try :
                outvalue = float(rd_aux0[:rm_mark])       		    # get parameter value
            except: outvalue = ''      								# parameter value is a string (e.g., 'NaN')
        else :
            P_aux1 = rd_aux00[(rm_mark0+len(queryString)):]         #when parameter is at the end of x[7]
            try :
                outvalue = float(P_aux1)       		                # get parameter value
            except: outvalue = ''      								# parameter value is a string (e.g., 'NaN')
    else :
        outvalue = ''
    return outvalue


### end function





#open output file 
 
#output files
straux1=str(sys.argv[1])
straux1 = straux1[:(len(straux1)-4)]

#print('straux1',straux1)
 
outfile1_name = straux1 + '-CLEAN.vcf'
outfile1 = open(outfile1_name, 'w')         




                



#flags and counters
counterALL = 0          # counts all sites (variant and non-variant)
ResultsList = list()    # stores output results
FLAG_PLOIDY = 1             # indicates whether ploidy has been detected




# check each individual variant/invariant site
for line in fhand:
    #
    # detect comment lines
    xx = re.findall('^#', line)
    #
    if len(xx) > 0 : 
        # comment lines
        ResultsList.append(line)      # store comment lines
        #
    else :
        # variants
        # counter, all sites 
        counterALL = counterALL + 1
        #
        x = line.split()
        #
        # detect genotype
        if FLAG_PLOIDY == 1 :
            PLY_aux1 = x[9]             # get genotype fields
            PLY_aux2 = PLY_aux1.count('/') 
            if PLY_aux2 == 1 :      # diploid (eg 0/0)
                GT_PASS = '0/0'
                GT_FAIL = './.'
            elif PLY_aux2 == 3 :   # tetraploid (eg 0/0/0/0)
                GT_PASS = '0/0/0/0'
                GT_FAIL = './././.'
            elif PLY_aux2 == 5 :   # hexaploid (eg 0/0/0/0/0/0)
                GT_PASS = '0/0/0/0/0/0'
                GT_FAIL = './././././.'
            elif PLY_aux2 == 7 :   # octaploid (eg 0/0/0/0/0/0/0/0)
                GT_PASS = '0/0/0/0/0/0/0/0'
                GT_FAIL = './././././././.'
            else :
                print('\n Error: ploidy not detected.')
                exit()
            #
            FLAG_PLOIDY = 0       
        #
        #print('GT_PASS=',GT_PASS)
        #print('GT_FAIL=',GT_FAIL)
        #bp()
        #
        ALT = x[4]
        #print(ALT)
        # 
        if (ALT != '.') :
            #
            # variant sites (do not change; store site info)
            ResultsList.append(line)      # store variant sites (it doesn't matter if they PASS or FAIL)
            #print(ResultsList[-1])
            #print()
            #print(x)
            #bp()
            #
        else :
            #
            # non-variant sites
            GTleg = x[8]                	        # GT fields legend
            GTleg_split = GTleg.split(':')          # split GT fields legend
            GTfield = x[9]                          # GT fields
            GT_split = GTfield.split(':')           # split GT field
            #
            GTindex = GTleg_split.index('GT')       # find 'GT' position
            #
            # find 'DP' position
            DP = -1
            try:
                DPindex = GTleg_split.index('DP')      # try genotype fields 
                DP = float(GT_split[DPindex])          # read depth
            except:
                DP = 0                              # if read depth not listed, set DP to zero
            #
            if DP == 0 :
                # site has not been genotyped 
                x[6] = 'Depth_filter_min'          # set filter to FAIL (use filter already in use so that VCF file header does need to be updated)
                GT_split[GTindex]=GT_FAIL           # update genotype
                x[9]=':'.join(GT_split)             # reformat genotype fields to VCF format
                line_aux = '\t'.join(x) + '\n'      # reverse modified x to string; add new line character
                ResultsList.append(line_aux)        # store modified variant
                #print(ResultsList[-1])
                #print()
                #print(x)
                #bp()
                #
            elif DP > 0 :
                #
                # site has been genotyped by at least one read
                x[6] = 'PASS'                       # set filter to PASS
                GT_split[GTindex]=GT_PASS           # update genotype
                x[9]=':'.join(GT_split)             # reformat genotype fields to VCF format
                line_aux = '\t'.join(x) + '\n'     # reverse modified x to string; add new line character
                ResultsList.append(line_aux)        # store modified variant
                #print(ResultsList[-1])
                #print()
                #print(x)
                #bp()
            else :
                print('\n Error: could not determine depth.')
                print(x)
                exit()
                
                        

# save results to file
for i in range(len(ResultsList)) :
    outfile1.write(ResultsList[i])
    #outfile1.write('\n')



# print  number of sites to slurm file
print()
print('TOTAL number of sites (variant + invariant):', counterALL)

print('\n Done!')

outfile1.close()
