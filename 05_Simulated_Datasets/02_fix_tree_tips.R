
## Fix tree tips
## Luis Leal 2021



## Script used to fix tree tips so that the two (haploid) individuals belonging to the same species are monophyletic  
## (as these will be recoded as one single diploid individual at a later stage)


##### Clear all states
if(!is.null(dev.list())) dev.off()
par(mfrow=c(1,1))
rm(list=ls(all=TRUE))



#library(BiocManager)
#BiocManager::install('ape')
#BiocManager::install('ips')
#BiocManager::install('phytools')
#BiocManager::install('castor')
library("ape")
library("ips")
#library("phytools")
#library(castor)

## Reference APE package
## Paradis E. et al.  (2004 ) APE: analyses of phylogenetics and evolution in R language . Bioinformatics , 20 , 289 –290.


#Verify versions of R and other loaded software:
#sessionInfo()



## read arguments
args <- commandArgs(TRUE)

## input phylogeny
Reference_Tree <- args[1]
mytr_ref <- read.tree(text = Reference_Tree)
#print(mytr_ref)
#plot(mytr_ref)


# output folder
OUTFolder <- args[2]


# output file
OUTFile <- args[3]


# number of species
NSPEC <- args[4]
#print(NSPEC)



# For each species in tree, check whether both individuals form a monophyletic group. If not, remove one of the individuals.

id_kept <- c()

for (i in c(1:NSPEC)){			

   # leaf name of the two individuals associated to each species
   IND1 <- paste(i,"_0_0",sep='')
   IND2 <- paste(i,"_0_1",sep='')
   #print(IND1)

   # determine sister group for each individual
   sister_ANC1 <- sister(mytr_ref, IND1, type = "terminal", label = TRUE)
   sister_ANC2 <- sister(mytr_ref, IND2, type = "terminal", label = TRUE)
   
   # if the two individuals do nor form a monophyletic group, remove one of the individuals; save id of individual kept
   if (( sister_ANC1 != IND2) || ( sister_ANC2 != IND1)){
     # remove "X_0_1" individual
      mytr_ref <- drop.tip(mytr_ref, IND2, trim.internal = TRUE, subtree = FALSE,root.edge = 0)
     # Save ID individual kept
      id_kept <- append(id_kept,IND1)
   } 

}


## Set mimimum branch length to 0.01
mytr_ref$edge.length[mytr_ref$edge.length < 0.01] <- 0.01


############ Save output to file

setwd(OUTFolder)
write.tree(mytr_ref, file = OUTFile, append = FALSE)

if (length(id_kept) > 0){
   write(id_kept, file = "id_kept.txt", append = FALSE)
}

