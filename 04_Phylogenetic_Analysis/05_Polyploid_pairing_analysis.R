
## Luis Leal 2021


## Script used to identify the sister species of the polyploid
## Input file is the ASTRAL input file containing a list of gene trees


##### Clear all states
if(!is.null(dev.list())) dev.off()
par(mfrow=c(1,1))
rm(list=ls(all=TRUE))


#library(BiocManager)
library("ape")
library("ips")

## Reference APE package
## Paradis E. et al.  (2004 ) APE: analyses of phylogenetics and evolution in R language . Bioinformatics , 20 , 289 –290.



#####Verify versions of R and other loaded software:
#sessionInfo()



##### read arguments
args <- commandArgs(TRUE)

# Import pre-ASTRAL file containing gene trees estimated using IQTREE2
gene_Trees_file <- args[1]
gene_Trees <- read.table(gene_Trees_file, header=FALSE, sep="", dec=".")

# Polyploid
POLY_ID <- args[2]

# Outgroup
outgroup_ID <- args[3]

# output folder
OUTFolder <- args[4]

# output file
OUTFile <- args[5]



##### for each  gene tree, identify sister species of polyploid

sister_ID <- c()

for (i in c(1:nrow(gene_Trees))){

   mytr <- read.tree(text = gene_Trees[i,1])
   mytr <- root(mytr , outgroup_ID)

   #### Find sister species of polyploid
   sister_ANC1 <- sister(mytr, POLY_ID, type = "terminal", label = TRUE)
   sister_ANC1 <- paste(sort(sister_ANC1), collapse = '_')
   #print(sister_ANC1)

   sister_ID <- append(sister_ID, sister_ANC1 )

}



##### frequency for sister species; order table, most frequent on top

counter_sister_df <- as.data.frame(table(sister_ID))
counter_sister_df <- counter_sister_df[order(-counter_sister_df$Freq), ]


##### normalized frequency

counter_sister_df$Freq_norm <- counter_sister_df$Freq/sum(counter_sister_df$Freq)



############ Save output to file

setwd(OUTFolder)

write.table(counter_sister_df, 
            file = OUTFile, sep = "\t", quote = FALSE, row.names=FALSE)


