################################
# Install packages (only once) #
################################
#install.packages('seqminer')
#install.packages("toprdata")
#install.packages('ggplot2')
#install.packages('stringr')
#install.packages('scales')
#install.packages("cutpointr")
#install.packages('plyr')
#install.packages('dplyr')


#################
# Load packages #
#################
library(seqminer)
library(toprdata)
library(ggplot2)
library(stringr)
library(scales)
library(cutpointr)
library(plyr)
library(dplyr)


#######################################
# Set your working dir and file paths #
#######################################
setwd("/Users/joeri/git/alphamissense-vkgl-threshold/outputs/")
# AlphaMissense, download from https://zenodo.org/records/8360242
alphaMissenseLoc <- "/Applications/AlphaFold2/AlphaMissense_hg19.tsv.gz"
# VKGL variant classifications, download from https://vkgl.molgeniscloud.org
vkglLoc <- "/Users/joeri/VKGL/VKGL-releases/VKGL_public_consensus_oct2023.tsv"


#######################################################
# Retrieve VKGL classifications and add AlphaMissense #
#######################################################
vkgl <- read.table(file=vkglLoc, sep = '\t',header = TRUE)
uniqGenes <- unique(vkgl$gene)
results <- data.frame()
for(i in 1:length(uniqGenes)){
  geneName <- uniqGenes[i]
  cat(paste("Working on ",geneName," (",i," of ",length(uniqGenes),")\n", sep=""))
  geneCoords <- subset(ENSGENES_37, gene_symbol==geneName)
  if(nrow(geneCoords)==0){
    next
  }
  geneChr <- gsub("chr","", geneCoords$chrom)
  geneTabix <- paste(geneCoords$chrom, paste(geneCoords$gene_start, geneCoords$gene_end, sep="-"), sep=":")
  geneAlphaMissenseData <- tabix.read(alphaMissenseLoc, geneTabix)
  alphaMissense <- read.table(text=geneAlphaMissenseData, sep = '\t', header = FALSE, fileEncoding = "UTF-16LE", col.names = c("CHROM", "POS", "REF", "ALT", "genome", "uniprot_id", "transcript_id", "protein_variant", "am_pathogenicity", "am_class"))
  vkglGeneWithAlph <- merge(alphaMissense, vkgl, by.x = c("POS","REF","ALT"), by.y = c("start","ref","alt"))
  vkglGeneWithAlph <- vkglGeneWithAlph[,c("am_pathogenicity","am_class","classification")]
  vkglGeneWithAlph <- subset(vkglGeneWithAlph, classification != "VUS")
  cat(paste("Adding ",nrow(vkglGeneWithAlph)," variants\n", sep=""))
  results <- rbind(results, vkglGeneWithAlph)
}
# Write results for next time
write.table(results, sep="\t",file="alphamissense-vkgl-results.txt", quote=FALSE, row.names =FALSE)


####################################################
# Determine optimal threshold using Youden's Index #
####################################################
# Read existing results
results <- read.table(file="alphamissense-vkgl-results.txt", sep = '\t',header = TRUE)
opt_cut <- cutpointr(results, am_pathogenicity, classification, direction = ">=", pos_class = "LP", neg_class = "LB", method = maximize_metric, metric = youden)
youdenIndex <- opt_cut$optimal_cutpoint
tp <- sum(results[results$classification=="LP",'am_pathogenicity'] >= youdenIndex)
fp <- sum(results[results$classification=="LB",'am_pathogenicity'] >= youdenIndex)
tn <- sum(results[results$classification=="LB",'am_pathogenicity'] < youdenIndex)
fn <- sum(results[results$classification=="LP",'am_pathogenicity'] < youdenIndex)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
sens <- opt_cut$sensitivity*100
spec <- opt_cut$specificity*100
cat(paste("The optimal AlphaMissense threshold based on ",dim(results)[1]," VKGL variant classifications is ",round(youdenIndex,5)," with PPV ",round(ppv),"%, NPV ",round(npv),"%, sensitivity ",round(sens),"% and specificity ",round(spec),"%.\n",sep=""))
# same but for AlphaMissense's own classification labels
resultsA <- subset(results, am_class == "benign" | am_class == "pathogenic")
opt_cut <- cutpointr(resultsA, am_pathogenicity, am_class, direction = ">=", pos_class = "pathogenic", neg_class = "benign", method = maximize_metric, metric = youden)
youdenIndexA <- opt_cut$optimal_cutpoint
tp <- sum(resultsA[resultsA$am_class=="pathogenic",'am_pathogenicity'] >= youdenIndexA)
fp <- sum(resultsA[resultsA$am_class=="benign",'am_pathogenicity'] >= youdenIndexA)
tn <- sum(resultsA[resultsA$am_class=="benign",'am_pathogenicity'] < youdenIndexA)
fn <- sum(resultsA[resultsA$am_class=="pathogenic",'am_pathogenicity'] < youdenIndexA)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
sens <- opt_cut$sensitivity*100
spec <- opt_cut$specificity*100
cat(paste("The optimal AlphaMissense threshold based on AlphaMissense's own classifications is ",round(youdenIndexA,5)," with PPV ",round(ppv),"%, NPV ",round(npv),"%, sensitivity ",round(sens),"% and specificity ",round(spec),"%.\n",sep=""))


#################
# Other outputs #
#################
# Confusion matrix between AlphaMissense label and VKGL label
table(results$am_class, results$classification)
# Density plot of AlphaMissense predictions vs VKGL label
ggplot(results, aes(am_pathogenicity, colour = classification, fill = classification)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title.x=element_text(size=10)) +
  geom_density(alpha = 0.1, adjust = 0.1) +
  scale_fill_manual(values=c("green", "red")) +
  scale_color_manual(values=c("green", "red")) +
  geom_vline(xintercept = youdenIndex) +
  ggtitle(paste("Optimal AlphaMissense threshold based on\nVKGL variant classifications: ",round(youdenIndex,5),sep="")) +
  xlab("AlphaMissense pathogenicity prediction") +
  ylab("Density") +
  guides(fill=guide_legend(title="VKGL\nvariant\nclassifi-\ncation")) +
  guides(color=guide_legend(title="VKGL\nvariant\nclassifi-\ncation"))
ggsave(paste("alphamissense-vkgl-threshold-density.png",sep=""), width=5, height=5)
# Density plot of AlphaMissense predictions vs AlphaMissense label
ggplot(resultsA, aes(am_pathogenicity, colour = am_class, fill = am_class)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title.x=element_text(size=10)) +
  geom_density(alpha = 0.1, adjust = 0.1) +
  scale_fill_manual(values=c("green", "red")) +
  scale_color_manual(values=c("green", "red")) +
  geom_vline(xintercept = youdenIndexA) +
  ggtitle(paste("Optimal AlphaMissense threshold based on\nAlphaMissense classifications: ",round(youdenIndexA,5),sep="")) +
  xlab("AlphaMissense pathogenicity prediction") +
  ylab("Density") +
  guides(fill=guide_legend(title="Alpha\nMissense\nvariant\nclassifi-\ncation")) +
  guides(color=guide_legend(title="Alpha\nMissense\nvariant\nclassifi-\ncation"))
ggsave(paste("alphamissense-amclass-threshold-density.png",sep=""), width=5, height=5)

