#### This script performs a log-ratio transformation on compositional data
# Generates an Aitchison distance dendrogram and biplot
# Generates a PhIRL (phylogenetic) distance dendrogram and PCA plot - I haven't figured out how to do a biplot but if you do please add to the script!
# NOTE: your representative sequence tree MUST be rooted - you can do this with ape or in a tree viewer program (e.g., figtree)

###########
#LIBRARIES
###########
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(compositions)
library(ape)
library(factoextra)
library(RColorBrewer)
library(philr)
library(phyloseq)

#################
#VARIABLE SET UP
#################
PATH="/Users/mann/github/triatomes" # working directory
MAP="/Users/mann/github/triatomes/map.txt" # metadata file
SEQTAB="/Users/mann/github/triatomes/sequence_table.16s.filtered.txt" # sequence table 
TAX="/Users/mann/github/triatomes/assigntax/rep_set_fix_tax_assignments.txt" # taxonomy table
TREE="/Users/mann/github/triatomes/rep_set.filt.tre" # representative sequence tree

##########
#ENV SETUP
##########
setwd(PATH)
cols <- brewer.pal(6, "Set2") # ColorBrewer palette

###########
#LOAD DATA
###########
#load sample data
rawmetadata <- read_delim(file = file.path(PATH, MAP), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries
seqtab.filtered <- read.table(SEQTAB, header=T, row.names=1)
taxa <- read.table(TAX, header=F, sep="\t", row.names=1)
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$SampleID) #record samples absent in either metadata or OTU table
notinraw <- setdiff(rawmetadata$SampleID, row.names(seqtab.filtered))
tree <- read.tree(TREE) #load representative tree

####OPTIONAL: clean up rownames in your sequence table -- e.g., if your sample names are fastq file names
# new.rownames <- row.names(seqtab.filtered)
# new.rownames <- gsub("_S.*.fastq.gz", "", new.rownames)
# new.rownames <- paste("Tri", new.rownames, sep="") # append characters to numerical samples names -- these can mess things up in the long run
# rownames(seqtab.filtered) <- new.rownames
# write.table(data.frame("row_names"=rownames(seqtab.filtered),seqtab.filtered),"sequence_table.16s.filtered.txt", row.names=F, quote=F, sep="\t")
# rownames(rawmetadata) <- rawmetadata$SampleID 
# new.samplenames <- paste("Tri", rownames(rawmetadata), sep="") #add character to sample names in metadata file
# rawmetadata$SampleID <- new.samplenames
# rownames(rawmetadata) <- rawmetadata$SampleID #set row names as the sample IDs for the metadata

################
#NORMALIZE DATA
################
#Centered log-ratio transformation
seqtab.clr <- clr(seqtab.filtered)
write.table(as.data.frame(seqtab.clr), "sequence_table.16s.clr.txt", quote=F, sep="\t", col.names=NA)

#############################
#AITCHISON DISTANCE (LINEAR)
#############################
# Heirarchical cluster dendrogram
hc <- hclust(dist(seqtab.clr), method="complete")
df2 <- data.frame(cluster=cutree(hc,6), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$cluster))) + geom_tile() + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("figs/aitchison_dendrogram.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

# Biplot
#have to read back in clr sequence table
comp <- read.table("sequence_table.16s.clr.txt", header=T, row.names=1)
pca <- prcomp(as.matrix(comp))
pdf("aitchison_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("aitchison_pca.pdf")
fviz_pca_biplot(pca, habillage=merge$cluster, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F) + labs(title="Aitchison Distance", ylab="PC1 (9.3%)", xlab="PC2 (11.3%)") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-12.5,11.5) + ylim(-12.5,11.5)
dev.off()

##############################
#PHILR DISTANCE (PHYLOGENETIC)
##############################
#create phyloseq object from "seqtab.filtered", "rawmetadata", "taxa", and "tree"
ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa[1])), tree)

philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
is.rooted(phy_tree(philr.dat)) #check that tree is rooted
is.binary.tree(phy_tree(philr.dat)) #check that multichotomies are resolved in tree
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
otu.table <- otu_table(philr.dat)
tree <- phy_tree(philr.dat)
metadata <- sample_data(philr.dat)
tax <- tax_table(philr.dat)
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")

# Heirarchical cluster dendrogram
hc <- hclust(dist(philr.t), method="complete")
df2 <- data.frame(cluster=cutree(hc,6), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
write.table(df2, "philr_cluster.txt", quote=F, sep="\t", col.names=NA)

hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$cluster))) + geom_tile() + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("figs/philr_dendrogram.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

# PCA
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))
pdf("philr_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("philr_pca.pdf")
fviz_pca_ind(pca, habillage=merge$cluster, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F) + labs(title="PhIRL Distance", xlab="PC1 (44.7%)", ylab="PC2 (38.8%)") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-50,50) + ylim(-50,50)
dev.off()










