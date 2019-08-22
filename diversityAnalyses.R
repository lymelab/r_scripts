####Diversity analyses (Weighted Unifrac, Aitchison Distance, Bray Curtis Dissimilarity)

###########
#LIBRARIES
###########
library(phylofactor)
library(ape)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(compositions)
library(ape)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(phangorn)
library(vegan)
library(factoextra)

PATH="/Users/mann/github/triatomes"

setwd(PATH)

#########################
#LOAD, FILTER, PREP DATA
#########################
#load sample data
rawmetadata <- read_delim(file = file.path(PATH, "map.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries

seqtab.filtered <- read.table("sequence_table.16s.filtered.txt", header=T, row.names=1)
seqtab.clr <- read.table("sequence_table.16s.clr.txt", header=T, row.names=1)
taxa <- read.table("assigntax/rep_set_fix_tax_assignments.txt", header=F, sep="\t", row.names=1)

#record samples absent in either metadata or OTU table
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$SampleID)
notinraw <- setdiff(rawmetadata$SampleID, row.names(seqtab.filtered))


#load representative tree
tree <- read.tree("rep_set.filt.tre")

#clean up rownames in filtered sequence table
# new.rownames <- row.names(seqtab.filtered)
# new.rownames <- gsub("_S.*.fastq.gz", "", new.rownames)
# new.rownames <- paste("Tri", new.rownames, sep="") # append characters to numerical samples names -- these can mess things up in the long run
# rownames(seqtab.filtered) <- new.rownames
# write.table(data.frame("row_names"=rownames(seqtab.filtered),seqtab.filtered),"sequence_table.16s.filtered.txt", row.names=F, quote=F, sep="\t")

#add character to sample names in metadata file
#set row names as the sample IDs for the metadata
rownames(rawmetadata) <- rawmetadata$SampleID
new.samplenames <- paste("Tri", rownames(rawmetadata), sep="")
rawmetadata$SampleID <- new.samplenames

rownames(rawmetadata) <- rawmetadata$SampleID
#create phyloseq object from "seqtab.filtered", "rawmetadata", "taxa", and "tree"
ps.dada2_join <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa[1])), tree)

# at this point I recommend saving your full unfiltered dataset in phyloseq format as a .RDS, so that you can pick up the analysis from this point easily if you decide to change your filtering criteria later on
saveRDS(ps.dada2_join, paste(PATH,"filtered_dataset.phyloseq_format.RDS",sep=""))
 
# OPTIONAL: Filter out samples or OTUs with low read counts                                     
# # Remove samples with fewer than N reads (N=100 in example. adjust per experiment.)
# ps.dada2_join <- prune_samples(sample_sums(ps.dada2_join) >= 100, ps.dada2_join)
# #Remove OTUs with fewer than N total reads. (N = 50 for example. adjust per experiment)
# ps.dada2_join <- prune_taxa(taxa_sums(ps.dada2_join) >= 50, ps.dada2_join)


#finally, generate Centered Log Ratio Transformed (CLR) dataset for compositional data
comp <- clr(seqtab.filtered, comp.filt)
write.table(as.data.frame(comp), "sequence_table.16s.clr.txt", row.names=F, quote=F, sep="\t")

#####################
#AITCHISON DISTANCE
#####################
# cdata <- acomp(seqtab.filtered) # convert frequency data to an Aitchison compositional dataset -- zeros are flagged as BDL
hc <- hclust(dist(seqtab.clr), method="complete")
df2 <- data.frame(cluster=cutree(hc,6), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$cluster))) + geom_tile() + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("aitchison.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

####Aitchison PCoA ordination
# write.table(cdata, "aitchison.txt", quote=F, sep="\t")
# aitch <- read.table("aitchison.txt", header=T, row.names=1)
pca <- prcomp(as.matrix(seqtab.clr))
pdf("aitchison_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("aitchison_pca.pdf")
fviz_pca_biplot(pca, habillage=merge$cluster, geom="point") + labs(title="Aitchison Distance") 
dev.off()

###########################
#BRAY CURTIS DISSIMILARITY
###########################
seqtab.filtered.rare <- rarefy(seqtab.filtered, 15000)
bcdist <- vegdist(seqtab.filtered.rare, method="bray")
hc <- hclust(bcdist, method="complete")
df2 <- data.frame(cluster=cutree(hc,10), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$cluster))) + geom_tile() + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("figs/braycurtis.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()
#pca plot of rarefied data
pca <- prcomp(bcdist)
pdf("figs/braycurtis_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/braycurtis_pca_tcruzi.pdf")
fviz_pca_ind(pca, habillage=merge$MachRes1, geom="point") + labs(title="Bray Curtis Distance")
dev.off()
pdf("figs/braycurtis_pca_species.pdf")
fviz_pca_ind(pca, habillage=merge$Sp_by_key, geom="point") + labs(title="Bray Curtis Distance")
dev.off()


##########################
#WEIGHTED UNIFRAC DISTANCE
##########################
ps.dada2_join.rare <- rarefy_even_depth(ps.dada2_join, sample.size=15000, rngseed=456)
wuni <- UniFrac(ps.dada2_join.rare, weighted=T)
hc <- hclust(wuni, method="complete")
df2 <- data.frame(cluster=cutree(hc,10), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,20) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$cluster))) + geom_tile() + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("wunifrac.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

#pca plot of rarefied data
pca <- prcomp(wuni)
pdf("figs/wunifrac_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/wunifrac_pca_tcruzi.pdf")
fviz_pca_ind(pca, habillage=merge$MachRes1, geom="point") + labs(title="Weighted Unifrac Distance")
dev.off()
pdf("figs/wunifrac_pca_species.pdf")
fviz_pca_ind(pca, habillage=merge$Sp_by_key, geom="point") + labs(title="Weighted Unifrac Distance")
dev.off()

########################
#RAREFACTION CURVE PLOT
########################
pdf("figs/rarefaction_curve.pdf")
rarecurve(seqtab.filtered, label=FALSE)
dev.off()