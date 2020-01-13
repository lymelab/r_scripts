#######
#DESEQ
#######
#after loading in your phyloseq object (here ps.dat), covert it to a deseq object
ps.deseq <- phyloseq_to_deseq2(ps.dat, ~ MYVARIABLE) #CHANGE ME to your desired grouping variable
desq <- DESeq(ps.deseq, test="Wald", fitType="parametric")
res <- results(desq, cooksCutoff=F)
sigtab <- res[which(res$padj < 0.01),] # set your desired p value cutoff (if you have a lot of taxa remaining, increase)
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps.dat)[rownames(sigtab),], "matrix"))

x <- tapply(sigtab$log2FoldChange, sigtab$V3, function(x) max(x)) #CHANGE ME to the taxonomic level you want to color your graph by, here V3==Phylum
x <- sort(x, TRUE)
sigtab$V3 <- factor(as.character(sigtab$V3), levels=names(x)) #CHANGE ME to the taxonomic level you want to color your graph by, here V3==Phylum

x <- tapply(sigtab$log2FoldChange, sigtab$V6, function(x) max(x)) #CHANGE ME to the taxonomic level you want to collapse your ASVs by, here V6==Genus
x <- sort(x, TRUE)
sigtab$V6 <- factor(as.character(sigtab$V6), levels=names(x)) #CHANGE ME to the taxonomic level you want to collapse your ASVs by, here V6==Genus
ggplot(sigtab, aes(x=V6, y=log2FoldChange, color=V3)) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip() #CHANGE ME to the two taxonomic levels you chose






