##############################
#RANDOM FOREST CLASSIFICATION 
##############################

###########
#LIBRARIES
###########
library(plyr)
library(randomForest)
library(rfUtilities)

#####################
#CHANGE THESE VALUES
#####################
OTUTABLE = "path/to/my/otu_table.txt" #should be tab separated, rows are samples, columns are ASVs
METADAT = "path/to/my/metadata.txt" #metadata file, rows are samples, columns are variables
CAT = "category" #what variable are you comparing, this should be the column header name of your test group
VAR1 = "variable_1" #first variable group, for example if you want to compare males vs females this could be males
VAR2 = "variable_2" #second variable group, for example females
SUBSAMP= "10" #how to subsample your groups. random forest models are sensitive to unequal sampling sizes so you should pick an appropriate subset number

otu_table <- read.table(OTUTABLE, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
otu_table <- t(otu_table)
metadata <- read.table(METADAT, sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

otu_nonzero_counts <- apply(otu_table, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_rare_removed <- remove_rare(table=otu_table, cutoff_pro=0.1)
otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)  

# can model predict tortoise group? random subset of 31 each (all solar, random 31 necropsy)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))  
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), CAT]
otu_table_var1 <- subset(otu_table_scaled_var, var==VAR1)
otu_table_var2 <- subset(otu_table_scaled_var, var==VAR2)
randomRows = function(df,n){
   return(df[sample(nrow(df),n),])
}

otu_table_var1sub <- randomRows(otu_table_var1, SUBSAMP)
otu_table_var2sub <- randomRows(otu_table_var2, SUBSAMP)
#bind together
subset_otu_table <- rbind(otu_table_var1sub, otu_table_var2sub)
#reset factors
subset_otu_table$var <- factor(subset_otu_table$var)
# run random forest
set.seed(151)  
rf <- randomForest(x=subset_otu_table[,1:(ncol(subset_otu_table)-1)] , y=subset_otu_table[ , ncol(subset_otu_table)] , ntree=10000, importance=TRUE, proximities=TRUE )  
rf
# calculate significance of rf model
rf.significance(x=rf, xdata=subset_otu_table[,1:(ncol(subset_otu_table)-1)], nperm=1000, ntree=1000)
