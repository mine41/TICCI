rm(list=ls())
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
getwd()

library(cluster)
library(Seurat)
library(dplyr)
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
library(scran)
library(DT)

source("./R_tools/slice.R")
source("./R_tools/cci.R")


# specify the path to the data directory
data.dir <- paste(getwd(),"/data/", sep="")
data.name <- "HSMM"
# a string describing the analysis and will be attached to the names of all the files generated during the execution
context_str = paste("SLICE-", data.name, "-", format(Sys.time(), "%b%d_%H_%M_%S"), "-", sep="") 


# 1.loading the HSMM data set————————————————————————————————————————————————————————————————————————————————————————————————————————————————
load(paste(data.dir, "HSMM.Rda", sep=""))
dim(HSMM@assayData[["exprs"]])

# the expressions, cell info, and gene info are stored as an Biobase::ExpressionSet object
es <- HSMM

# remove outlier cells with low number of expressed genes
exp.threshold=1
pData(es)$num_expressed <- colSums(exprs(es) > exp.threshold)

cat("\tNumber of Expressed Genes:\n")
cat("\tMean: ", mean(pData(es)$num_expressed),"; Median: ", median(pData(es)$num_expressed), "\n", sep="")
cat("\tMin: ", min(pData(es)$num_expressed),"; Max: ", max(pData(es)$num_expressed), "\n", sep="")
outlier.threshold <- 0.75 * as.numeric(quantile(pData(es)$num_expressed, probs=(0.1)))
outlier.cells <- which(pData(es)$num_expressed < outlier.threshold)
if (length(outlier.cells) > 0) {
  cat("\tThe following cells are detected as outliers and removed: ", paste(rownames(pData(es))[outlier.cells], collapse=" ", sep=" "), sep="\n")
  es <- es[, which(pData(es)$num_expressed >= outlier.threshold)]
} else {
  cat("\tNo outlier detected!\n")
}
# remove ERCC genes and Ribosomal genes
ercc.genes <- grep("^ERCC-", rownames(fData(es)), value = TRUE)
rb.genes <- grep("^RPL|^RPS|^MRPL|^MRPS", rownames(fData(es)), value = TRUE)
es <- es[which(!(rownames(fData(es)) %in% c(rb.genes, ercc.genes))), ]


# 2.Measuring cell differentiation state using scEntropy——————————————————————————————————————————————————————————————————————————————————————

# use the scRNA-seq of HSMM cells (n=266) to create a SLICE object
# use the predicted cell state information from the original analysis as the original cell identity
# set projname to describe this analysis
sc <- construct(exprmatrix=as.data.frame(exprs(es)), 
                cellidentity=factor(pData(es)$Category, levels=c("Proliferating cells","Differentiating myoblasts", "Interstitial mesenchymal cells")),
                projname=context_str)
# loading the pre-computed human gene-gene Kappa similarity matrix
load(paste(data.dir, "hs_km.Rda", sep=""))

# bootstrap calculation of scEntropy. 100 boostrapped samples, 1000 genes in each sample
sc <- getEntropy(sc, km=km, calculation="bootstrap", 
                 B.num=100, exp.cutoff=1, B.size=1000, 
                 clustering.k=floor(sqrt(1000/2)),  
                 random.seed=201602)
# visualize and compare differentiation states measured by 
# bootstrap and deterministic calculations of scEntropy
plotEntropies(sc)

# 3.processing————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# Convert to count data
counts <- as(es@assayData[["exprs"]], "matrix")
# Filter out genes present in less than 5% of cells
mincell <- as.integer(dim(counts)[2] * 0.05)
min.cells <- mincell
filter.index <- rowSums(counts > 0) >= min.cells
counts <- counts[filter.index, ]
# Calculate normalization factor for each sample to reach a sequencing depth of 1e4 per cell
size_factors <- 1e4 * colSums(counts) / sum(counts)
# Normalize gene expression
norm_exprs <- counts / size_factors
# Log1p transform the data
norm_exprs <- log1p(norm_exprs)
dim(norm_exprs)
# Save the preprocessed data
getwd()
save.image(paste(getwd(),"/Rdata/HSMM/1processed.RData", sep=""))


load(paste(getwd(),"/Rdata/HSMM/1processed.RData", sep=""))
# 4.cluster———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# Expression matrix for CellChat my_data
my_data = as.data.frame(t(norm_exprs))
# Gap statistic to calculate the best k value
gaps <- clusGap(my_data, FUN = kmeans, K.max = 10, B = 50)
plot(gaps)
print(gaps, method = "firstmax")
best_k <- which.max(gaps$Tab[, "gap"])
# K-means clustering
kmeans_result <- kmeans(my_data, centers = best_k)
my_meta <- data.frame(row.names = rownames(my_data), Classification = as.factor(kmeans_result$cluster))

# Plot PCA visualization
pca_result <- prcomp(my_data, center = TRUE, scale. = TRUE)
pca_components <- pca_result$x[, 1:2]
pca_df <- data.frame(x = pca_components[, 1], y = pca_components[, 2], meta = factor(my_meta$Classification))
p <- ggplot(pca_df, aes(x = x, y = y, color = meta)) +
  geom_point(size = 5) +
  ggtitle("HSMM_kmeans_Classification") +
  scale_color_discrete(name = "kmeans_Classification") +  # Set the color mapping name
  theme_minimal()
p
# Save the plot as SVG format
ggsave("HSMM_kmeans_Classification.svg", plot = p, width = 8, height = 4)
save.image(paste(getwd(),"/Rdata/HSMM/2clustered.RData", sep=""))

load(paste(getwd(),"/Rdata/HSMM/2clustered.RData", sep=""))
# 5.cellchat————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# Process data.input / meta
data.input <- as.matrix(do.call(rbind, my_data))
colnames(data.input) <- rownames(my_data)
data.input = as(data.input, "dgCMatrix")

meta <- my_meta
unique(meta$Classification)# View the cell annotation information stored in meta, which will be used as the grouping basis later

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Classification")
# Create CellChat object, group.by specifies the object for communication, using annotations in meta as the grouping basis
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "Classification")
groupSize <- as.numeric(table(cellchat@idents)) # Number of cells in each cell type

# View CellChatDB
# Choose the appropriate species, options are CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.human
# View the composition of the database
showDatabaseCategory(CellChatDB)
# View specific information in the database
# CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
cellchat@DB <- CellChatDB # Load the database content into the cellchat object, setting the reference database for the next steps

# Preprocess expression data
cellchat <- subsetData(cellchat) # Extract expression data
cellchat <- identifyOverExpressedGenes(cellchat) # Identify highly expressed genes
cellchat <- identifyOverExpressedInteractions(cellchat) # Identify highly expressed pathways
cellchat <- projectData(cellchat, PPI.human) # Project to PPI, store the results in cellchat@LR$LRsig

# Infer ligand-receptor communication network
# Calculate cell communication probability at the ligand-receptor level
gc()
cellchat <- computeCommunProb(cellchat, raw.use = T) # Default calculation method is type = "truncatedMean"
# Remove cells with very few communications, default cutoff is 20%, i.e., genes with an expression proportion below 25% will be considered 0, trim = 0.1 can adjust the proportion threshold
cellchat <- filterCommunication(cellchat)

# Calculate cell communication at the pathway level, CellChat calculates pathway-level communication probability by summarizing all ligand-receptor interactions related to each signaling pathway
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat) # Calculate the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
# Plot, 1 row and 2 columns
par(mfrow = c(1,2), xpd=TRUE)
# Use circle plot to display the number or total interaction strength (weight) between any two cell groups. Colors match the source, and circle size represents the number of cells in each cell group
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")

# Total interaction strength
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights")
# Manually add legend
legend("bottomright",y=0.1, legend = c("1", "2","3","4"), col = c("red", "blue","green","purple"), pch = 19)

# Check the signals emitted by each type of cell (interaction of each cell type with other cells)
mat <- cellchat@net$count
par(mfrow =c(3,3),xpd=T)

for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

# Extract the inferred results
df.net <- subsetCommunication(cellchat) # Extract the cell communication prediction results as a data frame
# install.packages("DT"), DT is a tool to display data frames, which can adjust display parameters, etc.
DT::datatable(df.net)
# write.csv(df.net,'../data/cci/net/3HSMM_filtered_normalized_kmeans4.net.csv')
save.image(paste(getwd(),"/Rdata/HSMM/3finished_cellchat.RData", sep=""))

load(paste(getwd(),"/Rdata/HSMM/3finished_cellchat.RData", sep=""))
# 6.CCI————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
source("./R_tools/cci.R")
gc()
# Generate CCI
my_cci <- get_all_cci(my_data,my_meta,df.net)
save.image(paste(getwd(),"/Rdata/HSMM/4finished_cci.RData", sep=""))

load(paste(getwd(),"/Rdata/HSMM/HSMM4finished_cci.RData", sep=""))

# Output
write.csv(my_data,file = "./output/HSMM/HSMM_norm_exprs.csv")
output_meta <- my_meta
output_meta$state <- es@phenoData@data$Category
output_meta$entropy <- sc@entropies$scEntropy.bootstrap
write.csv(output_meta,file="./output/HSMM/HSMM_meta.csv")
write.csv(my_cci,file="./output/HSMM/HSMM_cci.csv")
write.csv(HSMM@assayData$exprs,file="./output/HSMM/HSMM_raw_expr.csv")
