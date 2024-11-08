#-----Load libraries-----
library(NOISeq)                # Version: 2.48.0
library(SummarizedExperiment)  # Version: 1.34.0
library(TCGAbiolinks)          # Version: 2.0.0
library(tidyverse)             # Version: 1.1.4
library(biomaRt)               # Version: 2.60.1
library(edgeR)                 # Version: 4.2.1
library(EDASeq)                # Version: 2.38.0

#-----Load data-----
samples_data <- read.csv("3_Data/samples_data.tsv", 
                         sep = "\t", header = T)                     # .tsv file created in the Get_data.R script

#----------Get reads data----------
mirna_data <- GDCquery(project = "TCGA-LIHC",                        # Liver hepatocellular carcinoma
                       data.category = "Transcriptome Profiling",    # Self explanatory 
                       data.type ="miRNA Expression Quantification", # Self explanatory 
                       barcode = samples_data$Experiment_id)         # Barcodes used to filter the download files

# DO NOT RUN THIS LINE IF YOU'VE ALREADY DOWNLOADED YOUR DATA
GDCdownload(mirna_data, directory = "3_Data/")                       # Download data (will create GDCdata folder unless directory is set)
                                                                     # This folder will contain all the experiments' reads and genes (228 as 29-10-2024)

# Downloaded data has the read count, reads per million, and if they are cross mapped or not
mir_data <- GDCprepare(mirna_data, directory = "3_Data/")            # Reads the downloaded data and prepares it into an R-readable object
rownames(mir_data) <- mir_data$miRNA_ID                              # Set the miRNA id as the row names
mir_data <- mir_data[,grep("read_count",colnames(mir_data))]         # Only keep the "read_count" columns for each sample (1,881 genes, 228 samples)
colnames(mir_data) <- gsub("read_count_", "", colnames(mir_data))    # Only keep the sample name without the "read_count" part
mir_data <- mir_data[which(!rowSums(mir_data) == 0),]                # Remove miRNAs that were not detected in any single sample
dim(mir_data)                                                        # 1,490 miRNA to 228 samples (29-10-2024)

write.table(mir_data, "3_Data/miRNAseq.tsv", sep='\t', quote=F)      # Write table with reads

#----------Format dataframes and matrices----------
reduced_id <- substr(colnames(mir_data),1,19)                        # Reduced id data up to portion and analyte 
duplicated_id <- reduced_id[duplicated(reduced_id)]                  # Check for duplicated id's. None were detected (29-10-2024)

design_exp <- samples_data                                           # New object for manipulation
design_exp <- design_exp[order(match(design_exp$Experiment_id,       # Sort object by the experiment id
                                    substr(colnames(mirna_data),
                                           1,19))),]
design_exp$barcode <- colnames(mir_data)                                 # Add the reduced experiment id as a barcode
mir_data <- mir_data[,which(colnames(mir_data) %in% design_exp$barcode)] # Remove samples from expression matrix (29-10-2024: none removed)
mir_data <- mir_data[which(rowSums(mir_data) > 0),]                      # Remove features with 0 reads across all samples (29-10-2024: none removed)

#----------Annotation object----------
bmart_db <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")   # BioMart selection of ENSEMBL database. Also selected homo sapiens dataset
myannot <- getBM(attributes = c("ensembl_gene_id",                   # Get ENSEMBL gene id
                                "percentage_gene_gc_content",        # Get GC percentage for each gene
                                "mirbase_id",                        # Get the id from miRBase (most of the do not have one)
                                "start_position",                    # miRNA start position
                                "end_position",                      # miRNA end position
                                "chromosome_name"),                  # Chromosome name
                 mart = bmart_db)                                    # Object of class Mart

#----------Filtering----------
dim(myannot)                                                         # 86,627 genes
myannot <- myannot[myannot$mirbase_id%in%rownames(mir_data),]        # Only keep the annotation of the genes that are in the miRNA expression data
dim(myannot)                                                         # 1,590 genes remain
myannot$length <- myannot$end_position-myannot$start_position        # Get the transcript length
dim(myannot[!which(myannot$mirbase_id %in% rownames(mir_data)),])    # No genes are absent from the miRNA read data
myannot <- myannot[which(!duplicated(myannot$mirbase_id)),]          # Remove genes that have duplicated mirbase id (29-10-2024: 1,465 remain)
myannot <- myannot[which(!duplicated(myannot$ensembl_gene_id)),]     # Remove genes that have duplicated ensembl id (29-10-2024: 1,432 remain)
mir_data<-mir_data[which(rownames(mir_data)%in%myannot$mirbase_id),] # Remove genes that are not annotated
myannot <- myannot[order(myannot$mirbase_id),]                       # Sort by data frame by mirbase_id

#----------Annotation object----------
noiseqData <- NOISeq::readData(data = mir_data,                                                       # Create object based on counts matrix
                       factors = design_exp,                                                          # Assign conditions of interest
                       gc = myannot[,c("mirbase_id", "percentage_gene_gc_content")],                  # GC content for each feature
                       length = myannot[,c("mirbase_id", "length")],                                  # Length of every feature
                       chromosome = myannot[,c("mirbase_id", "start_position", "end_position")])      # Chromosome of each gene

mycountsbio <- dat(noiseqData, type ="countsbio", factor ="Subtype") # Generate expression data with feature annotation object

# Check expression values
explo.plot(mycountsbio, plottype = "boxplot", samples = 1:2)         # Exploratory boxplot demonstrating the expression values for each subtype

# Check for low count genes
explo.plot(mycountsbio, plottype = "barplot", samples = 1:2)         # Exploratory barplot showing CPM per subtype

# Mean of log CPM
cpm_mirna <-ggplot(data = mir_data,                                  # Exploratory histogram portrays mean log CPM to view low CPM values
                   aes(x=rowMeans(cpm(mir_data, log = T)))) + 
                   geom_histogram(color="black", fill="white", binwidth=.5) + 
                   labs(x="Mean of log CPM") +
                   geom_vline(xintercept = 0)
ggsave("2_Plots/cpm_mirna", cpm_exp, device = "png")

sum(rowMeans(cpm(mir_data, log = T)) > 0)/nrow(mir_data)             # ~28% have log CPM > 0                                       

# GC bias
GCcontent <- dat(noiseqData,  k = 0, type = "GCbias", factor = "Subtype")
GCcontent                                                            # View object for better understanding
par(mfrow=c(1,2))                                                    # Show the plot
sapply(1:2,function(x) explo.plot(GCcontent, samples = x))           # Show the GC content 

# Length bias
lenbias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "Subtype")
sapply(1:2,function(x) explo.plot(lenbias, samples = x))             # Show the length bias
par(mfrow=c(1,1)) 

# Transcript composition bias
mycd <- dat(noiseqData, type = "cd", norm = FALSE)                   # Reference sample is: TCGA-CC-5261-01A-01R-A130-13
table(mycd@dat$DiagnosticTest[,"Diagnostic Test"])                   # 227 failed (30-10-2024)
explo.plot(mycd,samples=sample(1:ncol(mir_data),10))                 # Data normalization is needed

# Check for batch effect
my_pca <- dat(noiseqData, type = "PCA", norm = F, logtransf = F)                   # Make PCA
explo.plot(my_pca, samples = c(1,2), plottype = "scores", factor = "Subtype")      # Plot to check if samples group by subtype
explo.plot(my_pca, samples = c(1,2), plottype = "scores", factor = "Variant")      # Plot to check if samples group by Variant
explo.plot(my_pca, samples = c(1,2), plottype = "scores", factor = "vital_status") # Plot to check if samples group by patients' vital status

#----------Reduce biases----------
# Filter genes with low read count
count_matrix_filtered <- filtered.data(mir_data, factor = "Subtype", norm = F,           # 227 features are to be kept (30-10-2024)
                                       depth = NULL, method = 1, cpm = 0, p.adj = "fdr") 
myannot <- myannot[myannot$mirbase_id %in% rownames(count_matrix_filtered),]             # Keep only the filtered genes (220 genes 30-10-2024)
count_matrix_filtered <- count_matrix_filtered[rownames(count_matrix_filtered) %in%      # Remove the genes that are not in my annotation
                                                 myannot$mirbase_id,]

# Create an object containing the new RNA expression set
new_exp_data <- newSeqExpressionSet(                                 # Create expression set function
  counts = as.matrix(count_matrix_filtered),                         # Gene count matrix
  featureData = data.frame(myannot, row.names=myannot$mirbase_id),   # Features of data, such as GC content, length, etc.
  phenoData = data.frame(design_exp, row.names=design_exp$barcode))  # Sample information assigned

# Correct GC and length biases
gcFull <- withinLaneNormalization(new_exp_data, "percentage_gene_gc_content", which = "full") # Correct for GC composition of reads
lFull <- withinLaneNormalization(gcFull, "length", which = "full")                            # Correct based on length

# Normalize data
normalized_count_matrix <- tmm(normCounts(lFull),                                             # Corrected read counts matrix
                               long = 1000,                                                   # No length correction is applied, already did it
                               lc = 0)   

# Create new expression set based on the corrected normalized counts
noiseqData <- NOISeq::readData(data = normalized_count_matrix, factors = design_exp)          # Reference sample is: TCGA-CC-5261-01A-01R-A130-13
mycd <- dat(noiseqData, type="cd", norm=TRUE)                                                 # Generate results of pre-processed data for exploration        
table(mycd@dat$DiagnosticTest[, "Diagnostic Test"])                                           # 217 PASSED, 10 FAILED (28-10-2024)

# Check for batch effect now that data has been pre-processed
myPCA <- dat(noiseqData, type = "PCA", norm = T, logtransf = F)                               # Perform PCA
explo.plot(myPCA, plottype = "scores", factor = "Variant")                                    # Visualize PCA (No clear clustering)
explo.plot(myPCA, plottype = "scores", factor = "Subtype")                                    # Plot to check if samples group by status

# Remove possible non identified batch effects 
full_data <- ARSyNseq(noiseqData,                                                             # Data to filter noise from
                      factor = "Subtype",                                                     # Factor of interest
                      batch = FALSE,                                                          # The factor SHALL NOT be used as the batch
                      norm = "n",                                                             # Data has already been normalized
                      logtransf = F)                                                          # False if we want to perform log-transformation of data

myPCA <- dat(full_data,                                                                       # Data to explore
             type = "PCA",                                                                    # PCA exploration
             norm = T,                                                                        # Data has already been normalized
             logtransf = T)                                                                   # Avoid performing a second log-transformation

explo.plot(myPCA, factor = "Subtype")                                                         # Visualize clustering (normal and cancer clearly separate)

#-----Final quality check-----
noiseqData <- NOISeq::readData(data = exprs(full_data),                                       # Create object based on counts matrix
                       factors = design_exp,                                                  # Assign conditions of interest
                       gc = myannot[,c("mirbase_id", "percentage_gene_gc_content")],          # GC content for each feature
                       length = myannot[,c("mirbase_id", "length")],                          # Length of every feature
                       chromosome = myannot[,c("chromosome_name",                             # Chromosome of each gene
                                               "start_position", "end_position")])

mycountsbio <- dat(noiseqData, type = "countsbio", factor = "Subtype")                        # Generate expression data with feature annotation object

# Check expression bias values
explo.plot(mycountsbio, plottype = "boxplot", samples= 1:2)                                   # Exploratory boxplot demonstrating the expression values for both types

# GC bias
GCcontent <- dat(noiseqData,type = "GCbias", factor = "Subtype", norm = T, logtransf = T)
GCcontent                                                                                     # View object for better understanding
par(mfrow=c(1,2))                                                                             # Show the plot
sapply(1:2,function(x) explo.plot(GCcontent, samples = x))                                    # Show the GC content 

# Length bias
lenbias <- dat(noiseqData, norm = T, type = "lengthbias", factor = "Subtype", logtransf = T)
sapply(1:2,function(x) explo.plot(lenbias, samples = x))                                      # Show the length bias
par(mfrow=c(1,1)) 

#----------Final matrix----------
final_mirna_mat <- exprs(full_data)                                                           # Extract final expression matrix 
write.table(final_mirna_mat, "3_Data/RNAseqnormalized.tsv", sep='\t', quote=F)                # Write table


