#----------Load libraries----------
library(SummarizedExperiment)  # Version: 1.34.0
library(TCGAbiolinks)          # Version: 2.32.0
library(dplyr)                 # Version: 1.1.4

#----------Expression data----------
exp_data <- GDCquery(project = "TCGA-LIHC",                        # Liver hepatocellular carcinoma
                     data.category = "Transcriptome Profiling",    # Self explanatory 
                     data.type = "Gene Expression Quantification", # Self explanatory
                     workflow.type="STAR - Counts")                # Read counts

exp_data <- getResults(exp_data)                                   # Get results from query (must be a GDCquery object)
exp_ids <- substr(exp_data$cases,1,19)                             # Extract all sample ids. Remove institute code
length(exp_ids)                                                    # There are 424 samples  (22-10-2024)

#----------Methylation data----------
meth_data <-  GDCquery(project = "TCGA-LIHC",                      # Liver hepatocellular carcinoma
                       data.category = "DNA Methylation",          # Self explanatory
                       platform="Illumina Human Methylation 450")  # Self explanatory

meth_data <- getResults(meth_data)                                 # Get results from query (must be a GDCquery object)
meth_ids <- substr(meth_data$cases,1,19)                           # Extract all sample ids. Remove institute code
length(meth_ids)                                                   # There are 1,290 samples   (22-10-2024)

#----------Micro RNA data----------
mirna_data <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma
                       data.category = "Transcriptome Profiling",     # Self explanatory
                       data.type = "miRNA Expression Quantification") # Self explanatory

mirna_data <- getResults(mirna_data)                                  # Get results from query (must be a GDCquery object)
mirna_ids <- substr(mirna_data$cases,1,19)                            # Extract all sample ids. Remove institute code
length(mirna_ids)                                                     # There are 425 samples  (22-10-2024)

#----------Extract shared data sample ids----------
samples_ids <- intersect(exp_ids, meth_ids) %>%                       # Intersect ids between the expression and methylation data ids
  intersect(., mirna_ids)                                             # Intersect the ids with the miRNA data ids
length(samples_ids)                                                   # There are 410 samples (22-10-2024)
samples_ids <- substr(samples_ids, 1, 16)                             # Format name

#----------Format data----------
tissues <- c()                                                        # Make empty vector
for(i in 1:length(samples_ids)){                                      # Iterate over every single sample id
  tissues[i] <- meth_data[which(meth_data$sample.submitter_id %in%    # Add the sample type to the vector 
                                  samples_ids[i]),]$sample_type
}

samples_data <- data.frame(cbind(samples_ids,                         # Make dataframe with first column being the whole experiment id
                                 substr(samples_ids,1,12),            # Second column has up to the patient id
                                 tissues),                            # Third column is sample types   
                           row.names = samples_ids)                   # The row names are the same as the whole experiment id

colnames(samples_data) <- c("Experiment_id",                          # Change column names
                            "Patient_id","Sample_type")               

# Rename samples' subtype from "Solid Tissue Normal" to "Normal"
samples_data[which(samples_data$Sample_type == "Solid Tissue Normal"),]$Sample_type <- 
  rep("Normal", nrow(samples_data[which(samples_data$Sample_type == "Solid Tissue Normal"),]))

samples_data$Sample_type %>% table()
# Normal     Primary Tumor     Recurrent Tumor 
# 40         367               3 

#----------HCC variants----------
# HCC is an hepatic cancer subtype by itself, therefore the subsequent divisions are variants
# Additional metadata (including variants) is necessary
sample_subtypes <- TCGAquery_subtype(tumor = "LIHC")                                 # Table with hepatocellular carcinoma variants
sum(samples_data$Patient_id %in% sample_subtypes$patient)                            # 228 samples metadata (22-10-2024)
samples_data <- samples_data[samples_data$Patient_id %in% sample_subtypes$patient,]  # Keep only annotated patient samples
table(sample_subtypes$`HCC subtypes`[sample_subtypes$patient %in%                    # View HCC subtypes (22-10-2024)
                                       samples_data$Patient_id])

# Cirrhotomimetic hepatocellular carcinoma      Clear cell hepatocellular carcinoma      Fibrolamellar carcinoma
# 3                                             10                                       4

# Lymphocyte rich hepatocellular carcinoma      Myxoid hepatocellular carcinoma          NA
# 3                                             1                                        13

# No specific subtype                           Sarcomatoid hepatocellular carcinoma     Scirrhous hepatocellular carcinoma 
# 135                                           2                                        8

# Steatohepatitic 
# 10 

samples_data$Variant <- sapply(samples_data$Patient_id,function(x)    # Assign each sample variant as metadata column
  sample_subtypes$`HCC subtypes`[sample_subtypes$patient==x])

# Change "NA" into "No specific subtype"
samples_data[which(samples_data$Variant == "NA"),]$Variant <- "No specific subtype"

# Assign normal samples as "Normal"
samples_data[which(samples_data$Sample_type == "Normal"),]$Variant <- rep("Normal", length(samples_data[which(
  samples_data$Sample_type == "Normal"),]$Variant))

# Add the subtype value even though every single sample is the same subtype of cancer
samples_data$Subtype <- rep("HCC", nrow(samples_data))

# Check for multiple samples in patients
# There are 39 samples that provide normal tissue along with primary tissue sample
table(samples_data[samples_data$Patient_id %in% samples_data$Patient_id[duplicated(samples_data$Patient_id)],2:3]) %>% table()

#----------Clinical Data----------
# There is additional clinical data available for some samples
# This type of information could be helpful to deepen research and add another layer of information
cli_data <- GDCquery_clinic("TCGA-LIHC","clinical")                   # Extract clinical data

#Keep specific clinical data
cli_data <- cli_data[, c("bcr_patient_barcode", "gender", "race", "vital_status", "ajcc_pathologic_t",
                         "ajcc_pathologic_stage", "primary_diagnosis", "prior_malignancy",
                         "ajcc_pathologic_m", "ajcc_pathologic_n",
                         "treatments_pharmaceutical_treatment_or_therapy",
                         "treatments_radiation_treatment_or_therapy")]

# Combine clinical data to expression sample data
samples_data <- cbind(samples_data,t(sapply(samples_data$Patient_id,function(x) 
  cli_data[cli_data$bcr_patient_barcode==x,2:ncol(cli_data)])))

# Basic data description
table(samples_data$Subtype)
table(as.character(samples_data$gender))
table(as.character(samples_data$race))
table(as.character(samples_data$vital_status))
table(as.character(samples_data$primary_diagnosis))
table(as.character(samples_data$ajcc_pathologic_t))
table(as.character(samples_data$ajcc_pathologic_stage))
table(as.character(samples_data$prior_malignancy))
table(as.character(samples_data$treatments_pharmaceutical_treatment_or_therapy))
table(as.character(samples_data$treatments_radiation_treatment_or_therapy))

#-----Make object into a tsv table for later pre-processing-----
samples_data <- apply(samples_data, 2, as.character)
write.table(samples_data,"3_Data/samples_data.tsv",sep='\t', row.names = F, quote = F)

