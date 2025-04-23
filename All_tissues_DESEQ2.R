
library(DESeq2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")

# Load data files
br_data <- read.csv('data_GA_br.csv', sep = '\t', row.names = 1)
ovaries_data <- read.csv('data_GA_O.csv', sep = '\t', row.names = 1)
testes_data <- read.csv('data_GA_Te.csv', sep = '\t', row.names = 1)
embryo_data <- read.csv('data_GA_E.csv', sep = '\t', row.names = 1)
tropho_data <- read.csv('data_GA_T.csv', sep = '\t', row.names = 1)

# Master colData
colData_master <- read.csv('coldata.csv')


# List of input countData DataFrames
countDataList <- list(
  data_GA_br = br_data,       
  data_GA_ov = ovaries_data,  
  data_GA_te = testes_data,   
  data_GA_em = embryo_data, 
  data_GA_tr = tropho_data    
)

colData_br <- subset(colData_master, sample == "Brain")
colData_ov <- subset(colData_master, sample == "Ovaries")
colData_te <- subset(colData_master, sample == "Testes")
colData_em <- subset(colData_master, sample == "Embryo")
colData_tr <- subset(colData_master, sample == "Trophotaenia")


colDataList <- list(
  colData_br = colData_br,
  colData_ov = colData_ov,
  colData_te = colData_te,
  colData_em = colData_em,
  colData_tr = colData_tr
)


countDataList <- lapply(countDataList, round)

dds_list <- list()

for (name in names(countDataList)) {
  # Extract corresponding colData
  colData <- colDataList[[paste0("colData_", substr(name, 9, 10))]]
  
  # Convert 'sex' to a factor and make sure "undefined" is a level
  colData$sex <- factor(colData$sex, levels = unique(c(levels(colData$sex), "undefined")))
  
  # Check if sex has at least two levels other than "undefined"
  if (length(levels(colData$sex)) > 1 && any(colData$sex != "undefined")) {
    design_formula <- ~ sex + species
  } else {
    design_formula <- ~ species
  }
  
  # Create DESeqDataSet and store it in dds_list
  dds_list[[name]] <- DESeqDataSetFromMatrix(
    countData = countDataList[[name]],
    colData = colData,
    design = design_formula
  )
}


dds_results_list <- lapply(dds_list, function(dds) {
  # Run DESeq
  DESeq(dds)
})

# Extract results names for each processed dataset
results_names_list <- lapply(dds_results_list, resultsNames)

# To view the results names for each dataset:
results_names_list



# Initialize a list to store normalized counts
normalized_counts_list <- lapply(dds_results_list, function(dds) {
  counts(dds, normalized = TRUE)
})

# Assign names to the normalized counts list using the names of dds_results_list
names(normalized_counts_list) <- names(dds_results_list)

# Apply the VST to each DESeq2 object in the list
vst_list <- lapply(dds_results_list, function(dds) {
  # Perform variance stabilizing transformation with blind=TRUE
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  # Extract the transformed counts matrix (genes x samples)
  assay(vsd)
})
# Specify the directory where to save the files
output_directory <- "~/Desktop/fish_counts"

# Loop through the normalized counts list and save each table
for (name in names(normalized_counts_list)) {
  # Define the file path
  file_path <- file.path(output_directory, paste0(name, "_normalized_counts.csv"))
  
  # Write the normalized counts to a file
  write.table(normalized_counts_list[[name]], 
              file = file_path, 
              sep = "\t", 
              col.names = TRUE,  # Include column names
              row.names = FALSE,  # Optionally include row names if needed
              quote = FALSE)
  
  # Print a message to confirm saving
  cat("Saved:", file_path, "\n")
}


for (name in names(vst_list)) {
  # Define the file path (e.g., "data_GA_br_vst_normalized_counts.csv")
  file_path <- file.path(output_directory, paste0(name, "_vst_normalized_counts.csv"))
  
  # Write the VST-transformed counts to a file
  write.table(vst_list[[name]], 
              file = file_path, 
              sep = "\t", 
              col.names = TRUE,  # Include column names
              row.names = TRUE,  # Include row names (gene IDs)
              quote = FALSE)
  
  # Print a message to confirm saving
  cat("Saved:", file_path, "\n")
}




