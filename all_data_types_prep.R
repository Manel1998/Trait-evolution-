## code for extracting the columns with specific samples 
##I am seprating them into 5 different datasets 
# Define the input and output file names
input_file <- "GA.isoform.counts.matrix"
output_file <- "GA.iso.O"
output_file <- "GA.iso.E"
output_file <- "GA.iso.T"
output_file <- "GA.iso.Te"

input_file <- "GM.isoform.counts.matrix"
output_file <- "GM.iso.O"
output_file <- "GM.iso.E"
output_file <- "GM.iso.T"
output_file <- "GM.iso.Te"

input_file <- "Xr.isoform.counts.matrix"
output_file <- "Xr.iso.O"
output_file <- "Xr.iso.E"
output_file <- "Xr.iso.T"
output_file <- "Xr.iso.Te"

input_file <- "Xc.isoform.counts.matrix"

# Define the delimiter 
delimiter <- "\t" 
# Read the file 
df <- read.table(input_file, sep = delimiter, header = FALSE, stringsAsFactors = FALSE)

# Print the first row 
cat("First row:\n")
print(df[1, ])


columns_with_O <- which(grepl("_O_", df[1, ], ignore.case = TRUE))
columns_with_E <- which(grepl("_E_", df[1, ], ignore.case = TRUE))
columns_with_T <- which(grepl("_T_", df[1, ], ignore.case = TRUE))
columns_with_Te <- which(grepl("_Te_", df[1, ], ignore.case = TRUE))

# Print identified column indices for debugging
cat("Identified columns:\n")
print(columns_with_O)
print(columns_with_E)
print(columns_with_T)
print(columns_with_Te)
# Check if there are any columns to process
if (length(columns_with_T) == 0) {
  cat("No columns found with 'O' in the first row.\n")
} else {
  # Extract the relevant columns (including all rows)
  df_selected <- df[, columns_with_T, drop = FALSE]
  
  # Print the extracted columns for debugging
  cat("Extracted columns:\n")
  print(head(df_selected))
  
  # Save the extracted columns to the output file
  write.table(df_selected, file = output_file, sep = delimiter, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("Columns containing 'O' have been saved to", output_file, "\n")
}


#extract list of genes in each table
library(dplyr)
genes_GA<- as.data.frame(pull(df, 1))
genes_GM<- as.data.frame(pull(df, 1))
genes_XR<- as.data.frame(pull(df, 1))
genes_XC<- as.data.frame(pull(df, 1))


write.table(genes_GA, "genes_GA", sep = "\t", col.names = NA, quote = FALSE,  eol = "\n")
write.table(genes_GM, "genes_GM", sep = "\t", col.names = NA, quote = FALSE,  eol = "\n")
write.table(genes_XC, "genes_XC", sep = "\t", col.names = NA, quote = FALSE,  eol = "\n")
write.table(genes_XR, "genes_XR", sep = "\t", col.names = NA, quote = FALSE,  eol = "\n")






