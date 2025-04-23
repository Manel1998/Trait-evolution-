import numpy as np
import pandas as pd

# Load RNA-seq data from a CSV file.
# The CSV should have genes as rows and samples as columns.
# For example, the first column can be gene IDs, and the first row contains sample names.
data = pd.read_csv('data_GA_br_vst_normalized_counts.csv', sep='\t', header=0)


# Calculate variance for each gene across samples.
gene_variances = data.var(axis=1)

# Calculate the overall variance (flatten all values).
overall_variance = data.values.flatten().var()

# Print the results
print("Variance for each gene:")
print(gene_variances)

print("\nOverall variance across all genes and samples:")
print(overall_variance)
