# Load necessary library
library(readr)

# Read the file with tab-separated values
sample_ids <- read_delim("sample_names.txt", delim = "\t", col_names = FALSE)

# Convert the data into a vector of IDs
sample_ids <- as.vector(t(sample_ids))  # Transpose and convert to vector

# Function to create NEXUS tree structure with correct branch lengths
generate_nexus_tree <- function(ids) {
  # Group IDs by species code
  ga_ids <- grep("_Ga_", ids, value = TRUE)
  gm_ids <- grep("_Gm_", ids, value = TRUE)
  xc_ids <- grep("_Xc_", ids, value = TRUE)
  xr_ids <- grep("_Xr_", ids, value = TRUE)
  
  # Check for empty groups and warn the user
  if (length(ga_ids) == 0) warning("No IDs found for group '_Ga_'")
  if (length(gm_ids) == 0) warning("No IDs found for group '_Gm_'")
  if (length(xc_ids) == 0) warning("No IDs found for group '_Xc_'")
  if (length(xr_ids) == 0) warning("No IDs found for group '_Xr_'")
  
  # Create each clade with individual tip and clade branch lengths
  gm_tips <- if (length(gm_ids) > 0) paste0(gm_ids, ":0.006726", collapse = ", ") else NULL
  gm_clade <- if (!is.null(gm_tips)) paste0("(", gm_tips, "):0.006726") else NULL
  
  xr_tips <- if (length(xr_ids) > 0) paste0(xr_ids, ":0.002579", collapse = ", ") else NULL
  xr_clade <- if (!is.null(xr_tips)) paste0("(", xr_tips, "):0.002579") else NULL
  
  xc_tips <- if (length(xc_ids) > 0) paste0(xc_ids, ":0.003169", collapse = ", ") else NULL
  xc_clade <- if (!is.null(xc_tips)) paste0("(", xc_tips, "):0.000880") else NULL
  
  ga_tips <- if (length(ga_ids) > 0) paste0(ga_ids, ":0.003409", collapse = ", ") else NULL
  ga_clade <- if (!is.null(ga_tips)) paste0("(", ga_tips, "):0.000880") else NULL
  
  # Combine Xc and Ga into an inner clade (if both exist)
  inner_clade <- if (!is.null(xc_clade) && !is.null(ga_clade)) {
    paste0("(", xc_clade, ",", ga_clade, "):0.000880")
  } else if (!is.null(xc_clade)) {
    xc_clade
  } else if (!is.null(ga_clade)) {
    ga_clade
  } else {
    NULL
  }
  
  # Construct the NEXUS tree text
  tree_text <- paste0(
    "#NEXUS\n",
    "begin trees;\n",
    "  tree tree_1 = [&R] ("
  )
  
  # Add non-empty clades to the tree
  clades <- c(gm_clade, xr_clade, inner_clade)
  clades <- clades[!sapply(clades, is.null)]  # Remove NULL clades
  tree_text <- paste0(tree_text, paste(clades, collapse = ","), ");\n", "end;")
  
  return(tree_text)
}

# Generate the NEXUS tree
nexus_tree <- generate_nexus_tree(sample_ids)

# Save to file
output_file <- "Ga_tr_tree.nex"
write(nexus_tree, file = output_file)
cat("NEXUS tree saved to", output_file, "\n")

