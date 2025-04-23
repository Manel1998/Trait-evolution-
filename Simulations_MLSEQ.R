library(ape)
library(compcodeR)
library(DESeq2)


my_tree <- read.tree('tree.nwk')
# Read the tree into R

# Confirm number of tips (should be 32)
print(Ntip(my_tree))  # Expected: 32

# Convert the tree to an ultrametric tree using chronos
tree_ultra <- chronos(my_tree, lambda = 1)

# Confirm ultrametric status
print(is.ultrametric(tree_ultra))  # Expected: TRUE


# (Optional) Visualize the tree with condition labels
plot(tree_ultra, label.offset = 0.01)
tiplabels(pch = 19, col = c("#D55E00", "#009E73")[id_cond])




# Assign each tip as a unique species (i.e. 32 unique species)
id_species <- factor(tree_ultra$tip.label)
names(id_species) <- tree_ultra$tip.label

# Now assign conditions in a balanced design (16 per condition)
id_cond <- factor(rep(1:2, each = length(tree_ultra$tip.label) / 2))
names(id_cond) <- tree_ultra$tip.label

# Check counts:
print(length(unique(id_species)))  # Should now return 32
print(table(id_species))           # Each species (tip) occurs once
print(table(id_cond))              # 16 per condition

# Then run generateSyntheticData() with samples.per.cond = 16:
samples_per_cond <- length(tree_ultra$tip.label) / 2  # 32 / 2 = 16
# Set random seed for reproducibility
set.seed(12890926)
# Define simulation parameters
n_genes <- 22183                # Number of genes
n_diffexp <- 2200               # Number of differentially expressed genes
seq_depth <- 1e7
effect_size <- 3
prop_var_tree <- 0.9           # Tree accounts for 90% of variance


neutral <- generateSyntheticData(
  dataset = "Neutral_Evolution",
  n.vars = n_genes,
  samples.per.cond = samples_per_cond,  # 16
  n.diffexp = 0,                        # No DE genes
  repl.id = 1,
  seqdepth = seq_depth,
  effect.size = 0,
  fraction.upregulated = 0,
  output.file = "Neutral_Evolution_repl1.rds",
  ## Phylogenetic parameters
  tree = tree_ultra,
  id.species = id_species,
  id.condition = id_cond,
  model.process = "BM",  # Brownian Motion
  prop.var.tree = prop_var_tree,
  lengths.relmeans = "auto",
  lengths.dispersions = "auto"
)

summarizeSyntheticDataSet(data.set = "Neutral_Evolution_repl1.rds", 
                         output.filename = "Neutral_Evolution_datacheck.html")

# Define simulation parameters for positive selection
n_diffexp_pos <- 5000      # Increase the number of DE genes under positive selection
effect_size_pos <- 5       # A larger effect size to mimic strong directional selection
fraction_up_pos <- 0.8     # Higher proportion of upregulated genes in the favored condition

positive_sel <- generateSyntheticData(
  dataset = "Positive_Selection",
  n.vars = n_genes,
  samples.per.cond = samples_per_cond,  # 16 in our example
  n.diffexp = n_diffexp_pos,            # More DE genes
  repl.id = 1,
  seqdepth = seq_depth,
  effect.size = effect_size_pos,        # Strong effect size
  fraction.upregulated = fraction_up_pos,
  output.file = "Positive_Selection_repl1.rds",
  ## Phylogenetic parameters
  tree = tree_ultra,
  id.species = id_species,  # Here, note: for simulation, each sample must be a unique species
  id.condition = id_cond,
  model.process = "BM",     # Replace with "OU" or another model if available for directional change
  prop.var.tree = prop_var_tree,
  lengths.relmeans = "auto",
  lengths.dispersions = "auto"
)


summarizeSyntheticDataSet(data.set = "Positive_Selection_repl1.rds", 
                          output.filename = "Positive_selection_datacheck.html")



stablising_sel <- generateSyntheticData(
  dataset = "Stabilising_Selection",
  n.vars = n_genes,
  samples.per.cond = samples_per_cond,  # 16 in our example
  n.diffexp = 2000,            # More DE genes
  repl.id = 1,
  seqdepth = seq_depth,
  effect.size = 2,        # Strong effect size
  fraction.upregulated = fraction_up_pos,
  output.file = "Stabilising_Selection_repl1.rds",
  ## Phylogenetic parameters
  tree = tree_ultra,
  id.species = id_species,  # Here, note: for simulation, each sample must be a unique species
  id.condition = id_cond,
  model.process = "OU",
  selection.strength = 0.2,
  prop.var.tree = 0.3,
  lengths.relmeans = "auto",
  lengths.dispersions = "auto"
)


summarizeSyntheticDataSet(data.set = "Stabilising_Selection_repl1.rds", 
                          output.filename = "Stabilising_Selection_datacheck.html")

divergent_evolution <- generateSyntheticData(
  dataset = "divergent_evolution ",
  n.vars = n_genes,
  samples.per.cond = samples_per_cond,  # 16 in our example
  n.diffexp = 4000,            # More DE genes
  repl.id = 1,
  seqdepth = seq_depth,
  effect.size = 2,        # Strong effect size
  fraction.upregulated = fraction_up_pos,
  output.file = "divergent_evolution_repl1.rds",
  ## Phylogenetic parameters
  tree = tree_ultra,
  id.species = id_species,  # Here, note: for simulation, each sample must be a unique species
  id.condition = id_cond,
  model.process = "BM",     # Replace with "OU" or another model if available for directional change
  prop.var.tree = 0.3,
  lengths.relmeans = "auto",
  lengths.dispersions = "auto"
)


summarizeSyntheticDataSet(data.set = "divergent_evolution_repl1.rds", 
                          output.filename = "divergent_evolution_datacheck.html")


BiocManager::install(c("MLseq", "DESeq2", "caret"))
library(caret)
library(MLSeq)

pos_counts <- positive_sel@count.matrix
neu_counts <- neutral@count.matrix
sta_counts <- stablising_sel@count.matrix
div_counts  <- divergent_evolution@count.matrix


coldata_all <- read.csv("../coldata_ga.csv")
coldata_brain <- subset(coldata_all, sample=='Brain')
# Create DESeq2 object
dds_sim <- DESeqDataSetFromMatrix(
  countData = div_counts,
  colData = coldata_brain,
  design = ~ sex + species
)


dds_sim <- DESeq(dds_sim)


vst_sim <- varianceStabilizingTransformation(dds_sim, blind = FALSE)
vst_counts_sim <- assay(vst_sim)  

write.csv(vst_counts_sim, 'sim_div_evo_counts_vst_ga_br.csv')



