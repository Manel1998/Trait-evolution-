#modifying phylogenetic tree
from ete3 import Tree

# Load the tree from a Newick format file
tree = Tree("goodeid_4sp_32.nwk")

# Print the tree before modifications
print(tree.get_ascii(show_internal=True))

# Add Ga samples
ga_node = tree.search_nodes(name="XCoutput_Ga_Br_M6")[0].up
ga_samples = [
    'XCoutput_Ga_O_R1', 'XCoutput_Ga_O_R3', 'XCoutput_Ga_O_R4', 'XCoutput_Ga_O_R5',
    'XCoutput_Ga_Te_M6', 'XCoutput_Ga_Te_M7', 'XCoutput_Ga_Te_M8', 'XCoutput_Ga_Te_M9',
    'XCoutput_Ga_T_R1', 'XCoutput_Ga_T_R2', 'XCoutput_Ga_T_R4', 'XCoutput_Ga_T_R5',
    'XCoutput_Ga_E_R1', 'XCoutput_Ga_E_R3', 'XCoutput_Ga_E_R4', 'XCoutput_Ga_E_R5'
]

for sample in ga_samples:
    ga_node.add_child(name=sample)

# Add Gm samples
gm_node = tree.search_nodes(name="XCoutput_Gm_Br_M1")[0].up
gm_samples = [
    'XCoutput_Gm_E_R10', 'XCoutput_Gm_E_R13', 'XCoutput_Gm_E_R8', 'XCoutput_Gm_E_R9',
    'XCoutput_Gm_O_R11', 'XCoutput_Gm_O_R1', 'XCoutput_Gm_O_R3', 'XCoutput_Gm_O_R6',
    'XCoutput_Gm_Te_M1', 'XCoutput_Gm_Te_M7', 'XCoutput_Gm_Te_M8', 'XCoutput_Gm_Te_M9'
]

for sample in gm_samples:
    gm_node.add_child(name=sample)

# Add Xc samples
xc_node = tree.search_nodes(name="XCoutput_Xc_Br_M6")[0].up
xc_samples = [
    "XCoutput_Xc_E__R10", "XCoutput_Xc_E__R11", "XCoutput_Xc_E__R12", "XCoutput_Xc_E__R7",
    "XCoutput_Xc_O_R1", "XCoutput_Xc_O_R2", "XCoutput_Xc_O_R3", "XCoutput_Xc_O_R7",
    "XCoutput_Xc_Te_M6", "XCoutput_Xc_Te_M7", "XCoutput_Xc_Te_M8", "XCoutput_Xc_Te_M9",
    "XCoutput_Xc_T_R10", "XCoutput_Xc_T_R11", "XCoutput_Xc_T_R15", "XCoutput_Xc_T_R7"
]

for sample in xc_samples:
    xc_node.add_child(name=sample)

# Add Xr samples
xr_node = tree.search_nodes(name="XCoutput_Xr_Br_M6")[0].up
xr_samples = [
    "XCoutput_Xr_O_R1", "XCoutput_Xr_O_R3", "XCoutput_Xr_O_R4", "XCoutput_Xr_O_R5",
    "XCoutput_Xr_Te_M6", "XCoutput_Xr_Te_M7", "XCoutput_Xr_Te_M8", "XCoutput_Xr_Te_M9"
]

for sample in xr_samples:
    xr_node.add_child(name=sample)

# Print the tree after adding the new samples
print(tree.get_ascii(show_internal=True))

# Save the modified tree to a new file
tree.write(outfile="xc_all_sp_tree.tre")
