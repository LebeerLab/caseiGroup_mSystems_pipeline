
library(tidyverse)
library(ggtree)
library(ape)
library(phangorn)
library(adephylo)

# tree: read and root
tree = read.tree("/media/harddrive/caseiGroup/Trees/RAxML_bipartitions.caseigroup")
tree = midpoint(tree)

# genomeTable: read, fix some stuff
annotation = read_tsv("data/genomeTableWithClades.tsv", col_names = F)
names(annotation) = c("genome", "species", "clade", "strain")
annotation[annotation$genome == "AMBR2", "strain"] = "AMBR2"
speciesLevels = c("isolate", "rhamnosus", "paracasei", "casei", "zeae", "sp")
speciesLabels = c("isolate", "Lactobacillus rhamnosus", "Lactobacillus paracasei", 
                  "Lactobacillus casei", "Lactobacillus zeae", "Lactobacillus sp.")
annotation$species = factor(annotation$species, 
                            levels = speciesLevels, labels = speciesLabels)
annotation$clade = factor(annotation$clade, levels = c("cladeA", "cladeB", "cladeC"), 
                          labels = c("clade A", "clade B", "clade C"))
annotation$strain = paste(" ", annotation$strain, sep = "")
annotation$strain = paste(annotation$strain, " ", sep = "")

# define colors to use for the species and the clades
speciesColors = c("isolate" = "#ffff33","Lactobacillus rhamnosus" = "#ff7f00", 
                  "Lactobacillus paracasei" = "#f781bf", "Lactobacillus casei" = "#a65628",
                  "Lactobacillus zeae" = "#999999", "Lactobacillus sp." = "#984ea3")
cladeColors = c("clade A" = "#e41a1c", "clade B" = "#377eb8", "clade C" = "#4daf4a")

# plot and save tree
ggtree(tree, layout = "circular", size = 0.2) %<+% annotation + 
  geom_tiplab2(aes(angle = angle, col = species, label = strain), size = 1.5, 
               align = T, linesize = 0.2, offset = 0.1) +
  geom_cladelabel(node = 301, color = cladeColors["clade A"], label = "", barsize = 2, offset = 1, align = T) + 
  geom_cladelabel(node = 188, color = cladeColors["clade B"], label = "", barsize = 2, offset = 1, align = T) + 
  geom_cladelabel(node = 197, color = cladeColors["clade C"], label = "", barsize = 2, offset = 1, align = T) + 
  scale_color_manual(values = speciesColors) +
  theme(plot.background = element_blank())
ggsave("results/tree.svg")
ggsave("results/tree.emf", width = 20, height = 20, units = "cm")
ggsave("results/tree.png", bg = "transparent")

# plot and save legend for species
ggplot(data = annotation, aes(col = species)) +
  geom_point(x = 1, y = 1, size = 3) +
  scale_color_manual(values = speciesColors) 
ggsave("results/legendSpecies.svg")
ggsave("results/legendSpecies.png", bg = "transparent")

# plot and save legend for clades
ggplot(data = annotation, aes(col = clade)) +
  geom_point(x = 1, y = 1, size = 3) +
  scale_color_manual(values = cladeColors)
ggsave("results/legendClades.svg")
ggsave("results/legendClades.png")

# function to identify branches as "weak", "strong" or "tip"
getBranchTypes = function(tree) {
  Ntip = length(tree$tip.label)
  branchType = rep("weak", length(tree$edge.length))
  branchType[1:Ntip] = "tip"
  branchType[Ntip + which(as.numeric(tree$node.label) > 70)] = "strong"
  return(branchType)
}

# plot clade A subtree
tree_cladeA = extract.clade(tree, node = 301)
save(tree_cladeA, file = "Robjects/tree_cladeA.Robject")
maxDist_cladeA = distRoot(tree_cladeA, 1:length(tree_cladeA$tip.label)) %>% max()
ggtree(tree_cladeA, layout = "rectangular", size = 0.2) %<+% annotation +
  geom_tiplab(aes(col = species, label = strain), size = 2, align = T, linesize = 0.2) +
  scale_color_manual(values = speciesColors) +
  geom_point2(aes(shape = getBranchTypes(tree_cladeA)), size = 2) +
  scale_shape_manual(values = c("weak" = 32, "tip" = 32, "strong" = 16)) +
  xlim(0, maxDist_cladeA * 1.1)
ggsave(file = "results/tree_cladeA.png", width = 18, height = 15, units = "cm")

# plot clade B subtree
tree_cladeB = extract.clade(tree, node = 188)
save(tree_cladeB, file = "Robjects/tree_cladeB.Robject")
maxDist_cladeB = distRoot(tree_cladeB, 1:length(tree_cladeB$tip.label)) %>% max()
ggtree(tree_cladeB, layout = "rectangular", size = 0.2)  %<+% annotation +
  geom_tiplab(aes(col = species, label = strain), size = 2, align = T, linesize = 0.2) +
  scale_color_manual(values = speciesColors) +
  geom_point2(aes(shape = getBranchTypes(tree_cladeB)), size = 2) +
  scale_shape_manual(values = c("weak" = 32, "tip" = 32, "strong" = 16)) +
  xlim(0, maxDist_cladeB * 1.1)
ggsave(file = "results/tree_cladeB.png", width = 18, height = 3, units = "cm")

# plot clade C subtree
tree_cladeC = extract.clade(tree, node = 197)
save(tree_cladeC, file = "Robjects/tree_cladeC.Robject")
maxDist_cladeC = distRoot(tree_cladeC, 1:length(tree_cladeC$tip.label)) %>% max()
ggtree(tree_cladeC, layout = "rectangular", size = 0.2)  %<+% annotation +
  geom_tiplab(aes(col = species, label = strain), size = 2, align = T, linesize = 0.2) +
  scale_color_manual(values = speciesColors) +
  geom_point2(aes(shape = getBranchTypes(tree_cladeC)), size = 2) +
  scale_shape_manual(values = c("weak" = 32, "tip" = 32, "strong" = 16)) +
  xlim(0, maxDist_cladeC * 1.1)
ggsave(file = "results/tree_cladeC.png", width = 18, height = 24, units = "cm")
