# libraries
library(ggplot2)
library(ggtree)
library(treeio)
library(data.table)
library(xlsx)
library(RColorBrewer)
library(dplyr)

# options
## provide metadata table as xlsx with seq, info and display_names columns and results of nextclade or pangolin as csv

option_root <- "NC_045512"
option_clade <- "P" ## choose if Nextclade (N) or Pangolin (P) nomenclature should be displayed
option_layout <- "fan" ## see ggtree layouts
option_text <- "" ##overwrite the display_names with no (N) or all sequences names (A), if not leave empty
option_color_tip <- F ##show info as tip colors
option_font_size <- 4
option_rename_and_recolor_tip <- T
option_rename_info <- "Info" ##rename legend
option_tip_color <- "Dark2" #RColorBrewer palette
option_pdf_name <- "Test" ##for pdf
option_height <- 10 ##for pdf
option_width <- 12 ##for pdf
option_only_specific_clade <- T
option_spec_clades <- c("B.1.1.7", "B.1.221")
option_color_branches_manual <- T
option_branch_color <- c("black", "red", "green")  ##first color is for uncolored branches, provide equal number of colors to number of clades

# format metadata

# nextclade/pangolin data

## read in nextclade or pangolin annotations
clade <- data.table()

if(option_clade == "N"){
  nextclade <- fread("nextclade.csv") 
  clade$seq <- nextclade$seqName
  clade$clade <- nextclade$clade
} else if(option_clade == "P") {
  pangolin <- fread("results.csv")
  clade$seq <- pangolin$`Sequence name`
  clade$clade <- pangolin$Lineage
} else {
  print("ERROR: N for Nextclade, P for pangolin")
}

## display secific lineages/clades only
if(option_only_specific_clade){
  clade_temp <- clade[clade == option_spec_clades[1]]
  
  for(i in 2:length(option_spec_clades)){
    clade_temp_temp <- clade[clade == option_spec_clades[i]]
    clade_temp <- rbind(clade_temp, clade_temp_temp)
  }
  clade$clade <- NULL
  clade <- full_join(clade, clade_temp, by = "seq")
  clade$clade[is.na(clade$clade)] <- "0"
}

## create list for OTU assignment
if(option_clade == "N" || option_clade == "P"){
  clade_names <- unique(clade$clade)
  cls <- list()
  
  for(i in 1:length(clade_names)){
    clade_temp <- clade[clade == clade_names[i]]
    clade_temp <- clade_temp$seq
    cls_temp <- list("clade_name_temp"=clade_temp)
    names(cls_temp) <- clade_names[i]
    cls <- append(cls, cls_temp)
  }
}

## read in iq tree and root (*.treefile)
tree <- read.iqtree("Tree.treefile")
tree <- as.phylo(tree)
tree <- root(tree, outgroup = option_root)

## group tree with OTU
if(option_clade == "N" || option_clade == "P"){
  tree <- groupOTU(tree, cls, "Clade", overlap = "overwrite")
}

#metadata for tips and heatmap

metadata <- read.xlsx("Metadata.xlsx", sheetName = "Sheet1", header = T)

if(option_text == "A"){
  metadata$display_names <- metadata$seq
} else if(option_text == "N"){
  metadata$display_names <- NA
} else {
}

# colors
tip_color <- brewer.pal(length(unique(metadata$info)), option_tip_color)

#tree visualization
if(option_clade == "N" || option_clade == "P"){
  tree_vis <- ggtree(tree, aes(color = Clade), layout = option_layout, show.legend = T)
} else {
  tree_vis <- ggtree(tree, layout = option_layout, show.legend = T)
}

if(option_color_tip){
  tree_vis <- tree_vis %<+% metadata +
    geom_tippoint(shape = 21, aes(fill = info) , size = 2, color = "black")
  
} else {
  tree_vis <- tree_vis %<+% metadata +
    geom_tippoint(shape = 21, fill = "grey", size = 2, color = "black")
}

if(option_rename_and_recolor_tip){
  tree_vis <- tree_vis + 
    scale_fill_manual(values = tip_color, name = option_rename_info)
}

tree_vis <- tree_vis + 
  geom_tiplab(aes(label = display_names), color = "black", hjust = -.1, size = option_font_size) +
  geom_treescale(fontsize = option_font_size, linesize = 0.2)

if(option_color_branches_manual){
  tree_vis <- tree_vis +
    scale_color_manual(values = option_branch_color)
}

# plot Tree
pdf(paste0(option_pdf_name, ".pdf"), height = option_height, width = option_width)

plot(tree_vis)

dev.off()
