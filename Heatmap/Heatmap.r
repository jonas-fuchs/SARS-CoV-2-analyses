

# libary import
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(fs)
library(tools)
library(stringr)
library(data.table)


# options 
name_plot <- "name"               ## name your plot
frequency <- 0.1                   ## adjust the variant frequency
display_as <- T                     ## display also AS changes 
percent <- F                        ## should variants be shown as percent?
rounded <- F                        ## should rounded values be displayed in the map
multiplier_width <- 0.6             ## pdf weight multiplier - adjust if needed
multiplier_height <- 0.6            ## pdf height multiplier - adjust if needed
color_gene_annotation <- c("Set3")  ## adjust color of gene annotation ("Set2" or "Paired")
date <- F                           ## do the file names contain a date (format: dd.mm.yyyy) and you want to sort
number <- T                         ## do the file names contain a number and you want to sort
clustering <- F                     ## should the samples be clustered?
clustering_method <- "ward.D2"      ## what clustering method -> see ?hclust for further information
number_of_clusters <- 7             ## do you assume a particular amount of clusters?


# data preparation

## import data
files<-list.files(path = ".",pattern = "*.tabular", recursive = F ,full.names = TRUE)

## only in case it contains a date:
if(date) {
  files<-files[order(as.Date(regmatches(files, regexpr("\\d{2}[[:punct:]]\\d{2}[[:punct:]]\\d{4}", files)), format="%d.%m.%y"))]
}

## create data.frame for first dataset
variants<-read.table(files[1], header = T, sep = "\t")
final<- paste0(variants$POS, variants$ALT, if(display_as){paste0(" ", "(", variants$EFF....AA, ")")})
final<- data.frame(cbind(final, variants$AF))
Sample_number<-path_file(files[1])
Sample_number = file_path_sans_ext(Sample_number)
colnames(final) <- c("Mutation", Sample_number)


## loop over all other datasets and join by Mutation
for(i in 2:length(files)) {
  variants<-read.table(files[i], header = T, sep = "\t")
  varpos<- paste0(variants$POS, variants$ALT, if(display_as){paste0(" ", "(", variants$EFF....AA, ")")})
  varpos<- cbind(varpos, variants$AF)
  Sample_number<-path_file(files[i])
  Sample_number = file_path_sans_ext(Sample_number)
  colnames(varpos) <- c("Mutation", Sample_number)
  final <- full_join(final, varpos, by="Mutation", copy=T)
}

## sort
final<-final[str_order(final$Mutation, numeric = T),]

## set row names
row.names(final)<-final$Mutation
final<-final[,-1]


## order samples after number in name
if(number){final<-final[,str_order(colnames(final), numeric = T)]}

## adjust the variant frequency:
final[final < frequency] <- NA

final <- final[rowSums(is.na(final)) !=ncol(final), ]
final <- t(final)
final[is.na(final)] <- 0
class(final) <- "numeric"

# convert values to percent
if(percent){final<-final*100}

# show rounded values in pheatmap
if(rounded){final_rounded <- round(final, 0)} else {final_rounded <- F}

# add annotations

## readout annotations
ann_final <- read.table(files[1], header = T, sep = "\t") 
ann_final <- cbind(paste0(ann_final$POS, ann_final$ALT, if(display_as){paste0(" ", "(", ann_final$EFF....AA, ")")}), ann_final$EFF....EFFECT, ann_final$EFF....GENE)


for(i in 2:length(files)) {
  ann <- read.table(files[i], header = T, sep = "\t")
  ann <- cbind(paste0(ann$POS, ann$ALT, if(display_as){paste0(" ", "(", ann$EFF....AA, ")")}), ann$EFF....EFFECT, ann$EFF....GENE)
  ann_final <- rbind(ann_final, ann)
}

ann_final<-data.frame(unique(ann_final))

## apply frequency filter
ann_final <- ann_final[ann_final$X1 %in% colnames(final),]
## sort
ann_final <- ann_final[str_order(ann_final$X1, numeric = T),]

## set rownames
row.names(ann_final) <- ann_final$X1
ann_final <- ann_final[,-1]

## set colnames
colnames(ann_final) <- c("effect", "gene")

## rename annotations
ann_final$effect <- sub("^$", "non-coding", ann_final$effect)
ann_final$gene <- sub("^$", "NCR", ann_final$gene)

ann_final$effect[ann_final$effect=="NON_SYNONYMOUS_CODING"] <- "non-syn"
ann_final$effect[ann_final$effect=="SYNONYMOUS_CODING"] <- "syn"
ann_final$effect[ann_final$effect=="CODON_CHANGE_PLUS_CODON_DELETION"] <- "deletion"
ann_final$effect[ann_final$effect=="CODON_DELETION"] <- "deletion"
ann_final$effect[ann_final$effect=="FRAME_SHIFT"] <- "frame shift"
ann_final$effect[ann_final$effect=="STOP_GAINED"] <- "stop gained"
ann_final$effect[ann_final$effect=="CODON_DELETION"] <- "deletion"
ann_final$effect[ann_final$effect=="FRAME_SHIFT+STOP_GAINED"] <- "stop gained"
ann_final$effect[ann_final$effect=="CODON_CHANGE_PLUS_CODON_INSERTION"] <- "insertion"
ann_final$effect[ann_final$effect=="FRAME_SHIFT+STOP_LOST+SPLICE_SITE_REGION"] <- "frame shift"
ann_final$effect[ann_final$effect=="INSERTION"] <- "insertion"
ann_final$effect[ann_final$effect=="START_LOST"] <- "non-syn"

# utomatically determine gaps for the heatmap

gap_vector <- c()

for (i in 2:length(ann_final$gene)){
  if (ann_final$gene[i]!=ann_final$gene[i-1]){
    gap_vector <- c(gap_vector, i-1)
  }
}

# colormanagement

## colormangment heatmap
my_colors <- colorRampPalette(c("grey93","brown","black"))

## colormanagement annotations (genes)
count <- length(unique(ann_final$gene))
gene_color <- c(brewer.pal(color_gene_annotation, n=count))
names(gene_color) = unique(ann_final$gene)

## colormanagement annotations (effect)

colors <- c()

if(c("non-coding") %in% ann_final$effect){
  colors <- c(colors, "white")
}
if(c("syn") %in% ann_final$effect){
  colors <- c(colors, "green")
}
if(c("non-syn") %in% ann_final$effect){
  colors <- c(colors, "orange")
}
if(c("deletion") %in% ann_final$effect){
  colors <- c(colors, "red")
}
if(c("frame shift") %in% ann_final$effect){
  colors <- c(colors, "black")
}
if(c("stop gained") %in% ann_final$effect){
  colors <- c(colors, "grey")
}
if(c("insertion") %in% ann_final$effect){
  colors <- c(colors, "blue")
}


all_colors<-data.frame(color=c("white", "green", "orange", "red", "black", "grey", "blue"), 
                       name = c("non-coding", "syn", "non-syn", "deletion", "frame shift", "stop gained", "insertion"))

subset_colors<-subset(all_colors, color %in% colors)

effect_color <- subset_colors$color
names(effect_color) = subset_colors$name

color_list <- list(gene_color = gene_color, effect_color = effect_color)
names(color_list) <- c("gene", "effect")


# visualize heatmap

## The pdf scales with the number of variants and samples
pdf(paste0(name_plot, "(", frequency, ")", ".pdf"), width = multiplier_width*ncol(final), height = multiplier_height*nrow(final))

pheatmap(final, 
         color = my_colors(100),
         cellwidth = 20,
         cellheight = 20,
         clustering_method = clustering_method,
         cluster_rows = clustering, 
         cluster_cols = F,
         cutree_rows = number_of_clusters,
         annotation_col = ann_final,
         annotation_colors = color_list,
         gaps_col = gap_vector,
         display_numbers = final_rounded
         )

dev.off()
