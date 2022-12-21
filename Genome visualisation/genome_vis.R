##### load library #####

list.of.packages <- c("data.table", "stringr", "RColorBrewer", "fs", "tools", "ggplot2", "BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!"trackViewer" %in% installed.packages()[,"Package"]){
  BiocManager::install("trackViewer")
}

lapply(c(list.of.packages, "trackViewer"), require, character.only = TRUE)

#### working dir ####

wd = "C:/Users/"
setwd(wd)

#### options ####

# see brewer pals
color_set = "GnBu" 
# false or vector with "non-syn", "other", "del", "syn", "frameshift", "ins", "stop" 
# combinations possible
hide_mutations = F
# intersect mutations present in another vcf: FLASE or relpath/*file_name*.vcf
ref_vcf = F
# intersect type: union, diff
intersect_type = "union"
# show nt changes
show_nt = F
# show amino acid changes, overwrites show_nt
show_as = F
# variant frequency threshold
freq_thres = 0.95
# false or provide nt positions c(start, stop)
show_part = F
# false or gene name  - will overwrite show_part
show_gene = F
# print density plot of mutation counts
print_dens = F

##### script #####

## define variables and colors
mutation_type_user = c("syn", "non-syn", "del", 
                       "other", "frameshift", "ins", 
                       "stop")
mutation_type = c("SYNONYMOUS_CODING","NON_SYNONYMOUS_CODING","DELETION",
                  "OTHER", "FRAME_SHIFT", "INSERTION", 
                  "STOP")
# colors for mutation type
all_colors = c("#9BC0E6", "#FAA41A", "#ED1E24", 
               "grey95", "olivedrab", "grey30", 
               "black")

## define features
# read in GFF
Gff_file =list.files(path = ".",pattern = "*.gff", recursive = F ,full.names = TRUE)
Gff_file = fread(file = Gff_file, header = F)

# extract genes and relevant ranges
Gff_file_temp = Gff_file[V3 == "gene"]
Genes_extracted = data.table()
Genes_extracted$gene = str_extract(Gff_file_temp$V9, "(?<=gene\\=)[:alnum:]{1,}")
Genes_extracted$start = Gff_file_temp$V4
Genes_extracted$end = Gff_file_temp$V5
Genes_extracted$length = Genes_extracted$end - Genes_extracted$start + 1

#adjust height of features
Genes_extracted$height = seq_len(nrow(Genes_extracted)) %% 2 
Genes_extracted$height[Genes_extracted$height == 1] = 8.5
Genes_extracted$height[Genes_extracted$height == 0] = 10

# define colors
all_pal = brewer.pal.info
pal_length = all_pal$maxcolors[rownames(all_pal)==color_set]

if (length(Genes_extracted$gene) > pal_length) {
  coul = brewer.pal(pal_length, color_set)
  coul = colorRampPalette(coul)(length(Genes_extracted$gene))
} else {
  coul = brewer.pal(length(Genes_extracted$gene), color_set)
}

# define genome name and range
genome_name = unique(Gff_file$V1)

if (isFALSE(show_part)) {
  genome_start = min(Gff_file$V4)
  genome_stop = max(Gff_file$V5)
} else {
  genome_start = show_part[1]
  genome_stop = show_part[2]
}

if (!isFALSE(show_gene)){
  genome_start = Genes_extracted$start[Genes_extracted$gene==show_gene]
  genome_stop = Genes_extracted$end[Genes_extracted$gene==show_gene]
}

genome_range = seq(genome_start, genome_stop, by= genome_stop/6)

# adjust x axis ranges depending on the size of the genome/zoom
if (max(genome_stop-genome_start) > 10) {
  genome_range = round(genome_range, 0)
}
if (max(genome_stop-genome_start) > 100) {
  genome_range = round(genome_range, -1)
}
if (max(genome_stop-genome_start) > 10000) {
  genome_range = round(genome_range, -2)
}
if (max(genome_stop-genome_start) > 100000) {
  genome_range = round(genome_range, -3)
}

genome_range = unique(c(genome_start, genome_range[genome_range>0], genome_stop))


# generate features for the plot
features = GRanges(genome_name, IRanges(start = Genes_extracted$start,
                                                 end = Genes_extracted$end,
                                                 width=Genes_extracted$length,
                                                 names=Genes_extracted$gene))
features$fill = coul
features$height = unit(Genes_extracted$height, "mm")

## define variants
vcf_list = list.files(path = ".",pattern = "*.vcf", recursive = F ,full.names = TRUE)

# read in ref vcf and generate mutation list to delete
if (!isFALSE(ref_vcf)){
  comparison_vcf = fread(ref_vcf)
  comparison_vcf$af = as.numeric(str_extract(comparison_vcf$INFO, "(?<=AF\\=)[^\\;]+"))
  comparison_vcf = comparison_vcf[comparison_vcf$af > freq_thres,]
  mutations_present = paste0(comparison_vcf$POS, comparison_vcf$ALT)
}

for (i in 1:length(vcf_list)) {
  
  # read in vcf
  vcf_file = fread(vcf_list[i])
  # create variant table
  variant_data = data.table()
  variant_data$position = vcf_file$POS
  variant_data$af = as.numeric(str_extract(vcf_file$INFO, "(?<=AF\\=)[^\\;]+"))
  variant_data$mutation = paste0(vcf_file$POS, vcf_file$ALT)
  variant_data = as.data.frame(variant_data)
  variant_data$effect = str_extract(vcf_file$INFO, "(?<=EFF\\=)[^\\;]+")
  variant_data$effect = str_extract(variant_data$effect, "[^\\(]+")
  # extract data from info field
  Info_temp = str_extract(vcf_file$INFO, "(?<=EFF\\=)[^\\;]+")
  Info_temp = str_extract(Info_temp, "\\([^\\)]+")
  Info_temp = str_extract(Info_temp, "[^\\(]+")
  Info_field_extracted = t(as.data.table(str_split(Info_temp, "\\|")))
  Info_field_extracted = as.data.frame(Info_field_extracted)
  variant_data$as_mutation = Info_field_extracted$V4
  # intersect with user defined vcf
  if (!isFALSE(ref_vcf)){
    if (intersect_type == "diff") {
      variant_data = variant_data[!variant_data$mutation %in% mutations_present,]
    }
    if (intersect_type == "union") {
      variant_data = variant_data[variant_data$mutation %in% mutations_present,]
    }
  }
  # sanitize data
  variant_data$as_mutation[is.na(variant_data$as_mutation)] = " "
  variant_data$effect[is.na(variant_data$effect)] = "OTHER"
  variant_data$effect[str_detect(variant_data$effect, "DELETION")] = "DELETION"
  variant_data$effect[str_detect(variant_data$effect, "FRAME_SHIFT")] = "FRAME_SHIFT"
  variant_data$effect[str_detect(variant_data$effect, "INSERTION")] = "INSERTION"
  variant_data$effect[str_detect(variant_data$effect, "CODON_CHANGE_PLUS_CODON_")] = "INSERTION"
  variant_data$effect[str_detect(variant_data$effect, "STOP")] = "STOP"
  variant_data$effect[str_detect(variant_data$effect, "START_LOST")] = "NON_SYNONYMOUS_CODING"
  variant_data$effect[str_detect(variant_data$effect, "GENE_FUSION")] = "NON_SYNONYMOUS_CODING"
  variant_data$as_mutation[variant_data$effect == "SYNONYMOUS_CODING"] = ""
  
  # filter by variant threshold
  variant_data = variant_data[variant_data$af > freq_thres,]
  # hide user defined mutation types
  if (!isFALSE(hide_mutations)) {
    if (any(hide_mutations %in% mutation_type_user)) {
      for (j in 1:length(mutation_type_user)) {
        if (mutation_type_user[j] %in% hide_mutations) {
          variant_data$mutation[variant_data$effect == mutation_type[j]] <- NA
        }
      }
      variant_data = variant_data[!is.na(variant_data$mutation),]
    }
  }

  # create legend
  labels = unique(variant_data$effect)
  label_colours = all_colors[is.element(mutation_type, labels)]
  labels = mutation_type[is.element(mutation_type, labels)]
  legend = label_colours
  names(legend) = labels 
  legend = list(labels=labels, 
                 col="grey30", 
                 fill=label_colours)
  
  # define colors in df
  for (j in 1:length(mutation_type)) {
    variant_data$colour[variant_data$effect == mutation_type[j]] = all_colors[j]
  }

  # create variants info for the plot
  if (show_nt) {
    frequency = GRanges(genome_name, 
                        IRanges(variant_data$position, 
                                width=1, 
                                names = variant_data$mutation))
  }
  if (show_as) {
    frequency = GRanges(genome_name, 
                        IRanges(variant_data$position, 
                                width=1, 
                                names = variant_data$as_mutation))
  }
  if (!any(c(show_as, show_nt))) {
    frequency = GRanges(genome_name, 
                        IRanges(variant_data$position, 
                                width=1))
  }
  
  frequency$border = "grey30"
  frequency$score = variant_data$af*100
  frequency$color = variant_data$colour
  
  # plot and create pdf
  Sample_number <- path_file(vcf_list[i])
  Sample_number = file_path_sans_ext(Sample_number)
  
  pdf(paste0(Sample_number, ".pdf"), 
      width = 21, 
      height = 5)
  
  print(lolliplot(frequency, 
            features,
            legend=legend,
            ranges = GRanges(genome_name, 
                             IRanges(genome_start, 
                                     genome_stop)),
            cex=.8,
            xaxis = genome_range
  ))
  
  dev.off()
  
  # create density plot for mutation counts
  if (print_dens) {
    bin = 10
    mutation_count_bin = data.table()
    
    mutation_count_bin$START = seq(1, genome_stop, by=bin)
    mutation_count_bin$STOP = c(seq(bin, genome_stop, by=bin),genome_stop)
    
    for (j in 1:length(mutation_count_bin$START)) {
      temp = variant_data[variant_data$position > mutation_count_bin$START[j],]
      temp = temp[temp$position < mutation_count_bin$STOP[j],]
      mutation_count_bin$COUNT[j] = nrow(temp)
    }
    
    mutation_count_bin = mutation_count_bin[!mutation_count_bin$COUNT == 0,]
    mutation_count_bin = mutation_count_bin[,c(1,3)]
    
    
    pdf(paste0(Sample_number,"_density", ".pdf"), 
        width = 5, 
        height = 3)
    
    print(
      ggplot(mutation_count_bin, aes(x = START)) + 
        geom_histogram(aes(y = after_stat(density)), binwidth = 500, colour = "black", fill= "white") +
        geom_density(lwd = 0.5, colour = "grey 30",
                     fill = "grey 70", alpha = 0.25) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_classic() +
        labs(y= "mutation density", x = "genome position") +
        scale_x_continuous(breaks = genome_range)
    )
    dev.off()
  }

}
