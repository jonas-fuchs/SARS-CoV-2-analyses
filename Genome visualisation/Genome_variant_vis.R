## required libraries
library(ggplot2)
library(gggenes)
library(data.table)
library(stringr)
library(ggforce)
library(fs)
library(tools)

## options
zoom <- c("S") #define gene name, range (start, end) or leave empty
show_position <- c("nt") ## nt, aa, or leave empty
hide_synonymous <- T

## extract gene map from gff file 

Gff_file <-list.files(path = ".",pattern = "*.gff", recursive = F ,full.names = TRUE)

if (length(Gff_file) > 1) {
  Gff_file <- fread(file = Gff_file[1], header = F)
} else {
  Gff_file <- fread(file = Gff_file, header = F)
}


Gff_file <- Gff_file[V3 == "gene"]

Genes_extracted <- data.table()
Genes_extracted$gene <- str_extract(Gff_file$V9, "(?<=gene\\=)[:alnum:]{1,}")
Genes_extracted$start <- Gff_file$V4
Genes_extracted$end <- Gff_file$V5
Genes_extracted$molecule <- 0

for (i in 1:length(Gff_file$V7)) {
  if (Gff_file$V7[i] == "+") {
    Genes_extracted$strand[i] <- "forward"
  } else  {
    Genes_extracted$strand[i] <- "reverse"
    }
}

for (i in 1:length(Genes_extracted$strand)) {
  if (Genes_extracted$strand[i] == "forward") {
    Genes_extracted$orientation[i] <- 1
  } else {
    Genes_extracted$orientation[i] <- -1
  }
}

Genes_extracted <- as.data.frame(Genes_extracted)

## define facet zoom and pdf height if zoom is true
if (is.character(zoom) && length(zoom) >=1) {
  zoom_parameter <- c(Genes_extracted$start[Genes_extracted$gene == zoom],
                      Genes_extracted$end[Genes_extracted$gene == zoom])
  pdf_height <- 8
} else if (is.numeric(zoom)) {
  zoom_parameter <- zoom 
  pdf_height <- 8
} else {
  pdf_height <- 4
}

## load vcf file and create variant data
vcf_list <-list.files(path = ".",pattern = "*.vcf", recursive = F ,full.names = TRUE)

for (i in 1:length(vcf_list)) {
  
  Sample_number <- path_file(vcf_list[i])
  Sample_number = file_path_sans_ext(Sample_number)
  
  
  vcf_file <- fread(vcf_list[i])
  variant_data <- data.table()
  variant_data$position <- vcf_file$POS
  variant_data$af <- as.numeric(str_extract(vcf_file$INFO, "(?<=AF\\=)[^\\;]+"))
  variant_data$mutation <- paste0(vcf_file$POS, vcf_file$ALT)
  variant_data <- as.data.frame(variant_data)
  variant_data$effect <- str_extract(vcf_file$INFO, "(?<=EFF\\=)[^\\;]+")
  variant_data$effect <- str_extract(variant_data$effect, "[^\\(]+")
  
  Info_temp <- str_extract(vcf_file$INFO, "(?<=EFF\\=)[^\\;]+")
  Info_temp <- str_extract(Info_temp, "\\([^\\)]+")
  Info_temp <- str_extract(Info_temp, "[^\\(]+")
  Info_field_extracted <-t(as.data.table(str_split(Info_temp, "\\|")))
  Info_field_extracted <- as.data.frame(Info_field_extracted)
  
  variant_data$as_mutation <- Info_field_extracted$V4  
  
  if (hide_synonymous) {
    variant_data$mutation[variant_data$effect == "SYNONYMOUS_CODING"] <- NA
    variant_data$as_mutation[variant_data$effect == "SYNONYMOUS_CODING"] <- NA
  }
  
  variant_data$effect[!variant_data$effect == "SYNONYMOUS_CODING"] <- "OTHER"
  variant_data$effect[is.na(variant_data$effect)] <- "OTHER"
  
  
  ## plotting stuff
  plot <- ggplot() +
    geom_gene_arrow(data = Genes_extracted, 
                    aes(xmin = start, xmax = end, y = molecule, fill = gene),
                    arrow_body_height = grid::unit(6, "mm"),
                    arrowhead_height = grid::unit(8, "mm"),
                    show.legend = FALSE) +
    geom_gene_label(data = Genes_extracted,
                    aes(xmin = start, xmax = end, label = gene, y = molecule), 
                    align = "left") +
    scale_fill_brewer(palette = "Set3") + 
    geom_segment(data = variant_data, 
                 aes(x=position, xend=position, y=0, yend=af), 
                 color="grey") +
    geom_point(data = variant_data, 
               aes(x=position, y=af), 
               fill="darkorange2", 
               size=4, 
               pch=21, 
               color="black") +
    labs(y= "variant frequency", 
         x = "genome position")
  
  if (show_position == "nt") {
    plot <- plot +
      geom_text(data = variant_data, 
                aes(x=position, y=af, label = mutation, angle = 90), 
                nudge_y = -0.02, 
                hjust = 1, 
                size = 3)
  } else if (show_position == "as") {
    plot <- plot +
      geom_text(data = variant_data, 
                aes(x=position, y=af, label = as_mutation, angle = 90), 
                nudge_y = -0.02, 
                hjust = 1, 
                size = 3)
  }
  
  
  if (is.character(zoom) && length(zoom) >=1 | is.numeric(zoom)) {
    plot <- plot +
      facet_zoom(xlim = zoom_parameter, zoom.data = z,
                 horizontal = FALSE,
                 zoom.size=1,
                 show.area = T) + 
      theme(zoom.y = element_blank(), 
            validate = FALSE, 
            panel.grid.major.y = element_line(colour="grey80", size=0.2),
            panel.background = element_blank())
  } else {
    plot <- plot +
      theme(panel.grid.major.y = element_line(colour="grey80", size=0.2),
            panel.background = element_blank())
  }

  ## pdf creation
  pdf(paste0(Sample_number, ".pdf"), width = 21, height = pdf_height)
  
  print(plot)
  
  dev.off()
  
}