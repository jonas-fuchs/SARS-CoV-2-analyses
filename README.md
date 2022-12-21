# SARS-CoV-2 tools for vizualization NGS and phylogenetic data


## Installation guide

The scripts can be run directly from source code. Options can be adjusted in the code. Dependent libraries are installed automatically if missing.

## Genome visualisation

Creates a genome map based on [gff3 files](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Genome%20visualisation/sequence.gff3) and annotes variants from SNPeff annotated [vcf files](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Genome%20visualisation/test.vcf). Comparison to [reference vcfs](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Genome%20visualisation/ref.vcf) is possible (difference/union).

```R
Genome visualisation/genome_vis.R
```

Multiple options let you explore the genome in more detail:

![Example output](https://raw.githubusercontent.com/jonas-fuchs/SARS-CoV-2-analyses/main/Genome%20visualisation/gif.gif)


## Variant frequency plot

The script takes two or more SNPeff annotated vcfs, extracted as tabular files to generate a heatmap of the variant frequencies: [example input](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Heatmap/example_input.rar). 

```R
Rscript Heatmap/Heatmap.R
```

![Example output](https://raw.githubusercontent.com/jonas-fuchs/SARS-CoV-2-analyses/main/Heatmap/Heatmap.png)

The script annotates genes and effects on the genes automatically based on the tabular files. Multiple **options** allow clustering, sorting and displaying the aminoacid changes in the mutation labels. The scripts can also selectively display variants that are over a desired frequency. Adjust options within the script if desired, otherwise default settings are used.

## Phylogenetic tree visualization

Vizualizes phylogenetic trees constructed with iqtree [Tree.treefile](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Tree%20clades/Tree.treefile). The script expects csv tables generated by pangolin (results.csv) or nextclade ([nextclade.csv](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Tree%20clades/nextclade.csv)) classification and colors the branches accordingly. Furthermore, the tips can be colored according to additional information (info field) and the sequences name can be selectively displayed (display_names field), both on the basis of [example input](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Tree%20clades/Metadata.xlsx).
Adjust options within the script if desired, otherwise default settings are used.

```R
Tree clades/Tree_vis.R
```

![Example output](https://raw.githubusercontent.com/jonas-fuchs/SARS-CoV-2-analyses/main/Tree%20clades/Phylogenetic_tree_1.png)

The script allows you to root to a certain sequence and clades can be selectively displayed and colored.

![Example output](https://raw.githubusercontent.com/jonas-fuchs/SARS-CoV-2-analyses/main/Tree%20clades/Phylogenetic_tree_2.png)

