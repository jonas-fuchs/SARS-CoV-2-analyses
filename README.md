# SARS-CoV-2 tools for vizualization NGS and phylogenetic data

## Variant frequency plot

This plot takes two or more SNPeff annotated vcfs, extracted as tabular files to generate a heatmap of the variant frequencies: [example input](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Heatmap/example.tabular). 

```R
Rscript Heatmap/Heatmap.R
```

![Example output](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Heatmap/Heatmap.png?raw=true)

The script annotates genes and effects automatically based on the tabular files. Multiple **options** allow clustering, sorting and displaying the aminoacid changes in the mutation labels. It is also to display selectively variants that are over a desired frequency.

## Phylogenetic tree visualization

Vizualizes phylogenetic trees constructed with iqtree (tree.treefile). It takes csv tables generated by pangolin (results.csv) or nextclade (nextclade.csv) classification and colors the branches accordingly. Furthermore, the tips can be colored according to additional information (info field) and the sequences name can be selectively displayed (display_names field), both on the basis of [example input](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Tree%20clades/Metadata.xlsx).

```R
Tree clades/Tree_vis.R
```

![Example output](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Tree%20clades/Phylogenetic_tree_1.png)

Multiple **options** let you customize the output. The script allows you to root to a certain sequence and clades can be selectively displayed and colored.

![Example output](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Tree%20clades/Phylogenetic_tree_2.png)
