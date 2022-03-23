# Figure1E
 
## 1 Introduction
The phylogenic tree of our SGBs. The 18,607 SGBs containing more than 80 marker genes are shown in the tree. The color of the branches indicates the classification information in the phylum level.The inner strip shows whether an SGB is unknown (novel identified) or known. The next six strip charts display those SGBs appearing in each of six host animal gut microbiomes. The outer strip chart shows the MAGs number that supports the SGBs.

## 2 Usage
### 2.1 Building a phylogenetic tree using phylophlan
For detailed steps, see https://github.com/QTP-team/build_phylogenetic_tree

### 2.2 Annotation file for the phylogenetic tree```annot.txt```
```
python graphlan.py
```

### 2.3 Visualizing phylogenetic trees using graphlan
```
python /ldfssz1/ST_META/P19Z10200N0314_LXP/2.feng_project/graphlan/graphlan/graphlan_annotate.py --annot annot.txt all_sgb_18607.tree all_sgb_18607.xml
python /ldfssz1/ST_META/P19Z10200N0314_LXP/2.feng_project/graphlan/graphlan/graphlan.py --dpi 360 --size 3.5 all_sgb_18607.xml all_sgb_18607.pdf
```