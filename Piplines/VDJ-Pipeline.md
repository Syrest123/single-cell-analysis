# VDJ Pipeline
This pipeline uses basic bioconductor packages to generate visuals for VDJ data. This VDJ data was combined with the GEX dataset to produce a .rds object which was used in this analysis.


## Installing required packages
```{r}
install.packages("fmsb")
BiocManager::install("immunarch")
BiocManager::install("ComplexHeatmap")
BiocManager::install("networkD3")
```

## Loading libraries
```{r}
# Getting the working directory
getwd()

# Set your working directory
setwd()

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(immunarch)
  library(ggpubr)
  library(tidyverse)
  library(pheatmap)
  library(UpSetR)
  library(RColorBrewer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(circlize)
  library(fmsb)
  library(ComplexHeatmap)
  library(circlize)
  library(igraph)
  library(ggraph)
})
```

## Reading in the .rds object
Object can be downloaded from the *data* directory
```{r}
# Load your Seurat object with HTO data
bcells <- readRDS("<< path to the .rds object >>")

# List of B cell identities you want to keep
bcell_types <- c("Class-switched memory B-cells", 
                 "Memory B-cells", 
                 "naive B-cells", 
                 "Plasma cells")

# Subset the Seurat object
bcells_only <- subset(bcells, 
                     subset = SingleR.labels %in% bcell_types)


# Verify HTO Metadata by checking for available metadata columns
head(bcells@meta.data)
head(bcells_only@meta.data)


# Check current hashtag labels
unique(bcells_only$hash.ID) 

# Rename the factor levels 'Hashtag-? to Patient??'
bcells_only$hash.ID <- factor(bcells_only$hash.ID,
                         levels = c("Hashtag-1", "Hashtag-2", "Hashtag-3", "Hashtag-4"),
                         labels = c("Patient1_Early", "Patient1_Late", "Patient2_Early", "Patient2_Late"))

# Verify changes
unique(bcells_only$hash.ID)
```

## Checking for cells with clonotypes
Only B cells have V(D)J rearrangements.
```{r}
#  Plot clusters vs. clonotype presence
bcells_only@meta.data %>%
  ggplot(aes(x = seurat_clusters, fill = is.na(raw_clonotype_id))) +
  geom_bar(position = "fill") +
  theme_classic() + 
  labs(title = "Missing Clonotypes by Cluster", y = "Proportion")

# Save the visual 
ggsave("<<path to visual>>")
```
<img width="663" height="411" alt="upload_c31b5e9a1d3f4537e56a64fad95aadb6" src="https://github.com/user-attachments/assets/2016b112-ffbd-4c9b-ba2c-259a53c3ecaf" />

```{r}
# Plot cell types vs. clonotype presence
bcells_only@meta.data %>%
  ggplot(aes(x = SingleR.labels, fill = is.na(raw_clonotype_id))) +
  geom_bar(position = "fill") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Missing Clonotypes by Cell type", y = "Proportion")

# Save the visual 
ggsave("<<path to visual>>")
```
<img width="668" height="408" alt="upload_4de2e4b52c17797064130544bd90c19a" src="https://github.com/user-attachments/assets/43547e0a-b7f7-465f-b831-56571c9ff89b" />

```{r}
# Plot cell types vs. clonotype presence
bcells_only@meta.data %>%
  ggplot(aes(x = hash.ID, fill = is.na(raw_clonotype_id))) +
  geom_bar(position = "fill") +
  labs(title = "Missing Clonotypes by Hashtag", y = "Proportion")

# Save the visual 
ggsave("<<path to visual>>")
```
<img width="659" height="407" alt="upload_8a58d5cfedbe09ae3b519096fe48a7a7" src="https://github.com/user-attachments/assets/512df2d2-67f0-4154-9dda-a1295ac484c7" />

## Filtering poor quality cells
```{r}
# Filter Seurat object to remove NA clonotypes and low-quality cells
bcells_only_filtered <- subset(
  bcells_only,
  subset = 
    !is.na(raw_clonotype_id) &        # Remove cells without clonotypes
    is_cell == "true" &               # Keep only cells marked as valid
    nCount_RNA > 500 &                # Minimum RNA counts threshold
    nFeature_RNA > 200 &              # Minimum detected genes threshold
    percent.mt < 10                   # Maximum mitochondrial threshold
)

bcells_only_filtered  # 17030 features across 903 samples

# Verify filtering
cat("Original cells:", ncol(bcells_only), "\n",
    "Filtered cells:", ncol(bcells_only_filtered))

bcells_only_filtered@meta.data %>%
  ggplot(aes(x = hash.ID, fill = SingleR.labels)) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d(
    name = "Cell types",  # Legend title
    option = "A",     # Try "A" (magma), "B" (inferno), "C" (plasma), "D" (viridis)
    begin = 0.1,      # Adjust start of color range (0-1)
    end = 0.9         # Adjust end of color range (0-1)
  ) +
  labs(title = "Clonotypes in annotated cell types in each patient", 
       y = "Proportion", 
       x = "Patients") +
  theme_classic(base_family = "sans") +
  theme(
    # Universal text size (12pt)
    text = element_text(size = 14, family = "sans"),
    # Axis labels
    axis.text.x = element_text(angle = 0, hjust = 0.5,  size = 14),
    axis.text.y = element_text(hjust = 0, size = 14),  # Align the lables to the left
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    #Layout
    panel.grid = element_blank(),
    legend.position = "right",              # Legend on the right
    # Move y-axis labels to the right (mirror axis)
    axis.text.y.left = element_text(hjust = 1, margin = margin(l = 10)),
    # Legend
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
    
  ) 

# Save the visual 
ggsave("<<path to visual>>")
```
<img width="1172" height="652" alt="upload_28538a4b9fa598cd50be244f8e36813d" src="https://github.com/user-attachments/assets/34572c91-7066-4d57-81cc-f222a02ee77d" />

## Prepare Clonotype-Hashtag Data
```{r}
 # Extract metadata with proper barcode handling
meta_bcells <- as.data.frame(bcells_only_filtered@meta.data)
head(meta_bcells)
# Verify barcode matching
head(meta_bcells %>% select(Run, hash.ID, raw_clonotype_id))
```
<img width="554" height="198" alt="upload_cf2b3da30f0d5880a23a77e52712b55a" src="https://github.com/user-attachments/assets/a216eeb2-88cf-47fb-a05c-6e2b4efb2d14" />

```{r}
# Check for proper barcode-clonotype matching
meta_bcells %>% 
  group_by(hash.ID) %>% 
  summarise(
    cells_with_clonotypes = sum(!is.na(raw_clonotype_id)),
    total_cells = n()
  ) %>%
  mutate(clonotype_coverage = cells_with_clonotypes/total_cells)

```
<img width="637" height="196" alt="upload_7f3b903a0d1dc3e89a522faf36ad3465" src="https://github.com/user-attachments/assets/004e3313-376d-475e-8f36-101134c839a0" />

```{r}
# Clonotype Visualization by Hashtag 
# Top clonotypes per hashtag group
top_clonos <- meta_bcells %>%
  group_by(hash.ID) %>%
  count(raw_clonotype_id) %>%
  top_n(10, n) %>%
  ungroup()

head(top_clonos)
```
<img width="474" height="246" alt="upload_614a0ef9ce8b76c1a490d28e967df7be" src="https://github.com/user-attachments/assets/313439d2-163c-45cd-aa33-a4892fb852ae" />

## Number of clonotypes in each patient
```{r}
library(UpSetR)

upset_data <- meta_bcells %>%
  distinct(raw_clonotype_id, hash.ID) %>%
  table() %>%
  as.data.frame.matrix()

UpSetR::upset(upset_data, nsets = 10, order.by = "freq")
```
<img width="665" height="410" alt="upload_73b6812f409cd7403d4b8457e89e03aa" src="https://github.com/user-attachments/assets/10fac397-8026-4ef9-98cd-a97cdfe8e73a" />

## Clonotype sizes
```{r}
# Step 1: Count cells per clonotype per sample
clone_sizes <- meta_bcells %>%
  group_by(hash.ID, raw_clonotype_id) %>%
  summarise(clone_size = n(), .groups = "drop")

# Step 2: Count how many clones of each size per sample
clone_distribution <- clone_sizes %>%
  group_by(hash.ID, clone_size) %>%
  summarise(count = n(), .groups = "drop")

# Optional: limit to expanded clones only (e.g. clone_size > 1)
clone_distribution <- clone_distribution %>%
  filter(clone_size > 1)

# Step 3: Compute percentage if desired
clone_distribution <- clone_distribution %>%
  group_by(hash.ID) %>%
  mutate(percentage = count / sum(count))

# Step 4: Plot
ggplot(clone_distribution, 
       aes(x = hash.ID, y = percentage, fill = factor(clone_size))) +
  geom_bar(stat = "identity") +
  labs(x = "Samples",
       y = "Percentage of Expanded Clones",
       fill = "Clone Size",
       title = "Expanded Clone Size Distribution per Sample") +
  theme_classic(base_family = "sans") +
  theme(
    # Universal text size (12pt)
    text = element_text(size = 12, family = "sans"),
    # Axis labels
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
    axis.text.y = element_text(hjust = 0, size = 12),  # Align V genes to the left
    axis.title.y = element_blank(),         # Remove y-axis title
    plot.title = element_text(hjust = 0.5, face = "bold"),
    #Layout
    panel.grid = element_blank(),
    legend.position = "right",              # Legend on the right
    # Move y-axis labels to the right (mirror axis)
    axis.text.y.left = element_text(hjust = 1, margin = margin(l = 10)),
    # Legend
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
    )
```
<img width="661" height="413" alt="upload_bbd21bb786d7cfe9418f7adcbfb8e309" src="https://github.com/user-attachments/assets/5290d4e7-66be-4eda-9145-fc6bb0987c9c" />

## V-J gene Pairing
The chord diagram visualizes the most frequent V-J gene pairings in the 
immune repertoire, specifically focusing on the top 10% of pairings by 
occurrence. The circular plot represents V genes and J genes as segments 
around the circumference, with connecting ribbons indicating pairing events.
The thickness of each ribbon corresponds to the relative frequency of that 
particular V-J combination, revealing which pairings dominate the immune 
repertoire. The color-coding system assigns distinct colors to each gene 
segment, making it easy to track individual V and J genes across multiple 
pairings. The directional arrows emphasize the biological flow from V to J
genes during recombination. By filtering to the top 10% most frequent pairs, 
the plot highlights preferential or recurrent V-J combinations that may
indicate biological constraints, antigen-driven selection, or technical 
biases in the sequencing approach. The minimalist design focuses attention 
on these dominant pairings while the legend provides essential decoding of
the gene segment colors. This visualization is particularly valuable for
identifying public clonotypes (shared across individuals) or detecting 
anomalous pairing patterns that might suggest sequencing artifacts or 
immune dysregulation. The circular layout effectively demonstrates both 
the diversity of V-J combinations and the presence of preferred pairings 
in the dataset.


```{r}
pdf("~/singlecell/new_joseph/bcells_only_vj_gene_chord_diagram.pdf", width = 14, height = 10)  # Adjust width and height as needed

library(circlize)
library(RColorBrewer)

# Prepare data and create color palettes
vj_pairs <- meta_bcells %>%
  count(v_gene, j_gene) %>%
  filter(n > quantile(n, 0.9))

v_genes <- unique(vj_pairs$v_gene)
j_genes <- unique(vj_pairs$j_gene)

# Color palettes
v_colors <- colorRampPalette(brewer.pal(9, "Reds"))(length(v_genes))
names(v_colors) <- v_genes
j_colors <- colorRampPalette(brewer.pal(9, "Blues"))(length(j_genes))
names(j_colors) <- j_genes
gene_colors <- c(v_colors, j_colors)

# Initialize plot
circos.par(gap.after = c(rep(3, length(v_genes)-1), 15, rep(3, length(j_genes)-1), 15))

# Create chord diagram
chordDiagram(vj_pairs,
             grid.col = gene_colors,
             transparency = 0.3,
             annotationTrack = c("grid", "axis"),
             directional = 1,
             direction.type = "arrows",
             link.arr.type = "big.arrow")

# Add gene type highlights
highlight.sector(v_genes, track.index = 1, col = "#FF000010", 
                 text = "V GENES", cex = 0.9, text.col = "darkred")
highlight.sector(j_genes, track.index = 1, col = "#0000FF10", 
                 text = "J GENES", cex = 0.9, text.col = "darkblue")



# Create unified legend in single column
legend("right",
       legend = c("GENE TYPES", "V genes", "J genes", "", "SPECIFIC GENES", names(gene_colors)),
       fill = c(NA, "indianred", "steelblue", NA, NA, gene_colors),
       border = c(NA, "indianred", "steelblue", NA, NA, rep(NA, length(gene_colors))),
       title.adj = 0,
       
       ncol = 1,
       bty = "n",
       cex = 0.8,
       x.intersp = 0.5,
       text.width = max(strwidth(c("V genes", "J genes", names(gene_colors))) * 0.8,
       title = NULL))  # Remove automatic title

title("V-J Gene Pairing (Top 10%)", cex.main = 1.3)
circos.clear()

# Close PDF device
dev.off()
```
<img width="1378" height="1052" alt="upload_4a6f17fee00615ca895bc68b86301ba6" src="https://github.com/user-attachments/assets/84c7b0a1-c81b-4618-b36a-efc1d35b0e95" />

The chord diagram titled **"V–J Gene Pairing (Top 10%)"** visualizes the 
most frequent pairings between variable (V) and joining (J) genes in the 
immune repertoire, focusing on the top 10% most abundant combinations. 
V genes (shaded in red tones) and J genes (shaded in blue tones) are 
arranged around the circle, with directional arrows indicating pairings. 
Thicker and more opaque links represent higher frequencies of specific V–J 
combinations. The diagram highlights prominent interactions such as between 
IGHV4-39 and IGLV2-14 or IGKV1-5 and IGHV3-23, suggesting potential biases 
or preferences in V–J recombination. This pattern may reflect underlying 
biological mechanisms such as clonal expansion, antigen-driven selection, 
or germline gene usage. The unified legend distinguishes gene types and 
specific gene names, aiding interpretation.

## V-D-J Recombination Frequency (Bar Plot)
```{r}
# Prepare recombination data
recomb_data <- meta_bcells %>%
  mutate(vdj_combo = paste(v_gene, d_gene, j_gene, sep = "-")) %>%
  count(vdj_combo, sort = TRUE) %>%
  head(20)  # Top 20 combinations

# Plot recombination frequencies
ggplot(recomb_data, aes(x = reorder(vdj_combo, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 V-D-J Recombination Events",
       x = "V-D-J Combination", y = "Count") +
  theme_minimal()
```
<img width="671" height="409" alt="upload_4166d19b7679ebf631d501f0a8d89963" src="https://github.com/user-attachments/assets/4df1dfa4-491d-489c-8cdb-d5068009b51a" />


## CDR3 Length Distribution
```{r}
# Calculate CDR3 lengths (assuming cdr3_aa column exists)
meta_bcells <- meta_bcells %>%
  mutate(cdr3_length = nchar(cdr3))  # Replace cdr3_aa with your column

ggplot(meta_bcells, aes(x = cdr3_length, fill = SingleR.labels)) +
  geom_density(alpha = 0.5) +
  labs(title = "CDR3 Length Distribution",
       x = "CDR3 Length (aa)", y = "Density") +
  theme_classic2() +
  scale_fill_brewer(palette = "Set2", name = "Sample")
```
<img width="665" height="413" alt="upload_7cd5408950e98eafec76a77c4c919635" src="https://github.com/user-attachments/assets/756b9cf9-8695-4ad7-bf38-84abc7053dd2" />

```{r}
ggplot(meta_bcells, aes(x = cdr3_length, fill = SingleR.labels)) +
  geom_histogram(binwidth = 1, alpha = 0.5, position = "identity") +
  scale_x_continuous(breaks = seq(0, max(meta_bcells$cdr3_length), by = 2)) +  # Custom x-axis
  labs(title = "CDR3 Length Histogram (Absolute Counts)") +
  theme_classic2()
```
<img width="666" height="410" alt="upload_4111a1550dc77f592c312943662bbafb" src="https://github.com/user-attachments/assets/65cb4c34-ef66-49ca-a1e0-753a696e53a8" />

```{r}
library(tidyverse)

# Convert to data.table
meta_dt <- as.data.table(meta_bcells)

# Count number of contigs per row (cell)
contig_counts <- meta_dt %>%
  mutate(num_contigs = lengths(strsplit(contig_id, ";\\s*"))) %>%
  count(num_contigs, name = "num_cells")

# Print the result
contig_counts
```
<img width="591" height="202" alt="upload_e4ecba96beb1fc30d906d9eec61e916b" src="https://github.com/user-attachments/assets/45bb753c-a648-4458-938e-000d21242a14" />

```{r}
# Calculate CDR3 lengths
vdj_long <- vdj_long %>%
  filter(chain == "IGH") %>%
  mutate(cdr3_length = nchar(cdr3))  # Replace cdr3_aa with your column

# Plot length distribution
ggplot(vdj_long, aes(x = cdr3_length, fill = hash.ID)) +
  geom_density(alpha = 0.5) +
  labs(title = "CDR3 Length Distribution",
       x = "CDR3 Length (aa)", y = "Density") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2", name = "Sample")
```
<img width="664" height="410" alt="upload_3c055247315e4ac5210414a5a15a4a9a" src="https://github.com/user-attachments/assets/f9c21902-36b9-4bce-b504-f386e28411e7" />

```{r}
# Calculating the patients
ggplot(vdj_long, aes(x = cdr3_length, fill = hash.ID)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  scale_x_continuous(breaks = seq(0, max(vdj_long$cdr3_length), by = 2)) +  # Custom x-axis
  scale_fill_viridis_d(
    name = "Sample",  # Legend title
    option = "D",     # Try "A" (magma), "B" (inferno), "C" (plasma), "D" (viridis)
    begin = 0.1,      # Adjust start of color range (0-1)
    end = 0.9         # Adjust end of color range (0-1)
  ) +
  labs(title = "CDR3 Length Histogram (Absolute Counts)") +
  theme_classic(base_family = "sans") +
      xlab("CDR3 Length") +
  theme(
    # Universal text size (12pt)
    text = element_text(size = 12, family = "sans"),
    # Axis labels
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14),
    axis.text.y = element_text(hjust = 0, size = 12),  # Align V genes to the left
    axis.title.y = element_blank(),         # Remove y-axis title
    plot.title = element_text(hjust = 0.5, face = "bold"),
    #Layout
    panel.grid = element_blank(),
    legend.position = c(0.95, 0.95),              # Legend on the right
    # Move y-axis labels to the right (mirror axis)
    axis.text.y.left = element_text(hjust = 1, margin = margin(l = 10)),
    # Legend
    legend.justification = c(1, 1),   # Anchors the legend at its top-right
    legend.box.margin = margin(6, 6, 6, 6), # Add some margin inside the box
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
```
<img width="667" height="410" alt="upload_882530ff7ecf7384f36715b35d5fa10e" src="https://github.com/user-attachments/assets/b22e17c2-8064-4633-b1a0-3901076a8f6a" />

```{r}
library(ggseqlogo)

# Extract and clean CDR3 amino acid data
cdr3_data <- vdj_long %>%
  filter(!is.na(cdr3), nchar(cdr3) == 11) %>%
  select(hash.ID, cdr3) %>%
  group_by(hash.ID) %>%
  filter(n() >= 5)  # Keep only samples with ≥ 5 sequences

# Split sequences into a named list by hash.ID
cdr3_split <- split(cdr3_data$cdr3, cdr3_data$hash.ID)

# Create the plot with increased font sizes
ggseqlogo(cdr3_split, method = "prob", facet = "wrap") +
  labs(title = "Motif plot of the heavy chain of the v-gene") +
  theme_classic2() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
    strip.text = element_text(size = 14),  # Facet labels
    axis.title = element_text(size = 14),  # Axis titles
    axis.text = element_text(size = 12),   # Axis text (sequence positions)
    legend.title = element_text( size = 14), # Legend title
    legend.text = element_text(size = 12)  # Legend items
  )

# Save with high resolution
ggsave("~/singlecell/Figure5c_Facet_motif_new.png", 
       width = 16, 
       height = 10, 
       dpi = 800)
```
<img width="665" height="399" alt="upload_4810f3f010de58bf59263f18cffd5a43" src="https://github.com/user-attachments/assets/b9e8be99-25f8-4882-9c65-918702389a73" />

```{r}
# Calculate V gene usage frequencies
v_usage <- unique_data %>%
  group_by(hash.ID, v_gene) %>%  # Replace v_gene with your actual V gene column name
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

# Heatmap of V gene usage
ggplot(v_usage, aes(x = hash.ID, y = v_gene, fill = freq)) +
  geom_tile(color = "white", linewidth = 0.2) +
  # Viridis color scale (option: try "magma", "plasma", or "inferno")
  scale_fill_viridis(
    name = "Frequency",
    option = "viridis",
    direction = 1  # Reversed for darker = higher
  ) +
  labs(
    title = "V-gene usage",
    x = NULL,
    y = NULL  # Remove y-axis title (labels on right)
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    # Universal text size (12pt)
    text = element_text(size = 12, family = "sans"),
    # Axis labels
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
    axis.text.y = element_text(hjust = 0, size = 12),  # Align V genes to the left
    axis.title.y = element_blank(),         # Remove y-axis title
    plot.title = element_text(hjust = 0.5, face = "bold"),
    #Layout
    panel.grid = element_blank(),
    legend.position = "right",              # Legend on the right
    # Move y-axis labels to the right (mirror axis)
    axis.text.y.left = element_text(hjust = 1, margin = margin(l = 10)),
    # Legend
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
    
  ) 
```
<img width="692" height="1044" alt="upload_daf92bf441fc60c8293315c32112906f" src="https://github.com/user-attachments/assets/2efd2066-87dd-4960-9baa-be0e31ef5396" />

```{r}
# Calculate J gene usage frequencies
j_usage <- unique_data %>%
  group_by(hash.ID, j_gene) %>%  # Replace v_gene with your actual V gene column name
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

# Heatmap of J gene usage
p2 <- ggplot(j_usage, aes(x = hash.ID, y = j_gene, fill = freq)) +
  geom_tile(color = "white", linewidth = 0.2) +
  # Viridis color scale (option: try "magma", "plasma", or "inferno")
  scale_fill_viridis(
    name = "Frequency",
    option = "viridis",
    direction = 1  # Reversed for darker = higher
  ) +
  labs(
    title = "J-gene usage",
    x = NULL,
    y = NULL  # Remove y-axis title (labels on right)
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    # Universal text size (12pt)
    text = element_text(size = 12, family = "sans"),
    # Axis labels
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
    axis.text.y = element_text(hjust = 0, size = 12),  # Align V genes to the left
    axis.title.y = element_blank(),         # Remove y-axis title
    plot.title = element_text(hjust = 0.5, face = "bold"),
    #Layout
    panel.grid = element_blank(),
    legend.position = "right",              # Legend on the right
    # Move y-axis labels to the right (mirror axis)
    axis.text.y.left = element_text(hjust = 1, margin = margin(l = 10)),
    # Legend
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
    
  ) 
```
<img width="601" height="384" alt="upload_f2e43c24ea83c0a1bb26bf493c0803fc" src="https://github.com/user-attachments/assets/298fec5b-a5b6-4324-baa2-bed5b73934d7" />

```{r}
# V & J
v_j_usage <- unique_data %>%
  group_by(v_gene, j_gene) %>%  # Replace v_gene with your actual V gene column name
  filter(grepl("^IGHJ", j_gene) & grepl("^IGHV", v_gene)) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

# Heatmap of V & J gene usage
ggplot(v_j_usage, aes(x = j_gene, y = v_gene, fill = freq)) +
  geom_tile(color = "white", linewidth = 0.2) +
  # Viridis color scale (option: try "magma", "plasma", or "inferno")
  scale_fill_viridis(
    name = "Frequency",
    option = "viridis",
    direction = 1  # Reversed for darker = higher
  ) +
  labs(
    title = "General IGH VJ pair",
    x = NULL,
    y = NULL  # Remove y-axis title (labels on right)
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    # Universal text size (12pt)
    text = element_text(size = 12, family = "sans"),
    # Axis labels
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
    axis.text.y = element_text(hjust = 0, size = 12),  # Align V genes to the left
    axis.title.y = element_blank(),         # Remove y-axis title
    plot.title = element_text(hjust = 0.5, face = "bold"),
    #Layout
    panel.grid = element_blank(),
    legend.position = "right",              # Legend on the right
    # Move y-axis labels to the right (mirror axis)
    axis.text.y.left = element_text(hjust = 1, margin = margin(l = 10)),
    # Legend
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
    
  ) 

write_csv(unique_data, "~/singlecell/unique_final_dataset.csv")

```
<img width="582" height="675" alt="upload_7d21042746db4557dff7e41e3789313e" src="https://github.com/user-attachments/assets/8634d3e1-596f-44e8-b8d3-18aa18b72f31" />

```{r}
library(ggalluvial)

# Filter and simplify for visualization (e.g., top 10 V/J genes)
#top_v <- meta_bcells %>% count(v_gene, sort = TRUE) %>% slice_head(n = 10) %>% pull(v_gene)
#top_j <- meta_bcells %>% count(j_gene, sort = TRUE) %>% slice_head(n = 10) %>% pull(j_gene)

top_v <- names(head(sort(table(unique_data$v_gene), decreasing = TRUE), 10))
top_j <- names(head(sort(table(unique_data$j_gene), decreasing = TRUE), 10))

top_vdj <- unique_data

# Filter to top V and J genes only
filtered_df <- top_vdj[top_vdj$v_gene %in% top_v & top_vdj$j_gene %in% top_j, ]

# Replace NAs (if any) in d_gene
filtered_df$d_gene[is.na(filtered_df$d_gene)] <- "Unassigned"

# Count combinations
plot_data <- as.data.frame(table(filtered_df$v_gene,
                                 filtered_df$d_gene,
                                 filtered_df$j_gene,
                                 filtered_df$hash.ID))

# Clean column names
colnames(plot_data) <- c("v_gene", "d_gene", "j_gene", "hash.ID", "n")

# Optionally convert to factors
plot_data$v_gene <- factor(plot_data$v_gene)
plot_data$d_gene <- factor(plot_data$d_gene)
plot_data$j_gene <- factor(plot_data$j_gene)
plot_data$hash.ID <- factor(plot_data$hash.ID)

# Drop zero-count rows
plot_data <- plot_data[plot_data$n > 0, ]

# Ploting

ggplot(plot_data,
       aes(axis1 = v_gene, axis2 = d_gene, axis3 = j_gene, axis4 = hash.ID, y = n)) +
  geom_alluvium(aes(fill = v_gene), width = 3/12) +
  geom_stratum(width = 3/12, fill = "lightgrey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("V", "D", "J", "Sample"), expand = c(.05, .05)) +
  labs(title = "V-J Segment Usage Across Samples",
       y = "Clone Count") +
  theme_minimal(base_family = "sans") +
  theme(
    # Universal text size (12pt)
    text = element_text(size = 12, family = "sans"),
    # Axis labels
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
    axis.text.y = element_text(hjust = 0, size = 12),  # Align V genes to the left
    plot.title = element_text(hjust = 0.5, face = "bold"),
    #Layout
    panel.grid = element_blank(),
    legend.position = "right",              # Legend on the right
    # Move y-axis labels to the right (mirror axis)
    axis.text.y.left = element_text(hjust = 1, margin = margin(l = 10)),
    # Legend
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
) 
```
<img width="976" height="626" alt="upload_e4cdafff3b8cc1dac934133a87921758" src="https://github.com/user-attachments/assets/dbafc9fa-3c1c-4180-a54a-df4540d477de" />




