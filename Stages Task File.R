#Stage One Task

#Task 1 - GC% Calculation

#Initial GC% calculator function
GC_percent <- function(input_gene) {{
  #split gene
  input_gene <- strsplit(x = input_gene, split = "")[[1]]
  print(length(input_gene))
  
  #core operations
  gc_counter <- 0
  for (nuc in input_gene) {
    
    #print(nuc)
    if(nuc == 'G' | nuc == 'C') { #if statement
      gc_counter = gc_counter+1}
  }
}
  (gc_counter / length(input_gene)) * 100
}


#Making it robust enough to handle nuc sequences written in upper and lower case
GC_percent <- function(input_gene) {{
  #normalize nucleotide case (converts everything to UPPERCASE)
  input_gene <- toupper(input_gene)
  
  #Split the string into indidual letters
  input_gene <- strsplit(x = input_gene, split = "")[[1]]
  
  #core operations
  gc_counter <- 0
  for (nuc in input_gene) {
    if(nuc == 'G' | nuc == 'C') {
      gc_counter = gc_counter+1}
  }
}
  #calculate the percentage
  (gc_counter / length(input_gene)) * 100
}

#Task 1 Example - Executing with given sequence
input_gene <- 'gcaTTTAT'
GC_percent(input_gene)



#Task 2 - Protein weights in KiloDalton

AminoA_weights <- c( A = 89.09, R = 174.20, N = 132.12, D = 133.10, C = 121.15, E = 147.13, Q = 146.15, G = 75.07, H = 155.16, I = 131.18, L = 131.18, K = 146.19, M = 149.21, F = 165.19, P = 115.13, S = 105.09, T = 119.12, W = 204.23, Y = 181.19, V = 117.15 )

Protein_MW <- function(protein_seq = "ESTHER") {
  protein_seq <- toupper(protein_seq)
  protein_seq <- strsplit(protein_seq, split = "")[[1]]
  
  #core operations
  total_mass <- 0
  for (AA in protein_seq) {
    if (AA %in% names(AminoA_weights)) {
      total_mass <- total_mass + AminoA_weights[AA]
    } else {
      total_mass <- total_mass + 0 #unknown amino acid will count as 0
    }
  }
  #convert to KiloDalton
  total_mass / 1000
}

#Task 2 question
result <- Protein_MW()
print(result)



#Stage 2 Task
#Installing Packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("igraph")
install.packages("dplyr")
install.packages("tidyr")
install.packages("gridExtra")
#Loading Packages
library(readxl)
library(ggplot2)
library(pheatmap)
library(igraph)
library(dplyr)
library(tidyr)
library(gridExtra)

#importing dataframes

#HBR vs UHR
fileOnline <- 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv'
#read in the file
read.csv(fileOnline, header = T)
#assign
HBR.vs.UHR <- read.csv(fileOnline, header = T)

#Deg_CHR22
fileOnline <- 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv'
#read in the file
read.csv(fileOnline, header = T)
#assign
Deg_CHR22 <- read.csv(fileOnline, header = T)

#Breast_CancerDS
fileOnline <- 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv'
#read in the file
read.csv(fileOnline, header = T)
#assign
Breast_CancerDS <- read.csv(fileOnline, header = T)


#Stage 2 Task 

#Part 1 -Gene Expression
#(a) Heatmap
pheatmap::pheatmap(mat = HBR.vs.UHR[ , 2:7],
                   border_color = 'black',
                   legend = T,
                   labels_row = HBR.vs.UHR$X,
                   fontsize_row = 5,
                   cluster_rows = T,
                   cluster_cols = T,
                   color = colorRampPalette(c('white', 'lightblue', 'navyblue'))(100)
                 
)

#(b) Volcano plot
#
ggplot(Deg_CHR22,
       aes(x = log2FoldChange, y = X.log10PAdj, color = significance)) +
  geom_point(size = 1.8, alpha = 1) +
  scale_color_manual(
    values = c("up" = "green",
      "down" = "orange",
      "ns" = "grey")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(
    x = "log2 Fold Change",
    y = "X.log10PAdj",
    color = "Significance"
  ) +
  theme_minimal(base_size = 12)


#Part 2 - Breast Cancer Data Exploration
#(c) Scatter plot - radius vs texture
#colour grade
col_vec <- c("B" = "blue", # Benign
             "M" = "red"  # Malignant
             )

plot(x = Breast_CancerDS$radius_mean,
     y = Breast_CancerDS$texture_mean,
     xlab = 'radius_mean',
     ylab = 'texture_mean',
     xlim = c(7, 28),
     ylim = c(9, 40),
     las = 1,
     main = 'Radius vs Texture',
     col = col_vec[Breast_CancerDS$diagnosis],
     pch = 19,
     cex = 0.5)


#(d) Correlation Heatmap
key_features <- Breast_CancerDS[, c("radius_mean", "texture_mean", 
                                    "perimeter_mean","area_mean", 
                                    "smoothness_mean", "compactness_mean")]

# Compute correlation matrix
corr_matrix <- cor(key_features, use = "complete.obs")

# Plot with pheatmap
pheatmap::pheatmap( mat = corr_matrix,
                    display_numbers = TRUE, # annotate correlation values
                    number_format = "%.1f", # round to 1 decimal place
                    fontsize_number = 11,
                    border_color = 'black',
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    color = colorRampPalette(c("white", "lightblue", "navyblue"))(100),
                    main = "Correlation Heatmap of Key Features")


#(e) Scatter plot - Smoothness vs Compactness
plot(x = Breast_CancerDS$smoothness_mean,
     y = Breast_CancerDS$compactness_mean,
     xlab = 'smoothness_mean',
     ylab = 'compactness_mean',
     xlim = c(0.05, 0.17),
     ylim = c(0.01, 0.37),
     las = 1,
     main = 'Smoothness vs Compactness',
     col = as.factor(Breast_CancerDS$diagnosis),
     pch = 19,
     cex = 0.7)
grid(col = "lightgray", lty = "dotted", lwd = 1)


#(f) Density plot
# Split the data
area_Mal <- Breast_CancerDS$area_mean[Breast_CancerDS$diagnosis == "M"]
area_Be <- Breast_CancerDS$area_mean[Breast_CancerDS$diagnosis == "B"]

# Compute densities
Dens_M <- density(area_Mal)
Dens_B <- density(area_Be)

# Plot the first density
plot(Dens_M,
     col = "red",
     lwd = 2,
     xlab = "Area Mean",
     ylab = "Density",
     main = "Kernel Density Estimates of Area Mean",
     ylim = c(0, 0.003)
)

# Add the second density
lines(Dens_B, col = "blue", lwd = 2)

# Add legend
legend("topright",
       legend = c("Malignant (M)", "Benign (B)"),
       col = c("red", "blue"),
       lwd = 1,
       title = "Diagnosis"
)


#part 3
#Task 0: Orientation and data hygiene (mandatory, graded)
# 1. Install and load packages
# install.packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("pheatmap") 
install.packages("igraph") 
install.packages("tidyr") 
install.packages("dplyr")
install.packages("ggrepel")
library(readxl)
library(ggplot2)
library(pheatmap)
library(igraph)
library(tidyr)
library(dplyr)

# 2. Define the Helper Function
transparent_color <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  return(t.col)
}

# 3. Initialize the Color Palette
hb_pal <- c("#4e79a7", "#8cd17d", "#e15759", "#fabfd2", "#a0cbe8", 
            "#59a14f", "#b07aa1", "#ff9d9a", "#f28e2b", "#f1ce63",
            "#79706e", "#d4a6c8", "#e9e9e9", "#ffbe7d", "#bab0ac",
            "#9d7660", "#d37295", "#86bcb6", "#362a39", "#cd9942")

# 4. Test the Palette
plot(1:length(hb_pal), 1:length(hb_pal), col = hb_pal, pch = 19, 
     main = "Project Color Palette", xlab = "Index", ylab = "Color")

#Task 1. Reproduce panel 2a: Cell-type ratio distributions
# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)

# 1. Setup Data Hygiene & Color Palette
file_path <- "C:\\Users\\user\\Downloads\\hb_stage_2.xlsx"

hb_pal <- c("#4e79a7", "#8cd17d", "#e15759", "#fabfd2", "#a0cbe8", 
            "#59a14f", "#b07aa1", "#ff9d9a", "#f28e2b", "#f1ce63",
            "#79706e", "#d4a6c8", "#e9e9e9", "#ffbe7d", "#bab0ac",
            "#9d7660", "#d37295", "#86bcb6", "#362a39", "#cd9942")

# 2. Read Sheet 'a' from the Excel file
# If you are running this locally, ensure the file is in your working directory.
df_a <- read_excel(file_path, sheet = "a")

# 3. Create the Boxplot
# We use reorder() to sort the x-axis by the median new_ratio (descending)
p2a <- ggplot(df_a, aes(x = reorder(cell_type, -new_ratio, FUN = median), 
                        y = new_ratio, 
                        fill = cell_type)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.6) +
  scale_fill_manual(values = hb_pal) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), # Rotated labels
    legend.position = "none",                                    # Hide legend for clarity
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Figure 2a: Cell-type ratio distributions",
    x = "Cell Type",
    y = "New Ratio (mRNA labeling)"
  )

# Display the plot
print(p2a)

# Optional: Save the plot
# ggsave("panel_2a_reproduction.pdf", plot = p2a, width = 7, height = 5)


#Task 2. Reproduce panel 2b: Half-life vs alpha-life scatter
# Load libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel) # For better label placement

# 1. Read the data from Sheet B
file_path <- "C:\\Users\\user\\Downloads\\hb_stage_2.xlsx"
df_b <- read_excel(file_path, sheet = "b")

# 2. Log-transformation and Thresholding
# We calculate thresholds based on medians for the quadrants
# (You can adjust these values if your paper uses specific cutoffs)
thresh_alpha <- median(df_b$alpha, na.rm = TRUE)
thresh_hl <- median(df_b$half_life, na.rm = TRUE)

df_b <- df_b %>%
  mutate(
    log_alpha = log2(alpha),
    log_hl = log2(half_life),
    regime = case_when(
      alpha > thresh_alpha & half_life > thresh_hl ~ "High/Stable",
      alpha > thresh_alpha & half_life <= thresh_hl ~ "High/Unstable",
      alpha <= thresh_alpha & half_life > thresh_hl ~ "Low/Stable",
      alpha <= thresh_alpha & half_life <= thresh_hl ~ "Low/Unstable"
    )
  )

# 3. Create the Scatter Plot
p2b <- ggplot(df_b, aes(x = log_alpha, y = log_hl, color = regime)) +
  geom_point(alpha = 0.4, size = 1) +
  # Add Cutoff lines
  geom_vline(xintercept = log2(thresh_alpha), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = log2(thresh_hl), linetype = "dashed", color = "grey40") +
  # Label Exemplar Genes (Camp, Ccr2)
  geom_label_repel(data = filter(df_b, cell %in% c("Camp", "Ccr2")),
                   aes(label = cell),
                   color = "black", fontface = "bold",
                   box.padding = 0.5, point.padding = 0.5) +
  # Custom Colors and Styling
  scale_color_manual(values = c("#e15759", "#4e79a7", "#f28e2b", "#76b7b2")) +
  theme_classic() +
  labs(
    title = "Figure 2b: Global Kinetic Regimes",
    x = "Transcription Rate [log2(alpha)]",
    y = "mRNA Half-life [log2(half-life)]",
    color = "Kinetic Regime"
  )

# Display plot
print(p2b)

# Save plot
# ggsave("panel_2b_reproduction.png", plot = p2b, width = 7, height = 6)

#Task 3. Reproduce panel 2c: Heatmap across cell types and time
# Load libraries
library(readxl)
library(pheatmap)
library(dplyr)
library(tidyr)

# 1. Read the data from Sheet C
file_path <- "C:\\Users\\user\\Downloads\\hb_stage_2.xlsx"
df_c <- read_excel(file_path, sheet = "c")

# 2. Convert to matrix
# Convert the 'genes' column to row names
mat_data <- as.matrix(df_c[,-1])
rownames(mat_data) <- df_c$genes

# 3. Create Column Annotations
# Column names look like "Macrophagen00h", "Bn12h", etc.
col_names <- colnames(mat_data)

# Extract Time (last 4 characters, e.g., n00h)
times <- gsub(".*(n[0-9]{2}h)$", "\\1", col_names)

# Extract Cell Type (everything before the time)
cell_types <- gsub("n[0-9]{2}h$", "", col_names)

# Create the annotation data frame
ann_col <- data.frame(
  CellType = cell_types,
  Time = times
)
rownames(ann_col) <- col_names

# 4. Define Annotation Colors (optional but matches Panel 2a)
# Mapping specific cell types to colors from our palette
ann_colors <- list(
  CellType = c(
    Macrophage = "#59a14f", 
    Monocyte = "#a0cbe8", 
    Neutrophil = "#8cd17d", 
    NK = "#b07aa1", 
    pDC = "#ff9d9a", 
    B = "#4e79a7", 
    T = "#e15759"
  )
)

# 5. Build the Heatmap
pheatmap(mat_data,
         cluster_cols = FALSE,           # DO NOT cluster columns 
         cluster_rows = TRUE,            # Cluster genes with similar patterns
         show_rownames = FALSE,          # Hide gene names if the list is long
         scale = "row",                  # Z-score normalization per gene
         annotation_col = ann_col,       # Add the Type/Time bars at the top
         annotation_colors = ann_colors, # Apply the project palette
         main = "Figure 2c: Temporal Gene Expression",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = NA)

#Task 4. Reproduce panel 2d: Pathway enrichment heatmap
# Load libraries
library(readxl)
library(pheatmap)
library(RColorBrewer)

# 1. Read the data from Sheet d_1
file_path <- "C:\\Users\\user\\Downloads\\hb_stage_2.xlsx"
df_d <- read_excel(file_path, sheet = "d_1")

# 2. Prepare the matrix
# Set pathway names as the index and remove the character column
mat_d <- as.matrix(df_d[,-1])
rownames(mat_d) <- df_d$pathway

# 3. Define a Diverging Color Palette
# Center the scale at 0 (white) with blue for negative and red for positive
my_palette <- colorRampPalette(c("#4e79a7", "white", "#e15759"))(100)

# 4. Generate the Heatmap
pheatmap(mat_d,
         cluster_cols = FALSE,          # DO NOT cluster timepoints
         cluster_rows = FALSE,          # DO NOT cluster pathways (maintain manual order)
         color = my_palette,
         display_numbers = FALSE,       # Keep it clean
         border_color = "white",        # Subtle grid
         main = "Figure 2d: Pathway Activity Over Time",
         angle_col = 45,                # Rotate time labels
         legend_labels = "Activity Score")

#Task 5. Reproduce panel 2e: Bubble plot of kinetic regimes
# Load libraries
library(readxl)
library(ggplot2)

# 1. Read the data from Sheet 'e'
file_path <- "C:\\Users\\user\\Downloads\\hb_stage_2.xlsx"
df_e <- read_excel(file_path, sheet = "e")

# 2. Create the Bubble Plot
p2e <- ggplot(df_e, aes(x = half_life, y = alpha, size = count, color = stage)) +
  geom_point(alpha = 0.7) +
  # Use log scales for better distribution (common in kinetic plots)
  scale_x_log10() + 
  scale_y_log10() +
  # Apply project palette for stages
  scale_color_manual(values = c("00h" = "#79706e", "02h" = "#8cd17d", 
                                "06h" = "#e15759", "72h" = "#4e79a7")) +
  # Customize size range for readability
  scale_size_continuous(range = c(2, 10)) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Figure 2e: Functional Kinetic Regimes",
    x = "mRNA Half-life (hours, log10)",
    y = "Transcription Rate (alpha, log10)",
    size = "Gene Count",
    color = "Stage"
  )

# Display plot
print(p2e)

# Optional: label the most prominent points
# library(ggrepel)
# p2e + geom_text_repel(aes(label = Description), size = 3, max.overlaps = 5)

#Task 6. Reproduce panel 2f: Stacked proportions
# Load libraries
library(readxl)
library(ggplot2)
library(dplyr)

# 1. Read the data from Sheet 'f'
file_path <- "C:\\Users\\user\\Downloads\\hb_stage_2.xlsx"
df_f <- read_excel(file_path, sheet = "f")

# 2. Subset the data
# Filter for start and end timepoints
df_f_sub <- df_f %>% 
  filter(stage %in% c("s00h", "s72h"))

# 3. Create the Stacked Barplot
p2f <- ggplot(df_f_sub, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  # Fixed y-axis limit as per requirements
  scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0)) +
  # Apply colors from the project palette (B and Plasma)
  scale_fill_manual(values = c("B" = "#4e79a7", "Plasma" = "#b07aa1")) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "Figure 2f: B-lineage Composition Shift",
    x = "Timepoint",
    y = "Proportion of Total Cells",
    fill = "Cell Type"
  )

# Display plot
print(p2f)

# Save plot
# ggsave("panel_2f_reproduction.png", plot = p2f, width = 5, height = 6)

#Task 7. Reproduce panel 2g: Directed cell–cell interaction network
# 1. Define the Helper Function (Must be run once per session)
transparent_color <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}

# 2. Define the Color Palette
hb_pal <- c("#4e79a7", "#8cd17d", "#e15759", "#fabfd2", "#a0cbe8", 
            "#59a14f", "#b07aa1", "#ff9d9a", "#f28e2b", "#f1ce63",
            "#79706e", "#d4a6c8", "#e9e9e9", "#ffbe7d", "#bab0ac",
            "#9d7660", "#d37295", "#86bcb6", "#362a39", "#cd9942")

# 3. Load libraries
library(readxl)
library(igraph)

# 4. Read the data from Sheet 'g'
file_path <- "C:\\Users\\user\\Downloads\\hb_stage_2.xlsx"
df_g <- read_excel(file_path, sheet = "g")

# 5. Convert to Adjacency Matrix
mat_g <- as.matrix(df_g[,-1])
rownames(mat_g) <- df_g[[1]]
colnames(mat_g) <- df_g[[1]]

# 6. Build Directed Graph
net <- graph_from_adjacency_matrix(mat_g, mode = "directed", weighted = TRUE, diag = FALSE)

# 7. Remove Zero-weight Edges
net <- delete.edges(net, which(E(net)$weight <= 0))

# 8. Define Layout and Styling
set.seed(42) 
layout_fr <- layout_with_fr(net)

# Map edge width and arrow size to the weight
E(net)$width <- E(net)$weight * 5
E(net)$arrow.size <- E(net)$weight * 1.5

# Map node colors from the palette
node_colors <- hb_pal[1:vcount(net)]

# 9. Plot the Network
plot(net, 
     layout = layout_fr,
     vertex.color = node_colors,
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     vertex.label.dist = 2,
     vertex.frame.color = "white",
     edge.color = transparent_color("grey", 30), # Now this function will work!
     edge.curved = 0.2,
     main = "Figure 2g: Cell-Cell Interaction Network")

#Task 8. Final assembly (mandatory)
# 1. Load All Necessary Libraries
library(readxl)
library(ggplot2)
library(pheatmap)
library(igraph)
library(dplyr)
library(gridExtra)
library(grid)

# 2. Global Setup
file_path <- "C:\\Users\\user\\Downloads\\hb_stage_2.xlsx"
hb_pal <- c("#4e79a7", "#8cd17d", "#e15759", "#fabfd2", "#a0cbe8", 
            "#59a14f", "#b07aa1", "#ff9d9a", "#f28e2b", "#f1ce63",
            "#79706e", "#d4a6c8", "#e9e9e9", "#ffbe7d", "#bab0ac",
            "#9d7660", "#d37295", "#86bcb6", "#362a39", "#cd9942")

# --- RECREATE PANELS AS OBJECTS ---

# Panel 2a: Boxplot
df_a <- read_excel(file_path, sheet = "a")
p2a <- ggplot(df_a, aes(x = reorder(cell_type, -new_ratio, FUN = median), y = new_ratio, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.5) + scale_fill_manual(values = hb_pal) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none") +
  labs(title = "a) Cell-type Ratios", x = "", y = "New Ratio")

# Panel 2b: Scatter
df_b <- read_excel(file_path, sheet = "b")
p2b <- ggplot(df_b, aes(x = log2(alpha), y = log2(half_life))) +
  geom_point(alpha = 0.2, color = "grey50") + theme_classic() +
  geom_text(data = filter(df_b, cell %in% c("Camp", "Ccr2")), aes(label=cell), color="red", fontface="bold") +
  labs(title = "b) Kinetic Landscape", x = "log2 Alpha", y = "log2 Half-life")

# Panel 2c: Gene Heatmap (Convert to Grob)
df_c <- read_excel(file_path, sheet = "c")
mat_c <- as.matrix(df_c[,-1]); rownames(mat_c) <- df_c[[1]]
p2c <- pheatmap(mat_c, cluster_cols = FALSE, scale = "row", silent = TRUE,
                color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                main = "c) Gene Expression")$gtable

# Panel 2d: Pathway Heatmap (Convert to Grob)
df_d <- read_excel(file_path, sheet = "d_1")
mat_d <- as.matrix(df_d[,-1]); rownames(mat_d) <- df_d[[1]]
p2d <- pheatmap(mat_d, cluster_cols = FALSE, cluster_rows = FALSE, silent = TRUE,
                color = colorRampPalette(c("#4e79a7", "white", "#e15759"))(50),
                main = "d) Pathway Activity")$gtable

# Panel 2e: Bubble Plot
df_e <- read_excel(file_path, sheet = "e")
p2e <- ggplot(df_e, aes(x = half_life, y = alpha, size = count, color = stage)) +
  geom_point(alpha = 0.6) + scale_x_log10() + scale_y_log10() + theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size=7)) +
  labs(title = "e) Functional Kinetics", x = "Half-life", y = "Alpha")

# Panel 2f: Stacked Bar
df_f <- read_excel(file_path, sheet = "f") %>% filter(stage %in% c("s00h", "s72h"))
p2f <- ggplot(df_f, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) + theme_classic() +
  labs(title = "f) Lineage Shift", x = "", y = "Prop.")

# Create the layout grid
# We arrange them in a 3-row layout
final_plot <- grid.arrange(
  p2a, p2b, 
  p2c, p2d, 
  p2e, p2f,
  ncol = 2,
  widths = c(1.2, 0.8),
  top = textGrob("Figure 2: Global Immune Response Kinetics", gp = gpar(fontsize=16, fontface="bold"))
)

# Export as high-resolution PNG
ggsave("Figure2_Publication_Ready.png", final_plot, width = 12, height = 16, dpi = 300)