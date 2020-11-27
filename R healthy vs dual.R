if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("flowCore")
BiocManager::install("CATALYST")

library(readxl)
library(SingleCellExperiment)
library(diffcyt)
library(HDCytoData)
library(DT)
library(ggplot2)
library(cowplot)
library(flowCore)
library(CATALYST)
library(uwot)



###### Input data: not activated samples
### Patients metadata excel import 
setwd("C:/Users/Kili/Desktop/Rechts der Isar/R/healthy vs dual")
md <- read_excel("meta_healthy_dual_baseline.xlsx")
DT::datatable(data.frame(md))

### Markers panel excel import
#Activation markers are indicated as "state", the rest - "type". This is a very important stratification as usually lineage markers (type) are used for clustering while functional markers (state) are used for the differential expression analysis. In the furter analysis we use CD63, CD107a, PAC1 and CD62P as state markers.
panel <- read_excel("panel_umap.xlsx")
datatable(data.frame(panel))

### .fcs files df import
files <- list.files(path = "C:/Users/Kili/Desktop/Rechts der Isar/R/healthy vs dual/input baseline/", pattern = "\\.fcs$", full.names = TRUE)
#files <- files[-c(7)]
fs <- read.flowSet(files, transformation = FALSE, truncate_max_range = FALSE, which.lines = 1)
fs = fs[,panel$fcs_colname]

### Saving the data for further analysis in Python
outDir = "C:/Users/Kili/Desktop/Rechts der Isar/R/healthy vs dual/output baseline/"
for (i in 1:length(fs)){
  write.flowSet(fs[i], outDir)
}

### Checking if panel file fits the column names in fcs files:
all(panel$fcs_colname %in% colnames(fs))

### Prepare metadata and check if the filenames match:
md$condition <- factor(md$condition, levels = c("healthy", "patient"))
#md$treatment <- factor(md$treatment, levels = c("no", "mono", "dual", "triple", "loaded"))
md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$condition)])
#md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$treatment)])

ids1 <- fsApply(fs, identifier)
md = subset(md,file_name %in% ids1)

### Constracting SingleCellExperiment object: 
sce <- prepData(fs, panel, md, features = panel$fcs_colname, cofactor = 5)
set.seed(1234)
library(uwot)
sce <- runDR(sce, "UMAP", cells = min(n_cells(sce)), features = "type")

 
### Saving the UMAP dimentions 
dr = "UMAP"
color_by = "condition"
dims = c(1,2)
xy <- reducedDim(sce, dr)[, dims]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(sce), xy)
df <- df[!(is.na(df$x) | is.na(df$y)), ]
outDir = "C:/Users/Kili/Desktop/Rechts der Isar/R/healthy vs dual/umap_baseline.csv"
write.csv(df,outDir, row.names = FALSE)

cs = c('#9AB8C8','#DBA794')
plotDR(sce, "UMAP", color_by = "condition", k_pal = cs) 


{r, fig.width = 10, fig.height = 8}
plot_grid(plotDR(sce, "UMAP", color_by = "CD63",facet_by = "condition"),
          plotDR(sce, "UMAP", color_by = "CD107a",facet_by = "condition"),
          plotDR(sce, "UMAP", color_by = "CD62P",facet_by = "condition"),
          plotDR(sce, "UMAP", color_by = "CD154",facet_by = "condition"),
          ncol = 2)








###########

## Activated samples
md <- read_excel("meta_healthy_dual_activated.xlsx")
DT::datatable(data.frame(md))

### Markers panel
#Activation markers are indicated as "state", the rest - "type". This is a very important stratification as usually lineage markers (type) are used for clustering while functional markers (state) are used for the differential expression analysis. In the furter analysis we use CD63, CD107a, PAC1 and CD62P as state markers and the rest as type markers.
panel <- read_excel("panel_umap.xlsx")
datatable(data.frame(panel))

### .fcs files
files <- list.files(path = "C:/Users/Kili/Desktop/Rechts der Isar/R/healthy vs dual/input activated/", pattern = "\\.fcs$", full.names = TRUE)
fs <- read.flowSet(files, transformation = FALSE, truncate_max_range = FALSE)
fs = fs[,panel$fcs_colname]


#Saving the data for further analysis in Python
outDir = "C:/Users/Kili/Desktop/Rechts der Isar/R/healthy vs dual/output activated/"
for (i in 1:length(fs)){
  write.flowSet(fs[i], outDir)}


#Checking if panel file fits the column names in fcs files:
all(panel$fcs_colname %in% colnames(fs))


#Prepare metadata and check if the filenames match:
md$condition <- factor(md$condition, levels = c("healthy", "patient"))
md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$condition)])


ids1 <- fsApply(fs, identifier)
md = subset(md,file_name %in% ids1)


#Constracting SingleCellExperiment object: 
sce <- prepData(fs, panel, md, features = panel$fcs_colname, cofactor = 5)
set.seed(1234)
sce <- runDR(sce, "UMAP", cells = min(n_cells(sce)), features = "type")

#Saving the UMAP dimentions 
dr = "UMAP"
color_by = "condition"
dims = c(1,2)
xy <- reducedDim(sce, dr)[, dims]
colnames(xy) <- c("x", "y")
df <- data.frame(colData(sce), xy)
df <- df[!(is.na(df$x) | is.na(df$y)), ]
outDir = "C:/Users/Kili/Desktop/Rechts der Isar/R/healthy vs dual/umap_activated.csv"
write.csv(df,outDir, row.names = FALSE)

plotDR(sce, "UMAP", color_by = "condition")


{r, fig.width = 10, fig.height = 8}
plot_grid(plotDR(sce, "UMAP", color_by = "CD63",facet_by = "condition"),
          plotDR(sce, "UMAP", color_by = "CD107a",facet_by = "condition"),
          plotDR(sce, "UMAP", color_by = "CD62P",facet_by = "condition"),
          plotDR(sce, "UMAP", color_by = "CD154",facet_by = "condition"),
          ncol = 2)


