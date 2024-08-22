
# Filtration and normalization --------------------------------------------

# 1. Create Seurat object and Filtering raw matrix .h5

# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DropletUtils")
# install.packages('Seurat') 

library(Seurat)
library(DropletUtils)
library(tidyverse)

# Настройка рабочей дирректории (ввести папку, в которой лежат файлы и будут сохраняться новые)
getwd()

path <- "FB_autophagy/"
path_fig <- "FB_autophagy/figures/"

# (1) Read 10X matrix
install.packages("hdf5r")
options(Seurat.object.assay.version = "v3")


Pdgfrb_Sham <- Read10X_h5('h5 matrixes/Pdgfrb_Sham_raw_gene_bc_matrix.h5', use.names = T)
Pdgfrb_TAC <- Read10X_h5('h5 matrixes/Pdgfrb_TAC_raw_gene_bc_matrix.h5', use.names = T)
Col1a1_Sham <- Read10X_h5('h5 matrixes/Col1a1_Sham_raw_gene_bc_matrix.h5', use.names = T)
Col1a1_TAC <- Read10X_h5('h5 matrixes/Col1a1_TAC_raw_gene_bc_matrix.h5', use.names = T)
Gli1_Sham <- Read10X_h5('h5 matrixes/Gli1_Sham_raw_gene_bc_matrix.h5', use.names = T)
Gli1_TAC <- Read10X_h5('h5 matrixes/Gli1_TAC_raw_gene_bc_matrix.h5', use.names = T)


# (2) Empty drops - так как используются нефильтрованный матрицы (raw_matrix), то 
# надо тщательно выполнить QC и удалить пустые капли, мертвые клетки и т.д.

EmtyDrops_h5 <- emptyDropsCellRanger(Pdgfrb_Sham, 
                                     n.expected.cells = 16000, # сколько клеток загружали в анализ
                                     max.percentile = 0.99,
                                     max.min.ratio = 10,
                                     umi.min = 500,
                                     umi.min.frac.median = 0.01,
                                     cand.max.n = 20000,
                                     ind.min = 45000,
                                     ind.max = 90000,
                                     round = TRUE,
                                     niters = 10000,
)
sum(EmtyDrops_h5$FDR <= 0.001, na.rm=TRUE)

Pdgfrb_Sham <- Pdgfrb_Sham[,which(EmtyDrops_h5$FDR <= 0.001)] # получаем новую матрицу, 
# отфильтрованную по пустым каплям

EmtyDrops_h5 <- emptyDropsCellRanger(Pdgfrb_TAC, n.expected.cells = 16000)
Pdgfrb_TAC <- Pdgfrb_TAC[,which(EmtyDrops_h5$FDR <= 0.001)]

EmtyDrops_h5 <- emptyDropsCellRanger(Gli1_Sham, n.expected.cells = 16000)
Gli1_Sham <- Gli1_Sham[,which(EmtyDrops_h5$FDR <= 0.001)]

EmtyDrops_h5 <- emptyDropsCellRanger(Gli1_TAC, n.expected.cells = 16000)
Gli1_TAC <- Gli1_TAC[,which(EmtyDrops_h5$FDR <= 0.001)]

EmtyDrops_h5 <- emptyDropsCellRanger(Col1a1_Sham, n.expected.cells = 16000)
Col1a1_Sham <- Col1a1_Sham[,which(EmtyDrops_h5$FDR <= 0.001)]

EmtyDrops_h5 <- emptyDropsCellRanger(Col1a1_TAC, n.expected.cells = 16000)
Col1a1_TAC <- Col1a1_TAC[,which(EmtyDrops_h5$FDR <= 0.001)]

rm(EmtyDrops_h5)


# (3) добавить свои матрицы с каунатми и клетками и сделать их объектом Seurat
Pdgfrb_Sham <- CreateSeuratObject(counts = Pdgfrb_Sham, project = "Pdgfrb_Sham", min.features = 100)
Pdgfrb_TAC <- CreateSeuratObject(counts = Pdgfrb_TAC, project = "Pdgfrb_TAC", min.features = 100)
Gli1_Sham <- CreateSeuratObject(counts = Gli1_Sham, project = "Gli1_Sham", min.features = 100)
Gli1_TAC <- CreateSeuratObject(counts = Gli1_TAC, project = "Gli1_TAC", min.features = 100)
Col1a1_Sham <- CreateSeuratObject(counts = Col1a1_Sham, project = "Col1a1_Sham", min.features = 100)
Col1a1_TAC <- CreateSeuratObject(counts = Col1a1_TAC, project = "Col1a1_TAC", min.features = 100)

gc()

# (4) Объединим образцы в один объект, чтобы работать сразу со всеми = merged Seurat object
# install.packages("devtools")
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)

merged_seurat_FB <- merge(x = Pdgfrb_Sham, y = c( Pdgfrb_TAC, Gli1_Sham, Gli1_TAC, Col1a1_Sham, Col1a1_TAC), 
                       add.cell.id = c("Pdgfrb_Sham", "Pdgfrb_TAC", "Gli1_Sham", "Gli1_TAC", "Col1a1_Sham", "Col1a1_TAC")) 
# добавим идентификаторы к баркодам

# Проверим, что идентификаторы применились правильно
head(merged_seurat_FB@meta.data)
tail(merged_seurat_FB@meta.data)

saveRDS(merged_seurat_FB, file= paste0(path, "merged_seurat_FB_unfiltered.rds"))


rm(Pdgfrb_Sham, Pdgfrb_TAC, Gli1_Sham, Gli1_TAC, Col1a1_Sham, Col1a1_TAC)

# 2. Оценка качества и фильтрация

# Generating quality metrics

# Посмотрим, что у нас в метаданных объединенного объекта

merged_seurat_FB <- readRDS("FB_autophagy/merged_seurat_FB_unfiltered.rds")
View(merged_seurat_FB@meta.data)

# (1) Для фильтрации низкокачественных данных надо добавить доп. метрики

# Добавим метрику Гены/UMI (чем больше генов на UMI, тем более сложные данные)
merged_seurat_FB$log10GenesPerUMI <- log10(merged_seurat_FB$nFeature_RNA) / log10(merged_seurat_FB$nCount_RNA)

# Добавим метрику процент митохондриальных генов и разделим на 100 (чтобы в долях было, чисто для удобства)
merged_seurat_FB$mitoRatio <- PercentageFeatureSet(object = merged_seurat_FB, pattern = "^mt-")

# Добавим метрику процент рибосомальных генов и разделим на 100
merged_seurat_FB$riboRatio <- PercentageFeatureSet(object = merged_seurat_FB, pattern = "^Rp[sl]")


# Для того, чтобы без риска добавить доп. метрики в метаданные, "выведем" отдельно эти метаданные из объекта сеурат в отдельный датафрэйм
metadata <- merged_seurat_FB@meta.data
head(metadata)

# Добавим ID клеток в метаданные
metadata$cells <- rownames(metadata)

# Переименуем колонки
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Создание колонки с образцами
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^Pdgfrb_Sham_"))] <- "Pdgfrb_Sham"
metadata$sample[which(str_detect(metadata$cells, "^Pdgfrb_TAC_"))] <- "Pdgfrb_TAC"
metadata$sample[which(str_detect(metadata$cells, "^Gli1_Sham_"))] <- "Gli1_Sham"
metadata$sample[which(str_detect(metadata$cells, "^Gli1_TAC_"))] <- "Gli1_TAC"
metadata$sample[which(str_detect(metadata$cells, "^Col1a1_Sham_"))] <- "Col1a1_Sham"
metadata$sample[which(str_detect(metadata$cells, "^Col1a1_TAC_"))] <- "Col1a1_TAC"


# Создание колонки Treatment
metadata$treatment <- NA
metadata$treatment[which(str_detect(metadata$sample, "Sham"))] <- "Sham"
metadata$treatment[which(str_detect(metadata$sample, "TAC$"))] <- "TAC"
head(metadata)

# Добавим метаданные обратно в объект сеурат
merged_seurat_FB@meta.data <- metadata

# (2) Визуализация метрик
# Визуализация числа клеток

png(paste0(path_fig, 'QC_unfiltered_Cells.png'), width = 1000, height = 1000)
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()
#Видим, что клеток в образце PDGFRB_TAC больше ожидаемого количества (>16000), 
# хотя функция EmptyDrops была применена

# Визуализация числа UMI/транскриптов на клетку (должно быть больше 500, чем больше, тем более глубокое прочтение)
png(paste0(path_fig, 'QC_unfiltered_UMI_per_cells.png'), width = 1000, height = 1000)
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) +
  ggtitle("UMIs per cell")
dev.off()
# Видим, что у нас преимущественно оч глубокое прочтение, отлично

# Визуализация распределения генов, детектированных на клетку
png(paste0(path_fig, 'QC_unfiltered_GENES_per_cells.png'), width = 1000, height = 1000)
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500) +
  ggtitle("Genes per cell")
dev.off()
# тут вопрос, какой трешхолд поставить? стандартно 300 < nFeature=nGenes < 5000.
# должен быть один большой пик, который представляет инкапсулированные клетки.
# если появляется плечо справа от пика или наблюдается бимодальное распределение, то это может быть из-за:
# - проблем с подачей клеток в прибор и др. причинами эксперимента
# - наличия иной клеточной популяции (покоящиеся клетки, клетки различающиеся по размеру)
# также мы видим, что образец Myh11_Sham имеет очень высокий уровень генов на клетку, возможно, дублеты?

# Визуализация между генами и числом UMI для определения присутствия клеток с небольшим количеством генов / UMIs
png(paste0(path_fig, 'QC_unfiltered_genes_UMI_mito.png'), width = 1000, height = 1000)
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample) +
  ggtitle("MitoRatio: Genes and UMIs")
dev.off()
# видим, что клетки, у которых низкое число генов (<800) и UMI имеют выский процент митохондриальных генов, т.е. клетки погибшие

png(paste0(path_fig, 'QC_unfiltered_genes_UMI_ribo.png'), width = 1000, height = 1000)
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=riboRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "yellow", high = "red") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample) +
  ggtitle("RiboRatio: Genes and UMIs")
dev.off()
# видим, что в некоторых образцах есть высокий процент рибосомальных генов >40 % (следовательно, пустые капли ушли не полностью)

# Визуализация экспрессии митохондриальных генов на клетку
png(paste0(path_fig, 'QC_unfiltered_mito_ratio.png'), width = 1000, height = 1000)
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 10) +
  ggtitle("mitoRatio")
dev.off()
# сделаем отсечку 10 = сохраним клетки, содержащие <10% мит. генов

# Визуализация сложности (complexity) экспрессии генов путем визуализации генов на UMI
png(paste0(path_fig, 'QC_unfiltered_complexity.png'), width = 1000, height = 1000)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  ggtitle("Compexity")
dev.off()
# видим, что наши данные обладают высокой "сложностью", однако есть аутлаеры

# Теперь можно фильтровать

#  (3) Doublets finder перед фильтрацией данных, но после EmtyDrops Removal

# BiocManager::install("scDblFinder")
library(scDblFinder)
library(SingleCellExperiment)
# BiocManager::install("BiocParallel")
library(BiocParallel)

## Convert object into singlecellexperiment
merged.sce <- as.SingleCellExperiment(merged_seurat_FB)
rm(merged_seurat_FB)
merged.sce <- scDblFinder(merged.sce, samples="sample", clusters=FALSE)
## Convert sce object back to seurat
merged_seurat_FB <- as.Seurat(merged.sce, counts = "counts", data = "logcounts")

table(merged.sce$scDblFinder.class) # всего было удалено 8764 дублетов, осталось 70208 синглетов, т.е. 12.5% дублетов
rm(merged.sce)

names(x = merged_seurat_FB[[]])
# удалим лишнюю информацию
merged_seurat_FB@meta.data <- merged_seurat_FB@meta.data[, !grepl("nCount_", colnames(merged_seurat_FB@meta.data))]
merged_seurat_FB@meta.data <- merged_seurat_FB@meta.data[, !grepl("nFeature_", colnames(merged_seurat_FB@meta.data))]
names(x = merged_seurat_FB[[]])


# заново пройтись по QC п. (2)

#  (4) Cell-level filtering
# Удаление низкокачественных ридов, используя выбранные трешхолды
filtered_seurat_FB <- subset(x = merged_seurat_FB, 
                          subset= (scDblFinder.class  == "singlet") & 
                            (nUMI >= 500) & 
                            (nGene >= 500) & # верхнюю границу пока ставить не будем, попробуем позже
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 10))
names(x = filtered_seurat_FB[[]])
filtered_seurat_FB@meta.data <- filtered_seurat_FB@meta.data[, !grepl("scDblFinder", colnames(filtered_seurat_FB@meta.data))]

rm(merged_seurat_FB)

# (5) Gene-level filtering
# Удаление генов, имеющих нулевую экспрессию во всех клетках и удаление генов, экспрессирующихся менее, чем в 10 клетках

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat_FB, layer = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat_FB <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat_FB@meta.data)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
png(paste0(path_fig, 'QC_immune_and_erythro.png'), width = 1000, height = 1000)
VlnPlot(filtered_seurat_FB,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
dev.off()

rm(counts)
rm(filtered_counts)
rm(merged_seurat)
rm(nonzero)
rm(keep_genes)
rm(metadata)


filtered_seurat_FB=subset(filtered_seurat_FB,cells = WhichCells(filtered_seurat_FB,expression = `Hba-a1` ==0))
filtered_seurat_FB=subset(filtered_seurat_FB,cells = WhichCells(filtered_seurat_FB,expression = `Hba-a2` ==0))
filtered_seurat_FB=subset(filtered_seurat_FB,cells = WhichCells(filtered_seurat_FB,expression = `Hbb-bs` ==0))
filtered_seurat_FB=subset(filtered_seurat_FB,cells = WhichCells(filtered_seurat_FB,expression = `Ptprc` ==0))

# Сохраним объект
saveRDS(filtered_seurat_FB, file=paste0(path, "seurat_FB_filtered.rds"))


## (6) Re-assess QC metrics (пройдемся заново, чтобы убедиться, что мы все правильно сделали)
# Save filtered subset to new metadata
metadata_clean <- filtered_seurat_FB@meta.data

# Визуализация числа клеток
png(paste0(path_fig, 'QC_filtered_CELLS.png'), width = 1000, height = 1000)
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

# Визуализация числа UMI/транскриптов на клетку (должно быть больше 500, чем больше, тем более глубокое прочтение)
png(paste0(path_fig,'QC_filtered_UMI_per_cell.png'), width = 1000, height = 1000)
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) +
  ggtitle("UMIs per cell")
dev.off()

# Визуализация распределения генов, детектированных на клетку
png(paste0(path_fig,'QC_filtered_GENES_per_cell.png'), width = 1000, height = 1000)
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 10000) +
  ggtitle("Genes per cell")
dev.off()
# Если есть клетки, содержащие больше 10000 генов, то это, вероятно, дублеты, тогда надо ставить второй фильтр <10000

# Визуализация между генами и числом UMI для определения присутствия клеток с небольшим количеством генов / UMIs
png(paste0(path_fig, 'QC_filtered_genes_UMI_mito.png'), width = 1000, height = 1000)
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample) +
  ggtitle("MitoRatio: Genes and UMIs")
dev.off()

# Визуализация экспрессии митохондриальных генов на клетку
png(paste0(path_fig, 'QC_filtered_mito_ratio.png'), width = 1000, height = 1000)
metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 10) +
  ggtitle("mitoRatio")
dev.off()
# сделаем отсечку 10 = сохраним клетки, содержащие 10% мит. генов

# Визуализация рибосомальных генов
png(paste0(path_fig, 'QC_filtered_genes_UMI_ribo.png'), width = 1000, height = 1000)
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=riboRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "yellow", high = "red") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample) +
  ggtitle("RiboRatio: Genes and UMIs")
dev.off()

# Визуализация сложности (complexity) экспрессии генов путем визуализации генов на UMI
png(paste0(path_fig,'QC_filtered_complexity.png'), width = 1000, height = 1000)
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  ggtitle("Compexity")
dev.off()

rm(metadata)
rm(metadata_clean)


# 3. Нормализация, PCA и трансформация ###

# Normalization, variance stabilization, and regression of unwanted variation for each sample
# ПО sctransform для нормализации и стабилизации вариаций "удаляет" шум глубины секвенирования (nUMI)
# другим шумом может быть разные стадии клеточного цикла, митохондриальные гены (зависит от эксперимента), 
# данное ПО не вычленяет это.
# Следовательно, перед дальнейшим анализом, избавимся от шумов путем regressed out

#---normalize and findvariable features-----------------------
# Разделим объект Сеурат по условиям и выполним нормализацию и необходимую регрессию
split_seurat_FB <- SplitObject(filtered_seurat_FB, split.by = "sample")
rm(filtered_seurat_FB)
split_seurat_FB <- split_seurat_FB[c("Pdgfrb_Sham", "Pdgfrb_TAC", 
                               "Gli1_Sham", "Gli1_TAC","Col1a1_Sham","Col1a1_TAC")]

# регрессия - затратный процесс, увеличим лимиты оперативной памяти
options(future.globals.maxSize = 4000 * 1024^2)


# BiocManager::install('glmGamPoi')
# [1] = Pdgfrb_Sham, [2] = Pdgfrb_TAC etc.
for (i in 1:length(split_seurat_FB)) {
  split_seurat_FB[[i]] <- SCTransform(split_seurat_FB[[i]],                  
                                   do.correct.umi = TRUE,
                                   vars.to.regress = c("riboRatio", "mitoRatio"),
                                   ncells = 5000, #
                                   do.scale = FALSE, # производитель не рекомендует шкалировать данные
                                   do.center = TRUE,
                                   variable.features.n = 3000) # стандартный параметр
}

saveRDS(split_seurat_FB, file=paste0(path, "seurat_normalized.RData"))

################################################################
################ 4. Интеграция образцов ########################
################################################################

#all datasets containing tdTomato tranced Fibroblasts cells
#perform integration of Pdgfrb, Col1a1, Gli1 Sham+TAC each (+28d gli1 tac)

integ_features <- SelectIntegrationFeatures(object.list = split_seurat_FB, 
                                            nfeatures = 3000) 

split_seurat_FB <- PrepSCTIntegration(object.list = split_seurat_FB, 
                                   anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat_FB, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

#exclude tdtomato from integration features
anchor.list = integ_anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
integ_anchors@anchor.features=anchor.list
Samples.combined_FB <- IntegrateData(anchorset = integ_anchors, dims = 1:20, verbose = TRUE)
saveRDS(Samples.combined_FB, file = paste0(path, "integrated_FB_after_bacis_filters.rds"))

