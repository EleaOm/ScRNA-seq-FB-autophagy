
# Fibroblast Analysis ----------------------------------------------------------------


library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(ggmap)
library(genesorteR, quietly = TRUE)
library(clustree)
library(reshape2)
library(scProportionTest)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

sample="Fibroblast_integrated_"
path <- "FB_autophagy/"
path_fig <- "FB_autophagy/figures/"

#colorblind friendly color panels for all figures
#for conditions/stim https://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
# for testing colors: https://davidmathlogic.com/colorblind/
# palette generator https://loading.io/color/random/
sham.col = "#FC8D62"
tac.col = "#7991C7"
stim.col=c(sham.col,tac.col)
stim.col.light=c("#E8BEAE","#AEB7CE")
#for genotypes #From Paul Tol: https://personal.sron.nl/~pault/ '#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'
col.pdgfrb="#E26876"
col.col1a1="#59AD66"
col.gli1="#6361A0"
geno.cl = c(col.col1a1,col.pdgfrb,col.gli1)
geno.cl.light = c("#8DB193","#DE9FA6","#8281A2")
#for subclustering #RColorBrewer
display.brewer.pal(name = "Set2",n = 6)
cluster.col=brewer.pal(name = "Set2",n=6)
#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
half.gradient.col = brewer.pal(n = 9, name = "YlOrRd")
col.ramp<-colorRampPalette(gradient.col)

umap_theme <- theme(axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.title = element_text(size = 60),
                    plot.caption=element_text(size=40, face="bold", hjust = 0), strip.text = element_text(size = 60),
                    legend.text =element_text(size=40))
bar_theme <- theme(axis.text.x = element_text(size=20),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(size=20),axis.title.y = element_text(size = 20)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 20)) +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 30))

################################################----start of sample processing before analysis----##################################################
Samples.combined = readRDS(file = paste0(path, "integrated_FB_after_bacis_filters.rds"))
#add treatment based on sham or tac
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Gli1_TAC','TAC',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='PDGFRb_TAC','TAC',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Col1a1_TAC','TAC',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Gli1_Sham','Sham',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='PDGFRb_Sham','Sham',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Col1a1_Sham','Sham',Samples.combined@meta.data$treatment)
Samples.combined$treatment <- factor(Samples.combined$treatment,levels=c("Sham","TAC"))

#add genotype to metadata
# add column for only the genotype
Samples.combined$orig.geno=Samples.combined$orig.ident
Samples.combined$orig.geno[Samples.combined@meta.data$orig.ident %in% c("Gli1_TAC","Gli1_Sham")]="Gli1"
Samples.combined$orig.geno[Samples.combined@meta.data$orig.ident %in% c("Pdgfrb_TAC","Pdgfrb_Sham")]="Pdgfrb"
Samples.combined$orig.geno[Samples.combined@meta.data$orig.ident %in% c("Col1a1_TAC","Col1a1_Sham")]="Col1a1"
Samples.combined$orig.geno=factor(Samples.combined$orig.geno,levels = c("Col1a1","Pdgfrb","Gli1"))

View(Samples.combined@meta.data)

#first clustering
#----function to perform clustering based on seurat cca integartion------
recluster <- function(object,res,rep.str=5,loc.con=1){
  DefaultAssay(object) <- "integrated"
  # Run the standard workflow for visualization and clustering
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, verbose = FALSE)
  # t-SNE and Clustering
  object <- RunUMAP(object, reduction = "pca", dims = 1:30,repulsion.strength = rep.str,local.connectivity = loc.con)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object, resolution = res)
  print(DimPlot(object))
  message("plot not saved")
  DefaultAssay(object) <- "RNA"
  return(object)
}

Samples.combined = recluster(Samples.combined,0.5) # function (см. выше)

#plot before tdTom filtering steps
p1 <- DimPlot(Samples.combined, reduction = "umap", group.by = "orig.geno", pt.size = 1, cols = geno.cl.light) + 
  theme(legend.text=element_text(size=30)) &
  theme(text=element_text(size=30)) & 
  theme(legend.position = "bottom")
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "treatment", pt.size = 1, cols = stim.col.light) +
  theme(legend.text=element_text(size=30)) &
  theme(text=element_text(size=30)) & 
  theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(path_fig, sample,"before_processing_clusters_res0.5.jpeg"), width=20 , height = 10)

#cluster plots - samples
DimPlot(Samples.combined, reduction = "umap", split.by = "orig.geno",group.by = "sample", pt.size = 1) +
  theme(text=element_text(size=30)) +
  theme(legend.position = "bottom")
ggsave(filename = paste0(path_fig, sample,"_samples_split.jpeg"), width=30 , height = 10)

#cluster plots - sample
#add celltypes + treatment with levels
lvl=NULL
for (levels in levels(Samples.combined$orig.geno)) {lvl=c(lvl,paste0(levels,"_",levels(Samples.combined$treatment)))}
Samples.combined$celltypes.treat=factor(paste0(Samples.combined$orig.geno,"_",Samples.combined$treatment),levels = lvl)
DimPlot(Samples.combined,split.by = "celltypes.treat", reduction = "umap", pt.size = 0.1,ncol = 2) + NoLegend() +
  theme(text=element_text(size=30))
ggsave(filename = paste0(path_fig, sample,"_sample_split.jpeg"), width=10 , height = 20)

#-------------ploting only tdtom
#quality check on tdTom expression
FeaturePlot(Samples.combined, features = c("tdTomato-Green"),pt.size = 0.5)
ggsave(filename = paste0(path_fig, sample,"tdTom_unfiltered.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)
ggsave(filename = paste0(path_fig, sample,"tdTom_Vln_unfiltered.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)+theme(axis.text.y = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.5,size = 20),axis.title.x = element_blank(),axis.title.y = element_text(size = 25))+NoLegend()
ggsave(filename = paste0(path_fig, sample,"tdTom_Vln_unfiltered_v2.jpeg"), width=10 , height =10)
RidgePlot(Samples.combined,features = "tdTomato-Green",group.by = "orig.ident", layer = "counts",assay = "RNA")
ggsave(filename = paste0(path_fig, sample,"tdTom_reads_unfiltered.jpeg"), width=10 , height = 10)

saveRDS(Samples.combined, file = paste0(path, sample,'_sample.combined1.RDS'))

Samples.combined = readRDS(file = paste0 (path, "Fibroblast_integrated__sample.combined1.RDS"))

Samples.combined <- Samples.combined %>% Seurat::NormalizeData(assay = 'RNA', verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE)

#run Genesorter to get GeneProb of tdtom per cluster
sg = sortGenes(Samples.combined@assays$RNA@data, Idents(Samples.combined))
colnames(sg$condGeneProb) = paste0(levels(Idents(Samples.combined)))
Samples.combined_geneProb = as.data.frame(sg$condGeneProb)["tdTomato-Green",]
View(Samples.combined_geneProb)
remove(sg)

#------------basic nfeature and ncount plots
FeaturePlot(Samples.combined, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(path_fig, sample,"_nfeature_ncount_unfiltered.jpeg"), width=10 , height = 5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(path_fig, sample,"_nfeature_ncount__unfilteredvln.jpeg"), width=10 , height =5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,size = 20),axis.title.x = element_blank(),axis.text.y = element_text(size = 25))+NoLegend()
ggsave(filename = paste0(path_fig, sample,"_nfeature_unfilteredvln_v2.jpeg"), width=10 , height =10)
DimPlot(Samples.combined, reduction = "umap", group.by = "seurat_clusters", pt.size = 1, label = T, label.size = 10) + NoLegend() +
  theme(text=element_text(size=30))
ggsave(filename = paste0(path_fig, sample,"_umap_before.jpeg"), width=10 , height = 10)


# Testing for Mural markers, because fibroblasts and mural cell share Gli1
FeaturePlot(Samples.combined, features = "Myh11", min.cutoff = "q9",pt.size = 1,order = T)
ggsave(filename = paste0(path_fig, sample,"Myh11_before.jpeg"), width=10 , height = 10)
FeaturePlot(Samples.combined, features = "Vtn", min.cutoff = "q9",pt.size = 1,order = T)
ggsave(filename = paste0(path_fig, sample,"Vtn_before.jpeg"), width=10 , height = 10)
FeaturePlot(Samples.combined, features = "Rgs5", min.cutoff = "q9",pt.size = 1,order = T)
ggsave(filename = paste0(path_fig, sample,"Rgs5_before.jpeg"), width=10 , height = 10)
# Clusters 5, 7 and probably 11 are mural cells 

#----------subsetting to filter out tdtom low and nfeature low clusters, also exclude mural cells
Samples.combined=subset(Samples.combined,idents = c("4","5","7","9","10","11"),invert=T) # subset + invert = True - убирает то, что указано

DefaultAssay(Samples.combined) <- "RNA"
Samples.combined <- subset(Samples.combined,cells = WhichCells(Samples.combined,slot = "counts",expression = `tdTomato-Green` == 0),invert=T)
Samples.combined <- recluster(Samples.combined,0.3,rep.str = 1,loc.con = 20)

# Для русской версии: Ложнооперированные (ЛОЖ), коарктация аорты (КА)
Samples.combined[['treatment']] <- Samples.combined[['sample']]
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Gli1_TAC','КА',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Pdgfrb_TAC','КА',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Col1a1_TAC','КА',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Gli1_Sham','ЛОЖ',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Pdgfrb_Sham','ЛОЖ',Samples.combined@meta.data$treatment)
Samples.combined[['treatment']] <- ifelse(Samples.combined@meta.data$treatment =='Col1a1_Sham','ЛОЖ',Samples.combined@meta.data$treatment)
Samples.combined@meta.data$treatment <- factor(Samples.combined@meta.data$treatment,levels=c("ЛОЖ","КА"))


#plot after tdTom filtering step and before nfeature filter
p1 <- DimPlot(Samples.combined, reduction = "umap", group.by = "seurat_clusters", cols = cluster.col, label = T, pt.size = 1, label.size = 12) + NoLegend() +
  theme(text=element_text(size=30))
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "treatment", pt.size = 1, cols = stim.col.light) +theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=20)) +
  theme(text=element_text(size=30))
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(path_fig, sample,"_filtered_clusters_res0.3.jpeg"), width=20 , height = 10)

#------------basic nfeature ncount plots and tdtom
FeaturePlot(Samples.combined, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(path_fig, sample,"_nfeature_ncount.jpeg"), width=10 , height = 5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(path_fig, sample,"_nfeature_ncount_vln.jpeg"), width=10 , height =5)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)
ggsave(filename = paste0(path_fig, sample,"tdTom_Vln.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0,cols=cluster.col)+theme(axis.text.y = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.5,size = 20),axis.title.x = element_blank(),axis.title.y = element_text(size = 25))+NoLegend()
ggsave(filename = paste0(path_fig, sample,"tdTom_Vln_v2.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0,cols=cluster.col)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,size = 20),axis.title.x = element_blank(),axis.text.y = element_text(size = 25))+NoLegend()
ggsave(filename = paste0(path_fig, sample,"_nfeature_v2.jpeg"), width=10 , height =10)

#plot after processing  FINAL
DimPlot(Samples.combined, reduction = "umap", group.by = "seurat_clusters", cols = cluster.col, label = TRUE, pt.size = 1, label.size = 12) + 
  NoLegend() +
  umap_theme +
  theme(plot.title = element_text(size = 50)) +
  ggtitle("Кластеры фибробластов")
ggsave(filename = paste0(path_fig, sample,"after_processing_clusters_res0.3.jpeg"), width=10 , height = 10)
#no label version
DimPlot(Samples.combined, reduction = "umap",group.by = "seurat_clusters",cols = cluster.col,  pt.size = 1) + NoLegend() +
  theme_nothing()
ggsave(filename = paste0(path_fig, sample,"after_processing_clusters_res0.3_noLabel.jpeg"), width=10 , height = 10)

#cluster plots - genotype
DimPlot(Samples.combined, reduction = "umap", split.by = "orig.geno",group.by = "orig.geno", pt.size = 1,cols = geno.cl) + NoLegend() +
  umap_theme +
  ggtitle("Исходный генотип")
ggsave(filename = paste0(path_fig, sample,"_genotype_clusters_res0.3_split.jpeg"), width=30 , height = 10)
#no label version
DimPlot(Samples.combined, reduction = "umap", split.by = "orig.geno",group.by = "orig.geno", pt.size = 1,cols = geno.cl) + NoLegend() +
  theme_nothing()
ggsave(filename = paste0(path_fig, sample,"_genotype_clusters_res0.3_split_noLabel.jpeg"), width=30 , height = 10)


#-------normalize and scale RNA data
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined <- NormalizeData(Samples.combined, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE)

#--------save filtered object
saveRDS(Samples.combined,file = paste0(path, sample,"_filtered_processed.rds"))
################################################----end of sample processing before analysis----##################################################
########################################################----start of sample analysis----###########################################################
#----------load filtered object---------
Samples.combined <- readRDS(file = paste0(path, "Fibroblast_integrated__filtered_processed.rds"))

#-------clustertree for optimazation of cluster resolution-----
sc.int=Samples.combined
DefaultAssay(sc.int)<-"integrated"
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.3, by=0.1), verbose = FALSE)
clustree(sc.int)
ggsave(filename = paste0(path_fig, sample,"_cluster_resolution_tree.jpeg"), width=10 , height = 10)
remove(sc.int)

#-------annotation for Fib clusters according to marker genes, GO terms and ECM scoring-------
Idents(Samples.combined)="seurat_clusters"
# перед переименованием кластеров надо сначала проверить, какие гены вооббще экспрессируют кластеры (top10)

#------- top10 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "wilcox",min.pct = 0.3,assay = "RNA")
save(all.markers, file = paste0(path, sample,'all_markers.RData'))
all.markers <- load(file = paste0(path, sample,'all_markers.RData'))

library(data.table)
all.markers_table <- data.table(all.markers)
all.markers_table$pct.diff = all.markers_table$pct.1 - all.markers_table$pct.2
all.markers_table <- all.markers_table[, c("cluster","gene","avg_log2FC","pct.1","pct.2",
                                           "pct.diff","p_val","p_val_adj")]
write_tsv(all.markers_table, file = paste0(path, sample,'all_markers.tsv'))

#get top5 marker genes per cluster
all.markers_table <- read_tsv(file = paste0(path, sample,'all_markers.tsv'))

top5 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) %>% pull(gene))

top5 <- unique(all.markers_table %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) %>% pull(gene))

DotPlot(Samples.combined, group.by = "seurat_clusters",features = top5,dot.scale = 10,assay = "RNA") &
  coord_flip() & theme(axis.text.x = element_text(angle = 45,hjust=1),axis.text.y = element_text(face = "italic")) & 
  scale_colour_gradientn(colours = gradient.col)

top5[top5 %in% c("Maff", "Dnajb1", "Nr4a3",
                 "Hsd11b1", "Gsta3", 
                 "Cst6", "Crlf1", 
                 "Ly6c1", "Tmem100","Dpp4",
                 "Gpc3", "Inmt", "Tmem176a")] = c("Fosb", "Fos","Jun",
                                                  "Fap", "Smoc2", 
                                                  "Ctgf", "Col8a1",
                                                  "Cd248", "Cd55", "Jpt1",
                                                  "Ptn", "Adm", "Fgl2")
# 0 - stressed
# 1 - active/proliferative
# 2 - ECM
# 3 - adhesive
# 4 - perivascular
# 5 - intereferon

# exchanged some markers with well defined marker genes for clear annotation

#
types <- list("0"="Стрессовый ответ","1"="Пролиферирующие", "2"="Синтетические","3"="Адгезивные","4"="Регуляторные","5"="Ответ на интерферон")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes <- Idents(Samples.combined)
cell.levels = c("Стрессовый ответ","Пролиферирующие","Синтетические","Адгезивные","Регуляторные","Ответ на интерферон")
Samples.combined$celltypes = factor(Samples.combined$celltypes,levels = cell.levels)
#short abroviation of annotation
Idents(Samples.combined)="seurat_clusters"
types <- list("0"="Стресс","1"="Пролиферация","2"="Синтез","3"="Адгезия","4"="Регуляция","5"="Интерферон")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes.short <- Idents(Samples.combined)
cell.levels = c("Стресс","Пролиферация","Синтез","Адгезия","Регуляция","Интерферон")
Samples.combined$celltypes.short = factor(Samples.combined$celltypes.short,levels = cell.levels)
Idents(Samples.combined)="seurat_clusters"
#add celltypes + stim with levels
lvl=NULL
for (levels in levels(Samples.combined$celltypes)) {lvl=c(lvl,paste0(levels,"_",levels(Samples.combined$treatment)))}
Samples.combined$celltypes.treat=factor(paste0(Samples.combined$celltypes,"_",Samples.combined$treatment),levels = lvl)
#save after renaming
saveRDS(Samples.combined,file = paste0(path, sample,"_filtered_processed.rds"))

Samples.combined <- readRDS(file = paste0(path, sample,"_filtered_processed.rds"))
#vertikal
DotPlot(Samples.combined, group.by = "celltypes.short",features = top5,dot.scale = 10,assay = "RNA") &
  coord_flip() & 
  theme(axis.text.x = element_text(angle = 45,hjust=1),axis.text.y = element_text(face = "italic"),
        axis.title.x = element_blank(), axis.title.y = element_blank()) &
  scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(path_fig, sample,"_top5_markers_per_cluster_v.jpeg"), width=5.2, height = 9.5)
# no label:
DotPlot(Samples.combined, group.by = "celltypes.short",features = top5,dot.scale = 10,assay = "RNA") &
  coord_flip() & 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_blank()) &
  scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(path_fig, sample,"_top5_markers_per_cluster_noLabel.jpeg"), width=3, height = 8)

#dimplot short labels
DimPlot(Samples.combined,group.by = "celltypes.short", reduction = "umap",cols = cluster.col, label = TRUE, pt.size = 1, 
        label.size = 12, repel=TRUE) & NoLegend() & 
  umap_theme &
  theme(plot.title = element_blank())
ggsave(filename = paste0(path_fig, sample,"after_processing_clusters_res0.2_labeled_short.jpeg"), width=10 , height = 10)

#dimplot long labels
DimPlot(Samples.combined,group.by = "celltypes", reduction = "umap", label = F, pt.size = 1, cols = cluster.col) & 
  umap_theme &
  theme(plot.title = element_blank()) &
  guides(color = guide_legend(override.aes = list(size = 10)))
ggsave(filename = paste0(path_fig,sample,"after_processing_clusters_res0.2_labeled.jpeg"), width=16 , height = 10)

#dimplot stim split
DimPlot(Samples.combined,group.by = "treatment",split.by="treatment", reduction = "umap",cols = stim.col, pt.size = 1) & NoLegend() &
  umap_theme &
  labs(title = "Воздействие", caption = "ЛОЖ - ложнооперированные\nКА - коарктация аорты")
ggsave(filename = paste0(path_fig, sample,"split_treatment.jpeg"), width=18 , height = 10)


# selective markers of fibroblasts population
FeaturePlot(Samples.combined, reduction = "umap", features = c('Fos', 'Fap', "Postn", 'Jpt1', 'Adm', 'Ifit3'), order = T,
            pt.size = 1, ncol=3, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  theme(axis.title.y.right = element_text(size = 12), axis.title = element_text(size = 15), axis.text = element_text(size = 10)) & 
  theme(legend.position = "right")
ggsave(filename = paste0(path_fig, sample,"_umap_ECM-Fib.jpeg"), width=10 , height = 5)
# no label:
FeaturePlot(Samples.combined, reduction = "umap", features = c('Fos', 'Fap', "Postn", 'Jpt1', 'Adm', 'Ifit3'), order = T,
            pt.size = 1, ncol=3, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  theme_nothing() & 
  theme(legend.position = "right")
ggsave(filename = paste0(path_fig, sample,"_umap_ECM-Fib_noLabel.jpeg"), width=10 , height = 5)


# selective markers of fibroblasts population
FeaturePlot(Samples.combined, reduction = "umap", features = c("Fap", "Ddr2", "Col1a1"), order = T,
            pt.size = 1, ncol=3, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  theme(axis.title.y.right = element_text(size = 12), axis.title = element_text(size = 15), axis.text = element_text(size = 10)) & 
  theme(legend.position = "right")
ggsave(filename = paste0(path_fig, sample,"_umap_FB.jpeg"), width=12 , height = 4)
# no label:
FeaturePlot(Samples.combined, reduction = "umap", features = c("Fap", "Ddr2", "Col1a1"), order = T,
            pt.size = 1, ncol=3, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  theme_nothing() & 
  theme(legend.position = "right")
ggsave(filename = paste0(path_fig, sample,"_umap_FB_noLabel.jpeg"), width=12 , height = 4)


#----barplot of genotype contribution per cluster (normalized to even input per genotype)-----
breakdown<-table(Samples.combined@meta.data$seurat_clusters, Samples.combined@meta.data$orig.geno) 
total.col1a1 = sum(breakdown[,1])
total.gli1 = sum(breakdown[,2])
total.pdgfrb = sum(breakdown[,3])
ratio.col1a1.gli1 = total.col1a1/total.gli1
ratio.col1a1.pdgfrb = total.col1a1/total.pdgfrb
breakdown[,2] = round(breakdown[,2]*ratio.col1a1.gli1)
breakdown[,3] = round(breakdown[,3]*ratio.col1a1.pdgfrb)
breakdown=t(breakdown)
breakdown <- round(apply(breakdown, 2, function(x){x*100/sum(x)}),2)
breakdown.df = as.data.frame(breakdown)
breakdown.df.new = breakdown.df[1,]
breakdown.df.new[2,] = breakdown.df[3,]
breakdown.df.new[3,] = breakdown.df[2,]
breakdown.df = breakdown.df.new
breakdown.df = melt(t(breakdown.df))
breakdown.df$Var2 = factor(breakdown.df$Var2,levels = levels(Samples.combined$orig.geno))

ggplot(data=breakdown.df, aes(x=Var1, y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  scale_x_continuous(breaks = as.numeric(levels(factor(breakdown.df$Var1))),labels = levels(factor(breakdown.df$Var1)),name="Cluster")+
  scale_y_continuous(name="Доля клеток, %") +
  scale_fill_manual( values = c("Col1a1" = col.col1a1,  "Pdgfrb" = col.pdgfrb,"Gli1" = col.gli1))+
  bar_theme +
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.title = element_blank(),legend.key.size = unit(1, "cm"))
ggsave(filename = paste0(path_fig, sample,"_Sample_contribution_per_custom_cluster_normalized.jpeg"), width=5 , height =5)

# no label:
ggplot(data=breakdown.df, aes(x=Var1, y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  scale_x_continuous(breaks = as.numeric(levels(factor(breakdown.df$Var1))),labels = levels(factor(breakdown.df$Var1)),name="Cluster")+
  scale_y_continuous(name="Доля клеток, %") +
  scale_fill_manual( values = c("Col1a1" = col.col1a1,  "Pdgfrb" = col.pdgfrb,"Gli1" = col.gli1))+
  theme_nothing() +
  theme(legend.position = "right") +
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.title = element_blank(),legend.key.size = unit(1, "cm"), legend.text = element_blank())
ggsave(filename = paste0(path_fig, sample,"_Sample_contribution_per_custom_cluster_normalized_noLabel.jpeg"), width=5 , height =5)


#-------barplot of condition contribution per cluster (input normalized)-----
breakdown<-table(Samples.combined@meta.data$seurat_clusters, Samples.combined@meta.data$treatment) 
breakdown[,1]=100*(breakdown[,1]/sum(breakdown[,1]))
breakdown[,2]=100*(breakdown[,2]/sum(breakdown[,2]))
breakdown=t(breakdown)
breakdown <- apply(breakdown, 2, function(x){x*100/sum(x)})
breakdown = round(as.data.frame(breakdown),digits = 2)
breakdown.df = as.data.frame(breakdown)
breakdown.df = melt(t(breakdown.df))
breakdown.df$Var2=factor(breakdown.df$Var2,levels = c("КА","ЛОЖ"))
celltypes.short = c("Стресс","Пролиферация","Синтез","Адгезия","Регуляция","Интерферон")

ggplot(data=breakdown.df, aes(x=rev(Var1), y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=c(celltypes.short,celltypes.short),y=1),hjust = 0, color="black", size=10)+
  scale_y_continuous(name="Доля клеток, %") +
  coord_flip() + 
  scale_fill_manual(values = rev(stim.col))+
  labs(caption = "ЛОЖ - ложнооперированные\nКА - коарктация аорты") +
  theme(axis.title.x = element_text(size = 20) ,axis.text.x = element_text(size=20), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 20)) +
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.caption = element_text(size = 20, hjust = 0), legend.text = element_text(size = 20)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_hline(yintercept=50, linetype="dashed", color = "black")
ggsave(filename = paste0(path_fig, sample,"_barplot_percentage_per_stim and cl.jpeg"), width=5 , height = 8)

# no label
ggplot(data=breakdown.df, aes(x=rev(Var1), y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  coord_flip() + 
  scale_y_continuous(name="Доля клеток, %") +
  scale_fill_manual(values = rev(stim.col))+
  theme_nothing() +
  theme(legend.position = "right") +
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.caption = element_text(size = 20, hjust = 0), legend.text = element_text(size = 20)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.text = element_blank(), legend.title = element_blank())+
  geom_hline(yintercept=50, linetype="dashed", color = "black")
ggsave(filename = paste0(path_fig, sample,"_barplot_percentage_per_stim and cl_noLabel.jpeg"), width=5 , height = 7)


#-------proportion analysis-------
prop_test <- sc_utils(Samples.combined)

prop_test <- permutation_test(
  prop_test, cluster_identity = "celltypes.short",
  sample_1 = "ЛОЖ", sample_2 = "КА",
  sample_identity = "treatment",
  n_permutations = 10000
)
permutation_plot(prop_test) & 
  theme(legend.position = "bottom", axis.title.y.left = element_blank(), 
                                    axis.title = element_text(size = 15), axis.text = element_text(size = 10)) & 
  ggtitle('Пермутационный тест: КА vs ЛОЖ')
ggsave(filename = paste0(path_fig, sample,"_proportionTest_sham_vs_tac.jpeg"), width=4.5, height = 4)
