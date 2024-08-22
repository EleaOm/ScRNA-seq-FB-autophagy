
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# Autophagy Analysis ------------------------------------------------------
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(ggmap)
library(reshape2)
library(scProportionTest)


sample="Fibroblast_autophagy_"
path <- "FB_autophagy/"
path_fig <- "FB_autophagy/figures/"

# colors
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

# graphic themes
feature_theme <-  theme(axis.title.y.right = element_text(size = 20), axis.title = element_text(size = 20), 
                        axis.text = element_text(size = 20), legend.position = "right", legend.text =element_text(size=15), 
                        plot.title = element_text(size = 40))
umap_theme <- theme(axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.title = element_text(size = 60),
                    plot.caption=element_text(size=40, face="bold", hjust = 0), strip.text = element_text(size = 60),
                    legend.text =element_text(size=40))
vinplot_theme <- theme(axis.title.y = element_text(size = 20), axis.title.x = element_blank(),
                       axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), 
                       title = element_text(size = 40),
                       legend.position = "right", legend.text = element_text(size = 30))
bar_theme <- theme(axis.text.x = element_text(size=20, face = "italic"),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(size=20),axis.title.y = element_text(size = 20)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 20)) +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 30))


# set function for counting cell percent
PrctCellExpringGene <- function(object, genes, group.by = "all", assay = "RNA", datatype = "counts", threshold = 0){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object, assay = assay, datatype=datatype, threshold=threshold))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    results = lapply(list, PrctCellExpringGene, genes=genes, assay = assay, datatype=datatype, threshold=threshold)
    results %>% reduce(full_join, by="Markers") %>% select(any_of("Markers")) -> genelist
    results %>% reduce(full_join, by="Markers") %>% select(!any_of("Markers")) %>% "colnames<-"(names(results)) -> percentages
    combined <- cbind(genelist,percentages)
    return(combined)
  }
}

calc_helper <- function(object,genes,assay,datatype,threshold){
  counts = slot(object[[assay]],datatype)
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    round((sum(counts[genes,]>threshold)/ncells)*100, 1)
  }else{return(NA)}
}



Samples.combined <- readRDS(file = paste0(path, "Fibroblast_integrated__filtered_processed.rds"))

target_genes=c('Plau', 'Plaur', "Serpine1",'Serpine2')
target_pathway='uPA_sys'

target_genes=c('Thy1', 'Ddr2', "Tcf21", 'Postn', 'S100a4', 'Wt1')
target_pathway='fib_pop'

target_genes=c('Atg3', 'Atg5', "Atg7",'Sqstm1', 'Becn1', "Ulk1")
target_pathway='autophagy'

# -------- plots of interesting genes (n = 4)----------------------------
FeaturePlot(Samples.combined, reduction = "umap", features = target_genes, order = T,
            pt.size = 1, ncol=2, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  feature_theme
ggsave(filename = paste0(path_fig, sample,target_pathway,"_umap_combine.jpeg"), width = 10, height = 12)
# no label
FeaturePlot(Samples.combined, reduction = "umap", features = target_genes, order = T,
            pt.size = 1, ncol=2, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  
ggsave(filename = paste0(path_fig, sample,target_pathway,"_umap_combine.jpeg"), width = 10, height = 12)

FeaturePlot(Samples.combined, reduction = "umap", features = target_genes, order = T, split.by = "treatment",  
            pt.size = 1, ncol=2, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  feature_theme
ggsave(filename = paste0(path_fig,sample,target_pathway,"_umap.jpeg"), width=10 , height = 20)


VlnPlot(object = Samples.combined,  group.by = "treatment", cols = stim.col,
        features = target_genes, pt.size = 0, ncol=2) &
  vinplot_theme & NoLegend()
ggsave(filename = paste0(path_fig, sample,target_pathway,"_Vln.jpeg"), width=10 , height = 10)

VlnPlot(object = Samples.combined,  group.by = "celltypes.short", cols = cluster.col,
        features = target_genes, pt.size = 0.1, ncol=2) &
  vinplot_theme & NoLegend()
ggsave(filename = paste0(path_fig, sample,target_pathway,"_Vln2.jpeg"), width=10 , height = 10)


#-------- plots of interesting genes (n = 6)----------------------------
FeaturePlot(Samples.combined, reduction = "umap", features = target_genes, order = T,
            pt.size = 1, ncol=3, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  feature_theme
ggsave(filename = paste0(path_fig, sample,target_pathway,"_umap_combine.jpeg"), width = 15, height = 10)
# no label
FeaturePlot(Samples.combined, reduction = "umap", features = target_genes, order = T,
            pt.size = 1, ncol=3, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  NoAxes() &
  theme(plot.title = element_blank(), legend.text = element_text(size=15))
ggsave(filename = paste0(path_fig, sample,target_pathway,"_umap_combine_noLabel.jpeg"), width = 15, height = 10)

FeaturePlot(Samples.combined, reduction = "umap", features = target_genes, order = T, split.by = "treatment",  
            pt.size = 1, ncol=2, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  feature_theme
ggsave(filename = paste0(path_fig,sample,target_pathway,"_umap.jpeg"), width=10 , height = 30)


VlnPlot(object = Samples.combined,  group.by = "treatment", cols = stim.col,
        features = target_genes, pt.size = 0, ncol=3) &
  vinplot_theme & NoLegend() & ylab("Отн. уровень экспрессии")
ggsave(filename = paste0(path_fig, sample,target_pathway,"_Vln.jpeg"), width=15 , height = 10)
# no label:
VlnPlot(object = Samples.combined,  group.by = "treatment", cols = stim.col,
        features = target_genes, pt.size = 0, ncol=3) &
  NoLegend() & NoAxes() & theme(plot.title = element_blank())
ggsave(filename = paste0(path_fig, sample,target_pathway,"_Vln_noLabel.jpeg"), width=15 , height = 10)

VlnPlot(object = Samples.combined,  group.by = "celltypes.short", cols = cluster.col,
        features = target_genes, pt.size = 0.1, ncol=3) &
  vinplot_theme & NoLegend() & ylab("Отн. уровень экспрессии")
ggsave(filename = paste0(path_fig, sample,target_pathway,"_Vln2.jpeg"), width=15 , height = 10)
#no label:
VlnPlot(object = Samples.combined,  group.by = "celltypes.short", cols = cluster.col,
        features = target_genes, pt.size = 0.1, ncol=3) &
  NoLegend() & NoAxes() & theme(plot.title = element_blank())
ggsave(filename = paste0(path_fig, sample,target_pathway,"_Vln2_noLabel.jpeg"), width=15 , height = 10)


# dot plot cell type and treatmrnt
DotPlot(Samples.combined, group.by = "celltypes.treat", features = target_genes, 
        dot.scale = 10, assay = "RNA") &
  theme(axis.text.x = element_text(face = "italic"),
        axis.title.x = element_blank(), axis.title.y = element_blank()) &
  scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(path_fig, sample,"_dot_plot.jpeg"), width=8, height = 6)
# No label:
DotPlot(Samples.combined, group.by = "celltypes.treat", features = target_genes, 
        dot.scale = 10, assay = "RNA") &
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_blank()) &
  scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(path_fig, sample,"_dot_plot_noLabel.jpeg"), width=4.5, height = 6)

# dotplot only treatment
DotPlot(Samples.combined, group.by = "treatment", features = target_genes, 
        dot.scale = 8, assay = "RNA", scale = FALSE) &
  coord_flip() &
  theme(axis.text.y = element_text(face = "italic"),
        axis.title.x = element_blank(), axis.title.y = element_blank()) &
  scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(path_fig, sample,"_dot_plot2.jpeg"), width=4, height = 4)
# No label:
DotPlot(Samples.combined, group.by = "treatment", features = target_genes, 
        dot.scale = 8, assay = "RNA", scale = FALSE) &
  coord_flip() &
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_blank()) &
  scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(path_fig, sample,"_dot_plot2_noLabel.jpeg"), width=2.5, height = 4)


# cell percent of autophagy
Cell_Percent <- PrctCellExpringGene(Samples.combined, genes = target_genes, 
                                    group.by = "treatment", assay = "RNA", datatype = "counts", threshold = 0)

rownames(Cell_Percent) <- Cell_Percent$Markers
names(Cell_Percent)
Cell_Percent <- Cell_Percent[-1]
Cell_Percent = as.data.frame(Cell_Percent)
Cell_Percent = melt(t(Cell_Percent))
Cell_Percent$Var1=factor(Cell_Percent$Var1,levels = c("ЛОЖ","КА"))
genes = target_genes
treatment =  c("ЛОЖ","КА")

ggplot(Cell_Percent, aes(x = Var2, y = value, fill = Var1)) + 
  geom_col(position = position_dodge()) + 
  geom_text(aes(label = value), vjust = -0.7, color="black", size=4, position = position_dodge(width = 0.9)) +
  scale_y_continuous(name="Доля клеток, %") +
  scale_fill_manual(values = stim.col)+
  bar_theme +
  theme(axis.text.x = element_text(face = "italic")) +
  geom_hline(yintercept=10, linetype="dashed", color = "black")
ggsave(filename = paste0(path_fig, sample,target_pathway,"_barplot_percentage.jpeg"), width=8 , height = 6)
# no label:
ggplot(Cell_Percent, aes(x = Var2, y = value, fill = Var1)) + 
  geom_col(position = position_dodge()) + 
  scale_fill_manual(values = stim.col)+
  theme_nothing()
ggsave(filename = paste0(path_fig, sample,target_pathway,"_barplot_percentage_noLabel.jpeg"), width=8 , height = 6)


# cell percent of uPA sys populations
target_genes=c('Plau', 'Plaur', "Serpine1",'Serpine2')
target_pathway='uPA_sys'

Cell_Percent <- PrctCellExpringGene(Samples.combined, genes = target_genes, 
                                    group.by = "treatment", assay = "RNA", datatype = "counts", threshold = 0)

rownames(Cell_Percent) <- Cell_Percent$Markers
names(Cell_Percent)
Cell_Percent <- Cell_Percent[-1]
Cell_Percent = as.data.frame(Cell_Percent)
Cell_Percent = melt(t(Cell_Percent))
Cell_Percent$Var1=factor(Cell_Percent$Var1,levels = c("Sham","TAC"))
genes = target_genes
treatment =  c("Sham","TAC")

ggplot(Cell_Percent, aes(x = Var2, y = value, fill = Var1)) + 
  geom_col(position = position_dodge()) + 
  geom_text(aes(label = value), vjust = -0.7, color="black", size=4, position = position_dodge(width = 0.9)) +
  scale_y_continuous(name="Percentage of total cells, %") +
  scale_fill_manual(values = stim.col)+
  bar_theme +
  geom_hline(yintercept=10, linetype="dashed", color = "black")
ggsave(filename = paste0(path_fig, sample,target_pathway,"_barplot_percentage_fibroblasts.jpeg"), width=6 , height = 6)



############## Thy1, Ddr2, Tcf21, Postn, S100a4, Wt1 positive fibroblasts ############
Samples.combined <- readRDS(paste0(path, "Fibroblast_integrated__filtered_processed.rds"))
population <- "Wt1"


FeaturePlot(Samples.combined, reduction = "umap", features = population, order = T, split.by = "treatment", keep.scale = "all",  
            pt.size = 1, ncol=3)  & 
  scale_colour_gradientn(colours = gradient.col) &
  feature_theme
ggsave(filename = paste0(path_fig, population,"before_umap.jpeg"), width=10 , height = 5)


VlnPlot(object = Samples.combined,  group.by = "celltypes.short", split.by = "treatment", assay = "RNA",
        features = population, pt.size = 0.1, cols = stim.col) &
  vinplot_theme & ylab("Отн. уровень экспрессии") &
  geom_hline(size = 1, yintercept=0.3, linetype="dashed", color = "red") # граница отсечки популяции
ggsave(filename = paste0(path_fig, population,"_before_subsetting_Vin.jpeg"), width=7 , height = 5)

VlnPlot(object = Samples.combined, group.by = 'treatment', cols = stim.col, 
        features = population, pt.size = 0.1)&
  vinplot_theme  & ylab("Отн. уровень экспрессии") &
  geom_hline(size = 1, yintercept=0.3, linetype="dashed", color = "red") 
ggsave(filename = paste0(path_fig, population,"_before_subsetting_Vin2.jpeg"), width=6 , height = 5)
# no label
VlnPlot(object = Samples.combined, group.by = 'treatment', cols = stim.col, 
        features = population, pt.size = 0.1)&
  NoLegend() & NoAxes() & theme(plot.title = element_blank()) &
  geom_hline(size = 1, yintercept=0.3, linetype="dashed", color = "red") 
ggsave(filename = paste0(path_fig, population,"_before_subsetting_Vin2_noLabel.jpeg"), width=6 , height = 5)


# Enter in subset your target gene!!!!!:
PosFB <- subset(x = Samples.combined, subset = Wt1 > 0.3)


DefaultAssay(PosFB) = "RNA"
FeaturePlot(PosFB, 
            reduction = "umap", 
            features = population, 
            order = T, split.by = "treatment",  
            pt.size = 1, ncol=2, repel = TRUE, keep.scale = "all") & 
  scale_colour_gradientn(colours = gradient.col) &
  feature_theme
ggsave(filename = paste0(path_fig, population,"_umap_subset.jpeg"), width=10 , height = 5)
# no label:
FeaturePlot(PosFB, 
            reduction = "umap", 
            features = population, 
            order = T, split.by = "treatment",  
            pt.size = 1, ncol=2) &
  scale_colour_gradientn(colours = gradient.col) &
  NoAxes() & theme(plot.title = element_blank(), axis.title.y.right = element_blank())
ggsave(filename = paste0(path_fig, population,"_umap_subset_noLabel.jpeg"), width=10 , height = 5)

VlnPlot(object = PosFB,  group.by = "treatment", cols = stim.col,
        features = population, pt.size = 0) &
  vinplot_theme & NoLegend() & ylab("Отн. уровень экспрессии") 
ggsave(filename = paste0(path_fig, population,"_subset_Vln.jpeg"), width=5 , height = 5)
# no label:
VlnPlot(object = PosFB,  group.by = "treatment", cols = stim.col,
        features = population, pt.size = 0) &
  NoLegend() & NoAxes() & theme(plot.title = element_blank())
ggsave(filename = paste0(path_fig, population,"_subset_Vln_noLabel.jpeg"), width=5 , height = 5)

Idents(PosFB) <- PosFB$treatment

saveRDS(PosFB, file = paste0(path, population,'.RDS'))

PosFB <-  readRDS(file = paste0(path, population,'.RDS'))


# visualize target genes
target_genes=c('Atg3', 'Atg5', "Atg7",'Sqstm1', 'Becn1', "Ulk1")
target_pathway='autophagy'

FeaturePlot(PosFB, 
            reduction = "umap", 
            features = target_genes,
            order = TRUE,
            repel = TRUE, ncol=3, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) &
  feature_theme
ggsave(filename = paste0(path_fig, population, "_", target_pathway,"_features.jpeg"), width=15 , height = 10)
#no label:
FeaturePlot(PosFB, 
            reduction = "umap", 
            features = target_genes,
            order = TRUE,
            repel = TRUE, ncol=3, keep.scale = 'all')  & 
  scale_colour_gradientn(colours = gradient.col) & 
  theme_nothing() &
  theme(legend.position = "right")
ggsave(filename = paste0(path_fig, population, "_", target_pathway,"_features_noLabel.jpeg"), width=15 , height = 10)

VlnPlot(object = PosFB, 
        features = target_genes, ncol=3, cols = stim.col, pt.size = 0)  &
  vinplot_theme & NoLegend() & ylab("Отн. уровень экспрессии") 
ggsave(filename = paste0(path_fig,population, "_", target_pathway,"_vln.jpeg"), width=10 , height = 8)
# no label:
VlnPlot(object = PosFB, 
        features = target_genes, ncol=3, cols = stim.col, pt.size = 0)  &
  NoLegend() & NoAxes() & theme(plot.title = element_blank())
ggsave(filename = paste0(path_fig,population, "_", target_pathway,"_vln_noLabel.jpeg"), width=10 , height = 8)

DotPlot(PosFB, group.by = "treatment", features = target_genes, 
        dot.scale = 8, scale = FALSE, assay = "RNA") &
  coord_flip() &
  theme(axis.text.y = element_text(face = "italic"),
        axis.title.y = element_blank(), axis.title.x = element_blank()) &
  scale_colour_gradientn(colours = gradient.col) &
  ggtitle(paste0(population, "+", " фибробласты"))
ggsave(filename = paste0(path_fig, population,"_", target_pathway,"_dot_plot.jpeg"), width=4, height = 4)
# no label:
DotPlot(PosFB, group.by = "treatment", features = target_genes, 
        dot.scale = 8, scale = FALSE, assay = "RNA") &
  coord_flip() &
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_blank()) &
  scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(path_fig, population,"_", target_pathway,"_dot_plot_noLabel.jpeg"), width=2.5, height = 4)

# Test for DEGs in subsetted fibroblasts ----------------------------------
population <- "Wt1"
PosFB <-  readRDS(file = paste0(path, population,'.RDS'))

TAC_vs_Sham <- FindMarkers(PosFB, ident.1 = 'КА', ident.2 = 'ЛОЖ')
head(TAC_vs_Sham, n=15)


# Reorder columns and sort by padj      
TAC_vs_Sham <- TAC_vs_Sham[, c(1, 3:5,2)]


TAC_vs_Sham <- TAC_vs_Sham %>% 
  dplyr::arrange(p_val_adj) %>%
  mutate(gene = rownames(.), .before = p_val) %>%
  mutate(FC = 2^avg_log2FC,
         FC_sig = case_when(avg_log2FC > 0 ~ FC, # FC=FoldCange
                            avg_log2FC < 0 ~ -(1/FC))) 

write_tsv(TAC_vs_Sham, file = paste0(path, population, "_TAC_vs_Sham_markers.tsv"), col_names = TRUE)


top10 <- TAC_vs_Sham %>% 
  top_n(n = 10, 
        wt = avg_log2FC)
top10




# cell percent of target genes --------------------------------------------
population <- "Wt1"
PosFB <-  readRDS(file = paste0(path, population,'.RDS'))

Cell_Percent <- PrctCellExpringGene(PosFB, genes = target_genes, 
                                    group.by = "treatment", assay = "RNA", datatype = "counts", threshold = 0.3)

rownames(Cell_Percent) <- Cell_Percent$Markers
names(Cell_Percent)
Cell_Percent <- Cell_Percent[-1]
Cell_Percent = as.data.frame(Cell_Percent)
Cell_Percent = melt(t(Cell_Percent))
Cell_Percent$Var1=factor(Cell_Percent$Var1,levels = c("ЛОЖ","КА"))
genes = target_genes
treatment =  c("ЛОЖ","КА")
Cell_Percent

ggplot(Cell_Percent, aes(x = Var2, y = value, fill = Var1)) + 
  geom_col(position = position_dodge()) + 
  geom_text(aes(label = value), vjust = -0.7, color="black", size=4, position = position_dodge(width = 0.9)) +
  scale_y_continuous(name="Доля клеток, %") +
  scale_fill_manual(values = stim.col)+
  bar_theme  +
  ggtitle(paste0(population, "+", " фибробласты"))
ggsave(filename = paste0(path_fig, population, "_",target_pathway,"_barplot_percentage.jpeg"), width=8 , height = 6)
# no label:
ggplot(Cell_Percent, aes(x = Var2, y = value, fill = Var1)) + 
  geom_col(position = position_dodge()) + 
  scale_y_continuous(name="Доля клеток, %") +
  scale_fill_manual(values = stim.col)+
  theme_nothing() +
  theme(legend.position = "right") +
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(), axis.ticks.y = element_blank(), plot.caption = element_text(size = 20, hjust = 0), legend.text = element_text(size = 20)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.text = element_blank(), legend.title = element_blank())
ggsave(filename = paste0(path_fig, population, "_",target_pathway,"_barplot_percentage_noLabel.jpeg"), width=8 , height = 6)

