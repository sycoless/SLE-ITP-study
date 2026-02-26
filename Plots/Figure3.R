library(Seurat)
library(patchwork)
library(cowplot)
library(tidyverse)
library(enrichR)
library(reshape2)
options(stringsAsFactors = F)



# ---------------------------------------
#       Figure 3A
# ---------------------------------------
setwd("./RNA/BMMC/")
# data <- read.csv("./tables/BMMC_inflamm.score.csv")
data <- read.csv("./tables/BMMC_inflamm.score.csv") %>% 
  mutate(major_cell_type_r1 = ifelse(cell_type_r1 %in% c('cM_CD36',
                      'cM_TMSB10',
                      'cM_NFKBIA_IL1B',
                      'cM_HLA-DR_CD74'),"cM",ifelse(cell_type_r1 %in% c('ncM_AIF1',
                      'Mono_C1Q'),"ncM",major_cell_type_r1)))
colnames(data)
gene_plot.df <- data.frame(
  inflammatory = data$GOBP_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS_gene_signature_ucell_score,
  Group = data$major_cell_type_r1) %>% melt
unique(gene_plot.df$Group)
gene_plot.df_ <- gene_plot.df %>% 
  filter(Group %in% c('Neutrophil','cM','ncM','EBM','DC','pDC'))

gene_plot.df_$Group <- factor(gene_plot.df_$Group ,levels = c('Neutrophil','cM','ncM','EBM','DC','pDC'))
ggplot(gene_plot.df_, aes(x = Group, y = value)) +
  geom_violin(aes(color = Group, fill = Group), scale = "width") +
  theme_cowplot(font_size = 12) +
  labs(y = "UCell score",title="Infammatory response \nto antigenic stimulus") +
  scale_fill_manual(values=c("#8fc0dc","#edeaa4","#bcacd3","#4e9592","#dbad5f","#ac5092"))+
  scale_color_manual(values=c("#8fc0dc","#edeaa4","#bcacd3","#4e9592","#dbad5f","#ac5092" ))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(size=12,hjust=0.5,face="plain"),
        strip.background = element_blank(),
        legend.position = "null")
ggsave("./figures/BMMC_myeloid_inflamm-2026.pdf",width = 4.2,height = 2.5)


# ---------------------------------------
#       Figure 3C
# ---------------------------------------
DEG_dir <- "./RNA/BMMC/subclustering/Myeloid/markers/Mono/"
DEG_tests <- dir(DEG_dir) 
options(stringsAsFactors = F)
library(tidyverse)
library(enrichR)
combined <- Reduce(rbind, lapply(dir(DEG_dir), function(file){
  read.csv(paste0(DEG_dir, file))
}))

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','MSigDB_Hallmark_2020','Reactome_Pathways_2024', 'WikiPathway_2021_Human', 'KEGG_2021_Human') #,'ENCODE_TF_ChIP-seq_2015'
combined$group <- factor(combined$group)
input_list <- combined %>% 
  group_by(group) %>% 
  top_n(n = 200, wt = avg_log2FC) 
input_list <- split(input_list$gene, input_list$group)

# size of lists
lapply(input_list, function(x){
  print(length(x))
})

enriched_df <- do.call(rbind, lapply(names(input_list), function(x){
  print(x)
  if(length(input_list[[x]]) > 0){
    cur_enrich <- enrichr(input_list[[x]], dbs)
  } else{return(data.frame())}
  cur_df <- do.call(rbind, lapply(dbs, function(cur_db){
    df <- cur_enrich[[cur_db]]
    if(nrow(df) > 1){df$degs <- x;  df$db <- cur_db}
    else{df <- data.frame()}
    df
  }))
  
}))


combined_output2 <- data.frame()
combined_output2 <- rbind(combined_output2, enriched_df) 
res_out <- combined_output2%>% 
  filter(Adjusted.P.value < 0.05) 

description_used <- c(
  "Neutrophil Degranulation",
  "GTPase Activator Activity (GO:0005096)",
  "Toll-like Receptor Cascades",
  "Specific Granule (GO:0042581)",
  
  "Defense Response To Fungus (GO:0050832)",
  "Transition Metal Ion Binding (GO:0046914)",
  "Arachidonic Acid Binding (GO:0050544)",
  
  "NF-kappa B signaling pathway",
  "IL-18 signaling pathway WP4754",
  "Interleukin-10 Signaling",
  "Osteoclast differentiation",
  "Interleukin-1 Signaling",
  "Regulation Of Angiogenesis (GO:0045765)",
  "Positive Regulation Of Angiogenesis (GO:0045766)",
  
  "Interferon Alpha Response",
  "Interferon Gamma Signaling",
  "MHC Protein Complex (GO:0042611)",
  "Regulation Of T Cell Activation (GO:0050863)",
  "Th17 cell differentiation",
  
  "Platelet Activation, Signaling and Aggregation",
  "Regulation Of Dendritic Cell Differentiation (GO:2001198)",
  "Positive Regulation Of Phagocytosis (GO:0050766)",
  "Regulation Of Macrophage Activation (GO:0043030)",
  "B Cell Mediated Immunity (GO:0019724)"
)

description_used <- c(
  "Neutrophil Degranulation",
  "Specific Granule (GO:0042581)",
  "Secretory Granule Membrane (GO:0030667)",
  "Toll-like Receptor Cascades",
  
  "Defense Response To Fungus (GO:0050832)",
  "Cellular Responses to Stress",
  "Oxidative phosphorylation",
  "Thermogenesis",
  
  "Inflammatory Response",
  "NF-kappa B signaling pathway",
  "Interleukin-10 Signaling",
  "Glucocorticoid Receptor Pathway WP2880",
  "Osteoclast differentiation",
  "IL1 and megakaryocytes in obesity WP2865",
  "Regulation Of Angiogenesis (GO:0045765)",
  "Positive Regulation Of Angiogenesis (GO:0045766)",
  
  
  "Antigen processing and presentation",
  "Interferon Signaling",
  "MHC Class II Protein Complex (GO:0042613)",
  "Interferon Alpha Response",
  "Th17 cell differentiation",
  
  "Platelet Activation, Signaling and Aggregation",
  "Platelet Alpha Granule (GO:0031091)",
  
  "Positive Regulation Of Phagocytosis (GO:0050766)",
  "Regulation Of Macrophage Activation (GO:0043030)",
  "Regulation Of Phagocytosis (GO:0050764)"
)
library(Matrix)
{
  cluster_go_plot <- res_out %>% 
    filter(Term %in% description_used)
  dim_x <- unique(as.character(cluster_go_plot$degs))
  dim_y <- description_used
  map_x <- setNames(seq_along(dim_x), dim_x)
  map_y <- setNames(seq_along(dim_y), dim_y)
  cluster_go_plot_mat <- sparseMatrix(
    i=map_x[as.character(cluster_go_plot$degs)], 
    j=map_y[as.character(cluster_go_plot$Term)], 
    x=-log10(cluster_go_plot$Adjusted.P.value), 
    dims=c(length(dim_x), length(dim_y)), 
    dimnames=list(dim_x, dim_y)
  )
  cluster_go_plot_mat <- as.matrix(cluster_go_plot_mat)
  cluster_go_plot_mat <- cluster_go_plot_mat[c('cM_CD36',
                                               'cM_TMSB10',
                                               # 'cM_NFKB2A_IL1B',
                                               'cM_NFKBIA_IL1B',
                                               'cM_HLA-DR_CD74',
                                               'ncM_AIF1',
                                               'Mono_C1Q'),]
  cluster_go_plot_mat_quantile <- quantile(cluster_go_plot_mat, c(0.1, 0.9))
  cluster_go_plot_mat <- pmax(cluster_go_plot_mat, cluster_go_plot_mat_quantile[1])
  cluster_go_plot_mat <- pmin(cluster_go_plot_mat, cluster_go_plot_mat_quantile[2])
  color_used <- circlize::colorRamp2(seq(min(cluster_go_plot_mat), max(cluster_go_plot_mat), length = 9), rev(RColorBrewer::brewer.pal(11,"RdYlBu")[2:10]))
  library(ComplexHeatmap)
  colnames(cluster_go_plot_mat) <- str_replace(colnames(cluster_go_plot_mat), " \\s*\\([^\\)]+\\)", "")
  colnames(cluster_go_plot_mat) <- str_replace(colnames(cluster_go_plot_mat), " \\s*WP[0-9]+", "")
  p <- Heatmap(cluster_go_plot_mat %>% t(),
               cluster_rows = F,
               cluster_columns = F,
               show_row_names = T,
               column_names_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               col = color_used,
               rect_gp = gpar(col = "white", lwd = 0.1),
               name = "-log10(adjusted P-value)",
               show_heatmap_legend=FALSE,
               row_names_max_width = unit(10, "cm"))
  p
}

lgd = Legend(col_fun = color_used, title = "-log10(adjusted P-value)" ,grid_width = unit(0.3, "cm"),title_gp = gpar(fontsize=10),labels_gp = gpar(fontsize=10))

pdf("./RNA/BMMC/subclustering/Myeloid/figures/Myeloid.marker.enrich.plot-onlyuseMono.pdf",height = 7,width = 6.2)
p
dev.off()
dev.off()

# ---------------------------------------
#       Figure 3D
# ---------------------------------------
setwd("./RNA/BMMC/subclustering/Myeloid/")
data <- read.csv("./tables/myeloid.signature.csv")
colnames(data)

gene_plot.df <- data.frame(
  Macrophage = data$Macro_gene_signature_ucell_score,
  Phagocytosis = data$Phagocytosis_gene_signature_ucell_score,
#   Angiogenesis = data$Angiogenesis_gene_signature_ucell_score,
  Cytokine = data$GOMF_CYTOKINE_ACTIVITY_gene_signature_ucell_score,
  Group = data$cell_type_r1) %>% melt
unique(gene_plot.df$Group)
gene_plot.df_ <- gene_plot.df %>% 
  filter(Group %in% c('cM_CD36',
                      'cM_TMSB10',
                      'cM_NFKBIA_IL1B',
                      'cM_HLA-DR_CD74',
                      'ncM_AIF1',
                      'Mono_C1Q'))

gene_plot.df_$Group <- factor(gene_plot.df_$Group ,levels = c('cM_CD36',
                                          'cM_TMSB10',
                                          'cM_NFKBIA_IL1B',
                                          'cM_HLA-DR_CD74',
                                          'ncM_AIF1',
                                          'Mono_C1Q'))
ggplot(gene_plot.df_, aes(x = Group, y = value)) +
  geom_violin(aes(color = Group, fill = Group), scale = "width") +
  labs(y = "UCell score") +
  theme_cowplot(font_size = 12) +
  scale_fill_manual(values=c("#82c785","#edeaa4","#cdaa9f","#794976","#bcacd3","#889b5d"  ))+
  scale_color_manual(values=c("#82c785","#edeaa4","#cdaa9f","#794976","#bcacd3","#889b5d"  ))+
  # facet_grid(variable ~ .,scales="free",) +
  facet_wrap(variable ~ .,
             scales="free_y",
             strip.position="top",ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中
        # plot.title.position = "plot",  # 标题位置相对于整个图表
        strip.background = element_blank(),
        legend.position = "null")
ggsave("./figures/C1Q_marophage-2026.pdf",width = 5,height = 4.8)


# ---------------------------------------
#       Figure 3E
# ---------------------------------------
{
  dat <- read.csv("./Myeloid_cellprop.csv") %>%
  column_to_rownames("X") %>% 
  mutate(Sample=rownames(.)) %>% 
  mutate(Disease=sub("_\\d+$", "", rownames(.))) %>% 
  as.data.frame()

# 
df_long <- dat %>%
  pivot_longer(
    cols = -c(Sample, Disease),   # 除了 Sample 和 Group 外的所有列
    names_to = "CellType",      # 新列：原来的列名
    values_to = "Value"         # 新列：原来的数值
  )

unique(df_long$CellType)
my_comparisons <- list(c("HC", "SLE-NP"),c("SLE-ITP", "SLE-NP"),c("SLE-ITP", "HC"))
df_long$Disease <- factor(df_long$Disease,levels=c("HC","SLE-NP","SLE-ITP"))
library(ggpubr)
}


for(i in unique(df_long$CellType)){
  print(i)
  ggplot(df_long %>% filter(CellType==i), aes(x = Disease, y = Value, fill = Disease)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
    geom_point(position = position_dodge(0.8), size = 2, shape = 21) +
    stat_compare_means(
      aes(group = Disease),
      comparisons=my_comparisons,
      method = "wilcox.test", label = "p.signif") +
    # method = "t.test", label = "p.format") +
    # theme_minimal() +
    theme_classic()+
    labs(title = i, y = "% of myeloid cells", x = "",fill="") +
    # scale_fill_manual(values = c("CR" = "#4B98C5", "PR" = "#BB2B34")) +
    scale_fill_manual(values = c("HC"="#B2CFE5","SLE-NP"="#D97351", "SLE-ITP"="#8A2A46" )) +
    scale_y_continuous(expand=c(0,2))+
    # guide_legend(fill=())
    theme(
      # axis.text.x = element_text(angle = 45, hjust = 1)
      # axis.text.x = element_text(angle=45,hjust=1,size=10,color = "black"),
      plot.title = element_text( hjust = 0.5,size=12),
      axis.text.x = element_text(size=12,color = "black"),
      axis.text.y = element_text(size=12,color = "black"),
      legend.position="none"
    )
  
  ggsave(sprintf("./RNA/BMMC/figures/Myeloid_%s_t.test.pdf",i),width = 2.5,height = 3.5)
}


# ---------------------------------------
#       Figure 3F
# ---------------------------------------

library(SCopeLoomR)
library(data.table)
library(ggplot2)
library(dplyr)
# library(reshape)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(data.table)
library(repr)
library(tidyverse)
options(stringsAsFactors = F)

{
  data2 <- data %>% 
    filter(cell_type_r1 %in% c('cM_CD36',
                               'cM_TMSB10',
                               'cM_NFKBIA_IL1B',
                               'cM_HLA-DR_CD74',
                               'ncM_AIF1',
                               'Mono_C1Q')) %>%

    filter(!(Sample %in% c("SLE-ITP_5P","SLE-ITP_6P","SLE-ITP_8P","SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P")))
  data2$Disease_ct <- paste(data2$cell_type_r1,data2$Disease,sep="_")
  data2$Sample_ct <- paste(data2$cell_type_r1,data2$Sample,sep="_")
  table(data2$Sample)
  data2$group1_ct <- paste(data2$cell_type_r1,data2$group1,sep="_")
  gene_plot.df <- data.frame(
    M1 = data2$M1.blood_gene_signature_ucell_score,
    M2 = data2$M2.blood_gene_signature_ucell_score,
    M1.M2 = (data2$M1.blood_gene_signature_ucell_score/data2$M2.blood_gene_signature_ucell_score),
    Group = data2$Disease_ct
  )
  
  dat <- gene_plot.df %>% group_by(Group) %>%
    summarize(across(where(is.numeric),mean)) %>%
    column_to_rownames("Group") 
  df <- t(scale(dat))
  
  library(scales)
  
  df <- dat
  df$M1 <- rescale(df$M1, to=c(0, 2))
  df$M2 <- rescale(df$M2, to=c(0, 2))
  df$M1.M2 <- rescale(df$M1.M2, to=c(0, 2))

  col_fun = colorRamp2(c(min(df),0,max(df)), c("#4B98C5",  "white","#BB2B34"))
  ha = HeatmapAnnotation(Disease = factor(c(rep(c("HC","SLE-NP","SLE-ITP"),1)),
                                          levels = c("HC","SLE-NP","SLE-ITP")), 
                         cell_type = factor(rep(c('Mono_C1Q'),each=3),levels = c('Mono_C1Q')),
                         annotation_name_side = "left",
                         annotation_legend_param = list(
                           Disease = list(direction = "horizontal"),
                           cell_type = list(nrow = 3)),
                         col = list( Disease=c("HC"="#B2CFE5","SLE-ITP"="#8A2A46","SLE-NP"="#D97351" ),
                         cell_type=c('Mono_C1Q'="#889b5d" )))
  colnames(df)[3] <- "M1/M2"
  ht <- Heatmap(t(df[paste(rep(c(
      'Mono_C1Q'),each=3),c('HC','SLE-NP','SLE-ITP'),sep="_"),]),name=" ",
                col = col_fun,
                # col = coul,
                border_gp = gpar(col = "black", lty = 1),
                rect_gp = gpar(col = "white", lwd = 0.2),
                cluster_rows = F,
                # cluster_rows = as.dendrogram(row_clust),
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                cluster_columns = F,
                # column_title = ,
                column_dend_side="bottom",
                column_dend_height = unit(8, "mm"),
                row_dend_width = unit(8, "mm"),
                show_column_names = FALSE,
                column_title_gp = gpar(fontface = "bold",fontsize = 14),
                ht_opt(heatmap_column_names_gp=gpar(fontsize = 12), #fontface = "bold",
                       heatmap_row_names_gp = gpar(fontsize = 14)), #fontface = "bold"
                row_names_side = "left",
                top_annotation = ha,
                heatmap_legend_param = (list(at = c(0,1, 2),
                                             labels = c("Min","","Max")     
                ))
  )
  # pdf("./RNA/BMMC/subclustering/Myeloid/figures/M1.M2.signature2.pdf",width = 3,height = 2.2)
  draw(ht,
       heatmap_legend_side="right",
       annotation_legend_side="bottom",
       # legend_grouping = "adjusted"
  )
  dev.off()
}

# ---------------------------------------
#       Figure 3H
# ---------------------------------------
library(enrichR)

res2 <- read.csv("./RNA/BMMC/tables/DEGs/MAST/fine/SLE-NP_ct/SLE-NP_ct_Mono_C1Q.csv")
res2 <- read.csv("./RNA/BMMC/tables/DEGs/MAST/fine/SLE-ITP_NP/SLE-ITP_ct_Mono_C1Q.csv")


dbs <- c('GO_Biological_Process_2021', 'WikiPathway_2021_Human', 'KEGG_2021_Human','Reactome_Pathways_2024','MSigDB_Hallmark_2020')


input_list <- list(
  up = res2 %>% filter(p_val<0.01 & avg_log2FC > 0) %>% .$gene,
  down = res2 %>% filter(p_val<0.01 & avg_log2FC < 0) %>% .$gene
)

enriched_df <- do.call(rbind, lapply(names(input_list), function(x){
  if(length(input_list[[x]]) > 0){
    cur_enrich <- enrichr(input_list[[x]], dbs)
  } else{return(data.frame())}
  cur_df <- do.call(rbind, lapply(dbs, function(cur_db){
    df <- cur_enrich[[cur_db]]
    if(nrow(df) > 1){df$degs <- x;  df$db <- cur_db}
    else{df <- data.frame()}
    df
  }))
}))


combined_output2 <- data.frame()
combined_output2 <- rbind(combined_output2, enriched_df) %>% 
  filter(Adjusted.P.value < 0.1) 
saveRDS(combined_output2,"./RNA/BMMC/tables/c1q.pathway_sletpct2.rds")
saveRDS(combined_output2,"./RNA/BMMC/tables/c1q.pathway_sletp-np.rds")

df1 <- readRDS("./RNA/BMMC/tables/c1q.pathway_sletpct2.rds") %>% mutate(group="vs. HC")
df2 <- readRDS("./RNA/BMMC/tables/c1q.pathway_sletp-np.rds") %>% mutate(group="vs. SLE-NP")
colnames(df1)
df <- rbind(df1,df2) %>% filter(degs=="up")
library(viridis)
unique(df$db)
library(Seurat)
selected.path <- c(
  "cellular response to type I interferon (GO:0071357)",
  "cellular response to interferon-gamma (GO:0071346)",
  "neutrophil degranulation (GO:0043312)",
  "interleukin-27-mediated signaling pathway (GO:0070106)",
  "regulation of tumor necrosis factor production (GO:0032680)",
  "regulation of interleukin-12 production (GO:0032655)",
  "regulation of interleukin-6 production (GO:0032675)",
  "Fc receptor signaling pathway (GO:0038093)",
  "regulation of stem cell differentiation (GO:2000736)"
)

wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
plot.df <- df %>% filter(Term %in% c(selected.path))
plot.df$Term <- str_replace(plot.df$Term, " \\s*\\([^\\)]+\\)", "")
plot.df$wrap <- wrapText(plot.df$Term, 50)
plot.df$wrap <- factor(plot.df$wrap,levels = unique(plot.df$wrap))
pdf("./RNA/BMMC/figures/Mono_C1Q_go_enrichment.pdf",width=5.5,height=4)
plot.df %>%
  filter(Adjusted.P.value< 0.1)%>%
  ggplot(aes(x = group, y = wrap, color = -log(Adjusted.P.value), size=log(Odds.Ratio))) +
  geom_point() +
  theme_bw()+
  scale_color_stepsn(colors=rev(magma(256))) +
  labs(title="SLE-ITP Mono_C1Q\nUp-regulated")+
  # scale_color_stepsn(colors=rev(mako(256))) +
  RotatedAxis() + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust=0.5,size=10,color='black'),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(size=1, color='black', fill=NA),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1,color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0.2,0.2,0.2,0.2),
    panel.grid = element_line(size=0.25, color='lightgrey'),
  )
dev.off()

