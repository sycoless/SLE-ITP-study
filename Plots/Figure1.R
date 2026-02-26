options(stringsAsFactors = F)
library(tidyverse)
library(cowplot)

# ---------------------------------------
#       Figure 1B
# ---------------------------------------
data <- read.csv("./data/clinical/clinical_sample.csv") %>% 
  filter(!Sample %in% c("HC_1","HC_2","HC_3","HC_4",
                          "ITP_1", "ITP_2","ITP_3","ITP_4"))


data$Sample <- factor(data$Sample,levels = c("HC_5" ,  "HC_6"  ,"HC_7" ,
                                             "SLE-NP_1"  , "SLE-NP_2" ,  "SLE-NP_3"  ,
                                             "SLE-ITP_1" , "SLE-ITP_2" , "SLE-ITP_3" , 
                                             "SLE-ITP_4" , "SLE-ITP_5" , "SLE-ITP_6",
                                             "SLE-ITP_7" , "SLE-ITP_8" , "SLE-ITP_9"  ,
                                             "SLE-ITP_10", "SLE-ITP_11" ,"SLE-ITP_12",
                                             "SLE-ITP_13" ,"SLE-ITP_14"))

data[,11:18] <- lapply(data[,11:18], function(x) {
  x[is.na(x) | x == ""] <- 'No Sample'
  x[x == "Yes"] <- 'No CAR-T Treatment'
  return(x)
})

data[,4:5] <- lapply(data[,4:5], function(x) {
  x[is.na(x) | x == ""] <- 'Not applicable'
  return(x)
})
data

{
    group_info_color <- 
    list(
      # "HC"="#B2CFE5","ITP"="#D2D1D2","SLE-NP"="#D97351","SLE-ITP"="#8A2A46",
      "HC"="#B2CFE5","ITP"="#9E6D7F","SLE-NP"="#D97351" ,"SLE-ITP"="#8A2A46",
      "Yes" = "#5384B6", 
      "No Sample"= "gray90",  #"#D3D3D3",
      # "Dataset" = "",
      # "Responder"="#417D8B",
      # "Non-responder"="#F2D2D5",
      "CR"="#417D8B",
      "PR"="#F2D2D5",
      "Not applicable" = "gray90",  # "#D3D3D3",
      "Not available"="gray90"     #"#D3D3D3"
      
    )
data$Sample
  group_info <- data %>% 
    select(Sample,Group,CAR.T_Response) %>% 
    pivot_longer(!Sample, names_to = "group_info_col", values_to = "value")
  # CAR.T_Tx
  group_info$value <- factor(group_info$value,levels = names(group_info_color))
  group_info$group_info_col <- str_replace(group_info$group_info_col,"CAR.T_Response","CAR-T Response")
  group_info$group_info_col <- factor(group_info$group_info_col,levels = c("CAR-T Response","Group"))
  clinical_group <-
    ggplot(group_info,aes(x=Sample,y=group_info_col,fill= value))+
    geom_tile(color = "white",lwd=0.5)+
    scale_fill_manual(values=unlist(group_info_color),na.value ="white") +
    #scale_fill_gradient(low="white",high = "black")+
    # scale_fill_manual(values=c("F"=alpha("#FF0000",alpha = 1),"M"=alpha("#5CACEE",alpha = 1)))+
    guides(fill=guide_legend(ncol = 1)) +
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,-5,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      # axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      # axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.text.y = element_text( color="black", size=8, angle = 0,hjust=1,vjust = 0.5), #face="bold",
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  clinical_group
}

{
    Baseline_charac_color <-
    list(
      "Positive" = "#5384B6",
      "Negative"="#FDF5EE",
      "Not available"= "gray90",#"#D3D3D3",
      "anti-GPIIb/IIIa"="#BB3E45"
      
    )
  colnames(data)
  Baseline_charac <- data %>%
    select(Sample,Anti.platelet.antibody,Anti.cMpl.antibody) %>%
    pivot_longer(!Sample, names_to = "group_info_col", values_to = "value")
  Baseline_charac$Sample <- factor(Baseline_charac$Sample,levels = c("HC_5" ,  "HC_6"  ,"HC_7" ,
                                                                     "SLE-NP_1"  , "SLE-NP_2" ,  "SLE-NP_3"  ,
                                                                     "SLE-ITP_1" , "SLE-ITP_2" , "SLE-ITP_3" , 
                                                                     "SLE-ITP_4" , "SLE-ITP_5" , "SLE-ITP_6",
                                                                     "SLE-ITP_7" , "SLE-ITP_8" , "SLE-ITP_9"  ,
                                                                     "SLE-ITP_10", "SLE-ITP_11" ,"SLE-ITP_12",
                                                                     "SLE-ITP_13" ,"SLE-ITP_14"))
  Baseline_charac$group_info_col <- str_replace(Baseline_charac$group_info_col,"Anti.platelet.antibody","Anti-platelet Ab")
  Baseline_charac$group_info_col <- str_replace(Baseline_charac$group_info_col,"Anti.cMpl.antibody","Anti-cMpl Ab")
  Baseline_charac$group_info_col <- factor(Baseline_charac$group_info_col,levels = rev(c("Anti-platelet Ab","Anti-cMpl Ab")))
  Baseline_charac$value <- factor(Baseline_charac$value,levels = names(Baseline_charac_color))
  clinical_Baseline_charac <-
    ggplot(Baseline_charac,aes(x=Sample,y=group_info_col,fill= value))+
    geom_tile(color = "white",lwd=0.5)+
    scale_fill_manual(values=unlist(Baseline_charac_color)) +
    #scale_fill_gradient(low="white",high = "black")+
    # scale_fill_manual(values=c("F"=alpha("#FF0000",alpha = 1),"M"=alpha("#5CACEE",alpha = 1)))+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(0,-5,0,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      #legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x= element_blank(),
      #axis.text.x = element_text( color="black", size=6, angle = 90,hjust = 1,vjust = 1),
      # axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.text.y = element_text( color="black", size=8, angle = 0,hjust=1,vjust = 0.5), #face="bold",
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  clinical_Baseline_charac
}

{
    sample_color <- 
    list(
      "Yes" = "#5384B6",
      "Not available" =  "gray90",
      "No Sample" = "gray90",
      "No CAR-T Treatment" = "#B8D5E9",
      "Pre CAR-T"= "#87CEFA",
      "Post CAR-T"= "#7587b1",
      "Pre/Post CAR-T"="#37538B"
      
    )
  colnames(data)
  sample <- data %>% 
    # select(Sample,WGS,HSPC.scRNA.seq,BMMC.scRNA.seq,BMMC.scBCR.seq,BMMC.scTCR.seq,BMMC.scATAC.seq,PBMC.scRNA.seq,PBMC.scTCR.seq,PBMC.scBCR.seq) %>% 
    # select(Sample,HSPC.scRNA.seq,BMMC.scRNA.seq,BMMC.scBCR.seq,BMMC.scTCR.seq,BMMC.scATAC.seq,PBMC.scRNA.seq,PBMC.scTCR.seq,PBMC.scBCR.seq) %>%
    select(Sample,HSPC.scRNA.seq,BMMC.scRNA.seq,BMMC.scBCR.seq,BMMC.scATAC.seq,PBMC.scRNA.seq,PBMC.scBCR.seq) %>% 
    pivot_longer(!Sample, names_to = "group_info_col", values_to = "value")
  unique(sample$value)
  sample$Sample <- factor(sample$Sample,levels = c("HC_5" ,  "HC_6"  ,"HC_7" ,
                                                   "SLE-NP_1"  , "SLE-NP_2" ,  "SLE-NP_3"  ,
                                                   "SLE-ITP_1" , "SLE-ITP_2" , "SLE-ITP_3" , 
                                                   "SLE-ITP_4" , "SLE-ITP_5" , "SLE-ITP_6",
                                                   "SLE-ITP_7" , "SLE-ITP_8" , "SLE-ITP_9"  ,
                                                   "SLE-ITP_10", "SLE-ITP_11" ,"SLE-ITP_12",
                                                   "SLE-ITP_13" ,"SLE-ITP_14"))
  sample$group_info_col <- factor(sample$group_info_col,levels = rev(c("WGS","HSPC.scRNA.seq","BMMC.scRNA.seq","BMMC.scBCR.seq",
                                                                       "BMMC.scTCR.seq","BMMC.scATAC.seq","PBMC.scRNA.seq",
                                                                       "PBMC.scTCR.seq","PBMC.scBCR.seq")))
  sample$value <- factor(sample$value,levels = names(sample_color))
  
  sample$group_info_col_display <- sample$group_info_col
  sample$group_info_col_display <- str_replace(sample$group_info_col_display, "^HSPC\\.", "")
  sample$group_info_col_display <- str_replace(sample$group_info_col_display, "^BMMC\\.", "")
  sample$group_info_col_display <- str_replace(sample$group_info_col_display, "^PBMC\\.", "")
  sample$group_info_col_group <- case_when(
    grepl("^HSPC", sample$group_info_col) ~ "HSPC",
    grepl("^BMMC", sample$group_info_col) ~ "BMMC",
    grepl("^PBMC", sample$group_info_col) ~ "PBMC",
    TRUE ~ sample$group_info_col
  )
  groups <- unique(sample$group_info_col_group)
  heatmaps <- lapply(groups, function(group) {
    # 过滤当前组的数据
    group_data <- sample %>% filter(group_info_col_group == group)
    
    # 绘制该组的热图
    ggplot(group_data, aes(x = Sample, y = group_info_col_display, fill = value)) +
      geom_tile(color = "white", lwd = 0.5, width = 1, height = 1) + # 设置单元格的宽度和高度为一致
      scale_fill_manual(values = unlist(sample_color), na.value = "white") +
      theme_bw() +
      theme(
        plot.margin = unit(c(0, -5, 0, 10), units = "pt"),
        #legend.position = "none", # 去除图例
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 8, angle = 0, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing = unit(1, "lines") # 行间隔
      ) +
      ggtitle(group) # 为每个热图加上标题
  })
  
  # combined_plot <- plot_grid(plotlist = heatmaps, ncol = 1, rel_heights = c(1,1,4,3), align = "hv")
  combined_plot <- plot_grid(plotlist = heatmaps, ncol = 1, rel_heights = c(1,4,3), align = "hv")
  combined_plot
}


  overview <- plot_grid(
    ncol = 1, align = "hv",
    clinical_group +  theme(plot.title = element_text(hjust = 1,size=6),axis.title.y = element_blank(),
                            axis.text = element_text(hjust = 1,size=6),plot.margin = unit(x=c(0,10,0,10),units="pt"),
                            axis.title.x = element_blank(),legend.position = "none"),
    clinical_Baseline_charac + theme(plot.title = element_text(hjust = 1,size=6),axis.title.y = element_blank(),
                                     axis.text = element_text(hjust = 1,size=6),plot.margin = unit(x=c(0,10,0,10),units="pt"),
                                     axis.title.x = element_blank(),legend.position = "none"),
                         # axis.title.x = element_blank(),legend.position = "none"),
    heatmaps[[1]] + theme(plot.title = element_blank(),axis.title.y = element_blank(),
                          axis.text = element_text(hjust = 1,size=6),plot.margin = unit(x=c(0,10,0,10),units="pt"),
                          axis.title.x = element_blank(),legend.position = "none"),
    heatmaps[[2]] + theme(plot.title = element_blank(),axis.title.y = element_blank(),
                          # panel.border = element_rect(color = "#CC0000", fill = NA, size = 1),
                          axis.text = element_text(hjust = 1,size=6),plot.margin = unit(x=c(0,10,0,10),units="pt"),
                          axis.title.x = element_blank(),legend.position = "none"),
    heatmaps[[3]] + theme(plot.title = element_blank(),axis.title.y = element_blank(),
                          # panel.border = element_rect(color = "#1f77b4", fill = NA, size = 1),
                          axis.text = element_text(hjust = 1,size=6),plot.margin = unit(x=c(0,10,0,10),units="pt"),
                          axis.title.x = element_blank(),legend.position = "none"),
    rel_heights = c(2,1.5,1.2,3,2)
    # rel_heights = c(2,1.5,1,1,3.3,2.2)
  )
overview
plot_dir <- "./Final_figures/"
overview
ggsave(paste0(plot_dir,"Fig1B.pdf"),overview,width = 6,height = 3)
legend
ggsave(paste0(plot_dir,"Fig1B-legend.pdf"),legend,width = 6,height = 5)



# ---------------------------------------
#       Figure 1C-E
# ---------------------------------------
ibrary(tidyverse)
scriptPath <- "./ATAC/utils/"
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/plotting_config.R"))
pointSize <- 0.005


disease_camp <- list("HC"="#B2CFE5","ITP"="#D2D1D2","SLE-ITP"="#8A2A46","SLE-NP"="#D97351")
disease_cmap <- unlist(disease_camp)

setwd("./RNA/")

# Figure 1C
# HSPC RNA labels UMAP
umapDF <- read.csv("./HSPC/tables/HSPC_umapDF.csv")
hspc_label <- read.csv("./HSPC/tables/HSPC_celltype_labels.csv") %>% 
  dplyr::select(c("barcode_id","cell_type"))
Cluster_levels <- c('HSC','MPP','MPP_cycling','LMPP','MPP_myeloid','GMP','pDCP','DCP','MDP','CLP','ProB','PreB','MEMP','EBMP','MKP','EryP')

umapDF <- umapDF %>% 
  left_join(.,hspc_label,by="barcode_id") 
colnames(umapDF) <- c("barcode_id","UMAP1","UMAP2","NamedClust")
rownames(umapDF) <- umapDF$barcode_id
umapDF <- umapDF[,-1]
umapDF$NamedClust <- factor(umapDF$NamedClust,levels = Cluster_levels)
# hspcCmap <- c('#CE5B5B','#FF9D9A','#935e5e','#79706E',"#AA9BC5",'#D68645','#F28E2B','#FFBE7D','#d26d08',
#               # "#AB2A4F"
#               "#BC5572","#CD7F97","#DEABBA",'#A0CBE8','#499894','#A9E0AB','#6FC6CE')
hspcCmap2 <- c('#CE5B5B','#FF9D9A','#EDA6A6','#79706E',"#AA9BC5",'#D68645','#F28E2B','#FFBE7D',
              "#AB2A4F","#BC5572","#CD7F97","#DEABBA",'#A0CBE8','#499894','#A9E0AB','#6FC6CE')

plotUMAP(umapDF, dataType = "qualitative",namedColors=TRUE, cmap=hspcCmap,point_size=pointSize, covarLabel="RNA_labels", useRaster=T)
ggsave("./HSPC/figures/HSPC_UMAP_r1_Raster.pdf",width=5,height = 5)
dev.off()

# Figure 1D
# BMMC RNA labels UMAP
umapDF <- read.csv("./BMMC/tables/BMMC_umapDF.csv")
hspc_label <- read.csv("./BMMC/tables/BMMC_celltype_labels.csv")
umapDF <- umapDF %>% 
  left_join(.,hspc_label,by="barcode_id") %>%  
  dplyr::select(c("barcode_id","X0","X1","major_cell_type_r1"))

colnames(umapDF) <- c("barcode_id","UMAP1","UMAP2","NamedClust")
rownames(umapDF) <- umapDF$barcode_id
umapDF <- umapDF[,-1]
bmmcLmap <-c("#D7523E","#dc8e97","#eebd85","#e98741","#d9b71a","#31599B" ,"#aedd2f" ,"#499894",
"#8fc0dc","#967568" ,"#7587b1" ,"#abd0a7" ,"#794976", "#bcacd3","#e5b5b5"   , "#82c785","#edeaa4","#e1a4c6","#dbad5f"  ) 
umapDF$NamedClust <- factor(umapDF$NamedClust,levels = c(
  'HSC.MPP', 'GMP',
    'Mono','Neutrophil',
    'pDC','DC',
    'EBM','CD4_T',
 'CD8_T', 'NK',
    'gdT','Pro.PreB',
    'B','ASC',
    'MEP', 'MKP.MK',
    'Erythroid','Stromal'))
pointSize <- 0.005   
plotUMAP(umapDF, dataType = "qualitative",namedColors=TRUE, point_size=pointSize,cmap=bmmcLmap, covarLabel="RNA_labels", useRaster=T)
ggsave("./BMMC/figures/BMMC_UMAP_r1_RasterT.pdf",width=5,height = 5)

# Figure 1E
suppressPackageStartupMessages({
    library(ArchR)
    library(dplyr)
    library(tidyr)
    library(rtracklayer)
  })

wd <- "./ATAC/filtered.cells.integrated.final"
setwd(wd)

data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

pointSize <- 0.25
barwidth <- 0.9

atac_proj <- loadArchRProject(wd,force = T)
plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")

# ATAC labels UMAP
pointSize <- 0.01
umapDF <- buildUMAPdfFromArchR(atac_proj, cellColData="NamedClust",embeddingName="UMAP")
umapDF$NamedClust <- factor(umapDF$NamedClust,
                            levels=factor(c('HSC.MPP','GMP','Mono','Neutrophil','pDC','EBM','T.NK','ProB','PreB','B','ASC','MEP.MKP','Erythroid'))) #'CLP'
atacClustCmap <- c("#CE5B5B","#D68645" ,"#eebd85","#e98741","#d9b71a","#aedd2f","#499894","#abd0a7", "#CD7F97","#794976", "#bcacd3","#82c785","#edeaa4")
plotUMAP(umapDF, dataType = "qualitative",namedColors=TRUE, point_size=pointSize,cmap=atacClustCmap, covarLabel="ATAC_labels", useRaster=T)
ggsave("../figures/UMAP_r1_Raster.pdf",width=5,height = 5)


# ---------------------------------------
#       Figure 1F-G
# ---------------------------------------
source("./cellular_analysis/roe_util.R")
options(stringsAsFactors = F)

#ROE - BMMC heatmap ------------
raw <- read.csv("./BMMC/tables/major_anno.update.csv") 
data <- raw %>%
  filter(!Sample %in% c("SLE-ITP_5P","SLE-ITP_6P" ,"SLE-ITP_8P","SLE-ITP_9P","SLE-ITP_10P" ,"SLE-ITP_11P")) %>% 
  mutate(type="HSPC")
Roe.res <- Roe(data, condition="type", cellType="cell_type_r1", samples="Sample", ctrl=c("HC_5","HC_6","HC_7"))

Cluster_levels <- c('HSC.MPP','GMP',
    'Mono','Neutrophil','pDC','DC','EBM',
    'CD4_T','CD8_T','NK', 'gdT',
    'Pro.PreB','B','ASC',
    'MEP','MKP.MK','Erythroid',
    'Stromal')
df <- Roe.res$samples %>% 
  dplyr::select("cellType","samples","O2E") %>% 
  pivot_wider(names_from = cellType, values_from = O2E) %>% 
  column_to_rownames("samples")

plot_ROE_heatmap(t(df),rev(Cluster_levels)) + coord_flip()

# ROE - HSPC heatmap ------------
raw <- read.csv("./HSPC/tables/HSPC_1st_Anno_meta.csv")
data <- raw %>% 
  filter(!Sample %in% c("SLE-ITP_6P","SLE-ITP_8P","SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P")) %>%
  mutate(type="HSPC")
Roe.res <- Roe(data, condition="type", cellType="cell_type", samples="Disease", ctrl=c("HC"))

df <- Roe.res$samples %>% 
  dplyr::select("cellType","samples","O2E") %>% 
  pivot_wider(names_from = cellType, values_from = O2E) %>% 
  column_to_rownames("samples")
plot_ROE_heatmap(t(df),rev(Cluster_levels)) + coord_flip()


# ---------------------------------------
#       Figure 1H-I
# ---------------------------------------

library(distdimscr)
library(tidyverse)
library(cowplot)
options(stringsAsFactors = F)


tmp <- read.csv("./HSPC/tables/HSPC_pca_harmonyDF.csv",check.names = F) %>% 
  column_to_rownames("barcode_id")
hspc_label <- read.csv("./HSPC/tables/HSPC_bc_meta.csv") %>% 
  filter(!Sample %in% c("SLE-ITP_6P" ,"SLE-ITP_8P","SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P")) #%>% 

### SLE-TP vs. HC
cellType <- unique(hspc_label$cell_type)
unique(hspc_label$Sample)
{
  bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
  names(bhatt.dist)<-cellType
  names(bhatt.dist.rand)<-cellType
  
  group_oi = "SLE-ITP"
  group_ctl = "HC"
}

for (CT in cellType){
  print(CT)
  
  for (j in 1:10){
    
    set.seed(j)
    
    n1=length(which(hspc_label$Disease==group_ctl & hspc_label$cell_type==CT))
    n2=length(which(hspc_label$Disease==group_oi & hspc_label$cell_type==CT))
    
    if (n1 >= n2){
      cells.S<-sample(hspc_label$barcode_id[which(hspc_label$Disease==group_oi & hspc_label$cell_type==CT)],n2)
      cells.H <- sample(hspc_label$barcode_id[which(hspc_label$Disease==group_ctl & hspc_label$cell_type==CT)],n2)
    }else{
      cells.S<-sample(hspc_label$barcode_id[which(hspc_label$Disease==group_oi & hspc_label$cell_type==CT)],n1)
      cells.H <- sample(hspc_label$barcode_id[which(hspc_label$Disease==group_ctl & hspc_label$cell_type==CT)],n1)
    }
    
    
    # tmp<-Sobj@reductions$harmony@cell.embeddings
    
    cells.H.pca <- tmp[cells.H,]
    cells.S.pca <- tmp[cells.S,]
    
    for (i in 1:10) {
      
      d<-(j-1)*10+i
      
      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.H.pca,embed_mat_y=cells.S.pca,dims_use=1:30,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=FALSE)
      
      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.H.pca,embed_mat_y=cells.S.pca,dims_use=1:30,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=TRUE)
      
    }
  }
}

{
  res <- as.data.frame(bhatt.dist-bhatt.dist.rand) 
  meas_df <- as.data.frame(apply(res, 2, median)) %>%
    mutate(cell_type = rownames(.))
  colnames(meas_df)[1] <- "meansV"
  sort <- res[,order(apply(res, 2, median))]
  res <- as.data.frame(bhatt.dist-bhatt.dist.rand) %>% 
    pivot_longer(everything(),names_to = c("cell_type"),values_to = ("value")) 
  res <- res %>% 
    left_join(.,meas_df,by="cell_type")
  res$cell_type <- factor(res$cell_type,levels = colnames(sort))
  ggplot(res, aes(x=fct_rev(cell_type), y=value,fill=meansV)) + 
    geom_violin(trim = T) +
    geom_boxplot(width=0.1,outlier.shape = NA) +
    labs(y="",x="",title="") + 
    theme_cowplot() +
    theme(
      # axis.text.x = element_text(angle=60, vjust=0.6,size=12,hjust = 1,family = "Arial", color="black")
      axis.text.x = element_text(size = 10,angle=45,hjust=1,vjust=1,color="black"), #,family = "Arial", 
      axis.text.y = element_text(size = 10 ,color="black"),
    )  +
    scale_fill_gradient2(low = "#889b5d",
                        mid = "#4e9592",
                        high = "#dbad5f",n.breaks=3,labels=c("","",""),name="") 
}


ggsave("./HSPC/figures/Bhatt_distance-SLE-ITP.vs.inhouseHC_100cells.pdf",height = 2.5,width = 6)


# ----------------------------------------
tmp <- read.csv("./BMMC/tables/BMMC_pca_harmonyDF.csv",check.names = F) %>% 
  column_to_rownames("barcode_id")
hspc_label <- read.csv("./BMMC/tables/BMMC_bc_meta.csv") %>% 
  filter(!Sample %in% c("SLE-ITP_5P" ,"SLE-ITP_6P" ,"SLE-ITP_8P","SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P")) 
### SLE-TP vs. HC
cellType<- unique(hspc_label$major_cell_type_r1)
table(hspc_label$major_cell_type_r1,hspc_label$Disease)
unique(hspc_label$Sample)
{
  bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
  names(bhatt.dist)<-cellType
  names(bhatt.dist.rand)<-cellType
  
  group_oi = "SLE-ITP"
  group_ctl = "HC"
}

for (CT in cellType){
  print(CT)
  
  for (j in 1:10){
    
    set.seed(j)
    
    n1=length(which(hspc_label$Disease==group_ctl & hspc_label$major_cell_type_r1==CT))
    n2=length(which(hspc_label$Disease==group_oi & hspc_label$major_cell_type_r1==CT))
    
    if (n1 >= n2 & n2>50){
      cells.S<-sample(hspc_label$barcode_id[which(hspc_label$Disease==group_oi & hspc_label$major_cell_type_r1==CT)],n2)
      cells.H <- sample(hspc_label$barcode_id[which(hspc_label$Disease==group_ctl & hspc_label$major_cell_type_r1==CT)],n2)
    }else if(n2 >=n1 & n1 > 50){
      cells.S<-sample(hspc_label$barcode_id[which(hspc_label$Disease==group_oi & hspc_label$major_cell_type_r1==CT)],n1)
      cells.H <- sample(hspc_label$barcode_id[which(hspc_label$Disease==group_ctl & hspc_label$major_cell_type_r1==CT)],n1)
    }else (next)
    
    
    # tmp<-Sobj@reductions$harmony@cell.embeddings
    
    cells.H.pca <- tmp[cells.H,]
    cells.S.pca <- tmp[cells.S,]
    
    for (i in 1:10) {
      
      d<-(j-1)*10+i
      
      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.H.pca,embed_mat_y=cells.S.pca,dims_use=1:30,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=FALSE)
      
      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.H.pca,embed_mat_y=cells.S.pca,dims_use=1:30,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=TRUE)
      
    }
  }
}

{
  res <- as.data.frame(bhatt.dist-bhatt.dist.rand)
  # sort <- res[,order(colMeans(res))]
  # meas_df <- as.data.frame(colMeans(sort)) %>% 
  meas_df <- as.data.frame(apply(res, 2, median)) %>%
    mutate(cell_type = rownames(.))
  colnames(meas_df)[1] <- "meansV"
  sort <- res[,order(apply(res, 2, median))]
  res <- as.data.frame(bhatt.dist-bhatt.dist.rand) %>% 
    pivot_longer(everything(),names_to = c("cell_type"),values_to = ("value")) 
  res <- res %>% 
    left_join(.,meas_df,by="cell_type")
  res$cell_type <- factor(res$cell_type,levels = colnames(sort))
  res <- na.omit(res)
  ggplot(res, aes(x=fct_rev(cell_type), y=value,fill=meansV)) +
    geom_violin(trim = T) +
    geom_boxplot(width=0.1,outlier.shape = NA) +
    labs(y="",x="",title="") +
    theme_cowplot() +
    theme(
      # axis.text.x = element_text(angle=60, vjust=0.6,size=12,hjust = 1,family = "Arial", color="black")
      axis.text.x = element_text(size = 10,angle=45,hjust=1,vjust=1,color="black"), #,family = "Arial",
      axis.text.y = element_text(size = 10 ,color="black"),
    )  +
    scale_fill_gradient2(low = "#889b5d",
                         mid = "#4e9592",
                         high = "#dbad5f",n.breaks=3,name="")
}
ggsave("./BMMC/figures/Bhatt_distance-SLE-ITP.vs.ITP_100cells.pdf",height = 2.5,width = 6)



table(hspc_label$Disease,hspc_label$major_cell_type_r1)

# atac

tmp <- read.csv("./atac_Harmony_LSIdf.csv",check.names = F,row.names=1) #%>%
head(tmp)
hspc_label <- read.csv("./atac_bc_meta.csv") %>% 
  column_to_rownames("X") %>% 
  mutate(Disease=str_extract(Sample, "^[^_]+")) %>%
  mutate(barcode_id=rownames(.)) %>%
  filter(!Sample %in% c("SLE-ITP_5-Post","SLE-ITP_6-Post","SLE-ITP_8-Post","SLE-ITP_9-Post","SLE-ITP_10-Post","SLE-ITP_11-Post"))

{
  res <- as.data.frame(bhatt.dist-bhatt.dist.rand)
  meas_df <- as.data.frame(apply(res, 2, median)) %>%
    mutate(cell_type = rownames(.))
  colnames(meas_df)[1] <- "meansV"
  sort <- res[,order(apply(res, 2, median))]
  res <- as.data.frame(bhatt.dist-bhatt.dist.rand) %>% 
    pivot_longer(everything(),names_to = c("cell_type"),values_to = ("value")) 
  res <- res %>% 
    left_join(.,meas_df,by="cell_type")
  res$cell_type <- factor(res$cell_type,levels = colnames(sort))
  res <- na.omit(res)
  ggplot(res, aes(x=fct_rev(cell_type), y=value,fill=meansV)) +
    geom_violin(trim = T) +
    geom_boxplot(width=0.1,outlier.shape = NA) +
    labs(y="",x="",title="") +
    theme_cowplot() +
    theme(
      # axis.text.x = element_text(angle=60, vjust=0.6,size=12,hjust = 1,family = "Arial", color="black")
      axis.text.x = element_text(size = 10,angle=45,hjust=1,vjust=1,color="black"), #,family = "Arial",
      axis.text.y = element_text(size = 10 ,color="black"),
    )  +
    scale_fill_gradient2(low = "#889b5d",
                         mid = "#4e9592",
                         high = "#dbad5f",n.breaks=3,name="")
}

ggsave("./figures/Bhatt_distance-majorclust-SLE-ITP.vs.HC_100cells.pdf",height = 2.5,width = 6)

