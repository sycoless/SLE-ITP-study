setwd("./RNA/")
# load library and function
source("./cellular_analysis/roe_util.R")
options(stringsAsFactors = F)
sample_order <- paste(rep("SLE-ITP_",14),1:14,sep="")

## -------------------------------------
##  MonoC1Q score
## -------------------------------------
df.myeloid.score <- read.csv("./RNA/BMMC/subclustering/Myeloid/myeloid.signature.csv") %>%
  filter(cell_type_r1 == "Mono_C1Q") %>%
  mutate(M1=M1.blood_gene_signature_ucell_score) %>%
  mutate(M2=M2.blood_gene_signature_ucell_score) %>%
  mutate("M1/M2"=M2.blood_gene_signature_ucell_score) %>%
  select(Sample,M1,M2,`M1/M2`)
df.myeloid.score.summ <- df.myeloid.score %>% group_by(Sample) %>%
  summarise(across(starts_with("M"), mean, na.rm = TRUE)) %>%
  filter(!Sample %in% c(
    # "SLE-ITP_1",
    "SLE-ITP_5P"
    , "SLE-ITP_6P","SLE-ITP_8P",
  "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P","SLE-ITP_8-Post","SLE-ITP_6-Post"
  )) %>%
  filter(!Sample %in% c("HC_1","HC_2","HC_3","HC_4","HC_5","HC_6","HC_7")) %>%
  filter(!Sample %in% c("ITP_1","ITP_2","ITP_3","ITP_4")) %>%
  filter(!Sample %in% c("SLE-NP_1","SLE-NP_2","SLE-NP_3")) %>%
  column_to_rownames("Sample") %>%
  t() %>%
  as.data.frame()
rownames(df.myeloid.score.summ) <- paste(rownames(df.myeloid.score.summ)," score",sep = "")


## -------------------------------------
##  bcr
## -------------------------------------
bcr <- read.csv("./RNA/BMMC/tables/Bcells.BCRclone.csv") %>%
  filter(cell_type_r1=="ASC") %>%
  select(Sample,Clone_state,Proportion) %>%
  group_by(Sample,Clone_state) %>%
  summarise(across(starts_with("Proportion"), mean, na.rm = TRUE)) %>%
  pivot_wider(names_from=Sample,values_from=Proportion,id_cols = Clone_state) %>%
  as.data.frame() %>%
  column_to_rownames("Clone_state")
rownames(bcr) <- paste("ASC",rownames(bcr),sep="-")
bcr[is.na(bcr)] <- 0 


## -------------------------------------
## ASC score
## -------------------------------------
library(Seurat)
asc.nmf <- read.csv("./RNA/BMMC/tables/ASC_NMF_score.csv") %>%
  select(Sample,MP1:MP5)
head(asc.nmf)
df.nmf <- asc.nmf %>% group_by(Sample) %>%
  summarise(across(starts_with("MP"), mean, na.rm = TRUE)) %>%
  filter(!Sample %in% c(
    # "SLE-ITP_1",
    "SLE-ITP_5P"
    , "SLE-ITP_6P","SLE-ITP_8P",
  "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P","SLE-ITP_8-Post","SLE-ITP_6-Post"
  )) %>%
  filter(!Sample %in% c("HC_1","HC_2","HC_3","HC_4","HC_5","HC_6","HC_7")) %>%
  filter(!Sample %in% c("ITP_1","ITP_2","ITP_3","ITP_4")) %>%
  filter(!Sample %in% c("SLE-NP_1","SLE-NP_2","SLE-NP_3")) %>%
  column_to_rownames("Sample") %>%
  t()
rownames(df.nmf) <- paste("ASC: ", rownames(df.nmf)," score",sep="")

## -------------------------------------
##  MKP maturation score
## -------------------------------------
mkp <- read.csv("./RNA/HSPC/tables/MKP.MK_geneset_Ucell_score.csv")
mkp <- mkp %>% filter(cell_type %in% c("MKP"))
mkp = mkp[,c('Sample','Polyploidization','MK.maturation','Thrombopoiesis')] #
mkp = mkp %>% group_by(Sample) %>%
    summarize(across(where(is.numeric),mean)) %>%
    column_to_rownames("Sample") %>%
    t() %>%
    as.data.frame()

## -------------------------------------
##  PRS score
## -------------------------------------
prs <- read.csv("./sleitp_prs3.csv")
# prs$PRS_score4 <- scale(prs$PRS_score4, center = T, scale = T)
prs = prs[,c('Sample','PRS_score4')]
colnames(prs) <- c("Sample","PLT Polygenic Score")
prs = prs %>%
    column_to_rownames("Sample") %>%
    t() %>%
    as.data.frame()


## -------------------------------------
##  proportion
## -------------------------------------
bmmc.roe <- read.csv("./BMMC/tables/BMMC_roe.table.csv") %>%
  column_to_rownames("X") %>%
  t() %>%
  as.data.frame()
rownames(bmmc.roe) <- paste(rownames(bmmc.roe), "Proportion", sep=" ")
hspc.roe <- read.csv("./HSPC/tables/HSPC_roe.table.csv") %>%
  column_to_rownames("X") %>%
  t() %>%
  as.data.frame()
rownames(hspc.roe) <- paste(rownames(hspc.roe), "Proportion", sep=" ")


## -------------------------------------
##  Figuer 7A
## -------------------------------------
final.df <- rbind(bmmc.roe[c("ASC Proportion"
          ),sample_order],
          df.nmf[,sample_order],
          df.myeloid.score.summ["M1/M2 score",sample_order],
          prs[,sample_order],
          mkp[c("MK.maturation"),sample_order],
          bcr[c("ASC-Clonal = 2","ASC-Clonal >= 3"),sample_order])

{
  df2 <- final.df
  mycol = c("#51574a","#e9d78e",ggsci::pal_d3("category20")(20)[1:20])
  df2 <- scale(t(df2))
  df2[df2 > 1] = 1
  df2[df2 < -0.5] = -0.5
  df2[is.na(df2)] = 0 
  clinial <- read.csv("./data/clinical/clinical_v1.csv")
  annotation <- data.frame(Sample=clinial$Sample_id,Sex=clinial$Sex,Age=clinial$Age,Platelet=clinial$Platelet.count,
                           # CAR.T=clinial$CAR.T.cohort,
                           Response=clinial$CAR.T_treatment..response) %>% 
    filter(Sample %in% rownames(df2)) %>% column_to_rownames("Sample")
  annotation$Age = ifelse(annotation$Age < 25, "<25", ifelse(annotation$Age > 50, ">50", "25~50"))
  annotation$Response <- ifelse(annotation$Response %in% c("CR","PR"),annotation$Response,"No CAR-T")
  annotation_col = list(Platelet=colorRampPalette(c("white", "#A92E29"))(50),
                        Sex = c(M = mycol[19],
                                F = mycol[12]),
                        Age = c("<25" = mycol[2], 
                                "25~50" = mycol[18], 
                                ">50" = mycol[8]),
                        Response = c(
                          "CR"="#417D8B",
                          "PR"="#F2D2D5",
                          "No CAR-T"="white"
                        )
                        
  )
  
  p=pheatmap::pheatmap(
    as.matrix(t(df2)),
    annotation_col = annotation,
    annotation_colors = annotation_col,
    clustering_distance_cols="euclidean",
    clustering_method ="ward.D2",
    # clustering_method ="complete",
    cutree_cols = 2,
    cutree_rows = 3,
    scale = "row",
    color = col,show_colnames = T)

}

## -------------------------------------
##  Figuer 7B, C
## -------------------------------------
{
    library(Seurat)
    library(ggplot2)
    library(UCell)
    library(patchwork)
    library(tidyr)
    library(dplyr)
    library(RColorBrewer)
    #install.packages("GeneNMF") #from CRAN 
    #remotes::install_github("carmonalab/GeneNMF")  # or from GitHub
    library(GeneNMF)
    setwd("./RNA/BMMC/subclustering/Bcells")
    seu <- readRDS("B_cellanno_r1.count.rds")
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    Idents(object = seu) <- seu$cell_type_r1
    seu_obj <- subset(x = seu, idents = "ASC")
    blockgene1 = grep(pattern = "^IG[HLK]", x = rownames(seu_obj),value = TRUE)
    DefaultAssay(seu_obj) <- "RNA"
    seu.list <- SplitObject(seu_obj, split.by = "Sample")
}

# geneNMF.programs <- multiNMF(seu.list, assay="RNA", k=2:6,min.exp=0.05,nfeatures = 2000,hvg.blocklist="") 
# GeneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
#                                             metric = "cosine",
#                                             specificity.weight=5,
#                                             weight.explained = 0.4,
#                                             min.confidence = 0.6,
#                                             # max.genes = 200,
#                                             nMP=6)

geneNMF.metaprograms <- readRDS("./geneNMF.metaprogramsk5mp6-final.rds")

{
    mp.genes <- geneNMF.metaprograms$metaprograms.genes
    seu_obj<- AddModuleScore_UCell(seu_obj, features=mp.genes,ncores=4, name="")
}


{
    seu_obj2 <- seu_obj
    Idents(seu_obj2) <- seu_obj2$Sample
    seu_obj2 <- subset(x = seu_obj2, idents = c('SLE-ITP_3','SLE-ITP_5','SLE-ITP_6','SLE-ITP_8','SLE-ITP_9','SLE-ITP_10','SLE-ITP_11'))
    seu_obj2$cart_Response <- ifelse(seu_obj2$Sample %in% c('SLE-ITP_3','SLE-ITP_6','SLE-ITP_8','SLE-ITP_11'),"PR","CR")
    df <- seu_obj2@meta.data %>%
    select(c("Sample","cart_Response",names(mp.genes))) %>% 
    filter(!Sample %in% c("SLE-ITP_5P","SLE-ITP_6P","SLE-ITP_8P","SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P"))
}

{
    df_avg <- df %>%
    group_by(Sample, cart_Response) %>%
    summarise(across(starts_with("MP"), mean, na.rm = TRUE)) %>%
    ungroup()
    df_long <- pivot_longer(df_avg, cols = MP1:MP5, names_to = "MP", values_to = "Value") %>%
        filter(MP=="MP3")
    library(ggpubr)

    ggplot(df_long, aes(x = MP, y = Value, fill = cart_Response)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
    geom_point(position = position_dodge(0.8), size = 2, shape = 21) +
    stat_compare_means(
        aes(group = cart_Response),  
        method = "wilcox.test", label = "p.format",method.args = list(alternative = "less")) +
    scale_y_continuous(limits=c(0.15,0.29),expand=c(0,0.02))+
    theme_classic()+
    labs(title = "MP3", y = "UCell score", x = "Meta Program (MP)",fill="") +
    scale_fill_manual(values = c("CR"="#417D8B","PR"="#F2D2D5")) +
    theme(
        plot.title = element_text(hjust=0.5,size=14),
        # axis.text.x = element_text(angle = 45, hjust = 1)
        axis.text.x = element_text(size=0,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        )

}


# ------------------------------------
# plot 2 gene
# ------------------------------------
{
    seu_obj2 <- seu_obj
    Idents(seu_obj2) <- seu_obj2$Sample
    seu_obj2 <- subset(x = seu_obj2, idents = c('SLE-ITP_3','SLE-ITP_5','SLE-ITP_6','SLE-ITP_8','SLE-ITP_9','SLE-ITP_10','SLE-ITP_11'))
    seu_obj2$cart_Response <- ifelse(seu_obj2$Sample %in% c('SLE-ITP_3','SLE-ITP_6','SLE-ITP_8','SLE-ITP_11'),"PR","CR")
    df <- seu_obj2@meta.data %>%
    select(c("Sample","cart_Response",names(mp.genes))) %>% 
    filter(!Sample %in% c("SLE-ITP_5P","SLE-ITP_6P","SLE-ITP_8P","SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P"))
}


vp_case1 <- function(gene_signature, test_sign,y_max=0.5,meth="wilcox.test"){
  plot_case1 <- function(signature){
    VlnPlot(seu_obj2, features = signature,
            # pt.size = 0.1, 
           group.by = "cart_Response",
           pt.size=0,
            y.max = y_max,cols=c(
                "#417D8B","#F2D2D5"
            )
            # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, 
    # label = "p.signif",method=meth,method.args = list(alternative ="greater"))
    label = "p.format",method=meth,method.args = list(alternative ="greater"))
  } +
  theme(
    plot.title = element_text(hjust=0.5,size=14),
    # axis.text.x = element_text(angle = 45, hjust = 1)
    axis.text.x = element_text(size=0,color = "black"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=12,color = "black"),
    )
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
#   file_name <- paste0(file_name, "_r.png")
#   ggsave(file_name, width = 14, height = 8)
}

comparisons <- list(c("PR", "CR"))
pdf("./NMF/202508_vlnplot_cart_JUN2-202622.pdf",width = 4,height=4)
# VlnPlot(seu_obj2,features=names(mp.genes),group.by="Response",pt.size=0,ncol=5)
vp_case1(gene_signature = "JUN", test_sign = comparisons,y_max=6) + labs(x="")
dev.off()

pdf("./NMF/202508_vlnplot_cart_DUSP1-202622.pdf",width = 4,height=4)
# VlnPlot(seu_obj2,features=names(mp.genes),group.by="Response",pt.size=0,ncol=5)
vp_case1(gene_signature = "DUSP1", test_sign = comparisons,y_max=5.5) + labs(x="")
dev.off()

# pdf("./NMF/202508_vlnplot_cart_CST3.pdf",width = 4,height=4)
# # VlnPlot(seu_obj2,features=names(mp.genes),group.by="Response",pt.size=0,ncol=5)
# vp_case1(gene_signature = "CST3", test_sign = comparisons,y_max=2) + labs(x="")
# #  scale_fill_manual(list("CR"="#417D8B",
# #       "PR"="#F2D2D5"))
# dev.off()

## -------------------------------------
##  Figuer 7D
## -------------------------------------
library(ggpubr)
library(tidyverse)
prs_ <- prs %>% filter(Sample %in% c("SLE-ITP_5P",
     "SLE-ITP_3"
      , "SLE-ITP_6P","SLE-ITP_8P",
      "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P",
      "SLE-ITP_5",
      "SLE-ITP_6",
      "SLE-ITP_8",
      "SLE-ITP_9","SLE-ITP_10","SLE-ITP_11")) %>% mutate(Time = ifelse(str_ends(Sample, "P"), "Post", "Pre")) %>%
    mutate(id = str_remove(Sample, "P$")) %>%
    mutate(Response = ifelse(id %in% c("SLE-ITP_5","SLE-ITP_9","SLE-ITP_10"),"CR","PR"))
prs_$PRS_score4 <- scale(prs_$PRS_score4, center = T, scale = T)
prs_$Response <- factor(prs_$Response,levels=c("CR","PR"))
my_comparisons <- list(c("CR","PR"))

getwd()
pdf("./plt_prs_response_scaled_psignif_score4.pdf",width=3,height=3)
ggplot(prs_, aes(x = Response, y = PRS_score4, fill = Response)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_point(position = position_dodge(0.8), size = 2, shape = 21) +
    stat_compare_means(
      aes(group = Response),
      comparisons=my_comparisons,
      method = "wilcox.test", 
    #   label = "p.signif",
      label = "p.format",
      method.args = list(alternative = "greater")
      ) +
  theme_classic()+
  scale_y_continuous(limits=c(-2,1.8),expand=c(0,0.02))+
  labs(title = "PLT PGS", y = "Polygenic Score ", x = "",fill="") +
  scale_fill_manual(values = c("CR"="#417D8B","PR"="#F2D2D5")) +
  theme(
    plot.title = element_text(hjust=0.5,size=14),
    axis.text.x = element_text(size=0,color = "black"),
    axis.text.y = element_text(size=12,color = "black"))
dev.off()
