library(ggpubr)
source("./cellular_analysis/roe_util.R")
options(stringsAsFactors = F)

# ---------------------------------------
#       Figure 5A
# ---------------------------------------

raw <- read.csv("./HSPC/tables/major_anno_meta_filter.csv")
data <- raw %>% 
  mutate(type="HSPC")
Roe.res <- Roe(data, condition="type", cellType="cell_type", samples="Sample", ctrl=c("HC_1","HC_2","HC_3","HC_4","HC_5","HC_6","HC_7"))

df <- Roe.res$samples %>% 
  select("cellType","samples","O2E") %>% 
  pivot_wider(names_from = cellType, values_from = O2E) %>% 
  column_to_rownames("samples")
roe_ <- df %>% 
  mutate(Sample = rownames(.)) %>% 
  select(c("MKP",'Sample')) %>% 
  filter(!Sample %in% c("SLE-NP_1","SLE-NP_2")) %>% 
  mutate(group = ifelse(Sample %in% c("ITP_1","ITP_2","ITP_3","ITP_4"),"ITP",
                        ifelse(Sample %in% c("SLE-ITP_1","SLE-ITP_3","SLE-ITP_8"),"SLE-ITP.MKP low","SLE-ITP.MKP high")
                        ))
roe_$group <- factor(roe_$group,levels = c("ITP","SLE-ITP.MKP high","SLE-ITP.MKP low"))
ggdotchart(roe_, x = "Sample", y = "MKP",
           color = "group",                                # Color by groups
           palette = alpha(c("#9E6D7F","#8A2A46","#e5b5b5"),0.8),
           # palette= c("#9E6D7FB2", "#8A2A46B2", "#E5B5B5B2"),
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 8,                                 # Large dot size
           label = round(roe_$MKP,2),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
)+
  geom_hline(yintercept = 1, linetype = 2, color = "black")+theme(legend.position='none')+ylab("Ro/e")+ggtitle("")+
  coord_flip()

# ---------------------------------------
#       Figure 5C
# ---------------------------------------
library(ArchR)
library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SCAVENGE)
library(parallel)
library(data.table)

{
  numThreads = detectCores()/4
  viridis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF",
              "#5DC863FF", "#AADC32FF", "#FDE725FF")


  subgroup <- "HSPC"
  wd <- sprintf("./ATAC/subclustered_%s", subgroup)
  atac_proj <- loadArchRProject(wd, force=TRUE)

  plotDir <- paste0(atac_proj@projectMetadata$outputDirectory, "/Plots")
}

{
  proj_PeakMatrix <- getMatrixFromProject(
  ArchRProj = atac_proj,
  useMatrix = "PeakMatrix",
  )

  peakbycellmat <- assay(proj_PeakMatrix)
  SE_gvar <- SummarizedExperiment(assays = list(counts = peakbycellmat),
                            rowRanges = rowRanges(proj_PeakMatrix), 
                            colData = DataFrame(names = colnames(peakbycellmat)))

  assayNames(SE_gvar) <- "counts"
  save(SE_gvar, file="./HSPC_SE_gvar.rda") 

  SE_gvar <- addGCBias(SE_gvar, genome = BSgenome.Hsapiens.UCSC.hg38) # use the right genome
  SE_gvar_bg <- getBackgroundPeaks(SE_gvar, niterations=200)
  save(SE_gvar, SE_gvar_bg, file="HSPC_SE_gvar.gc.bg.rda") 
  mutualknn30 <- getmutualknn(lsimat = atac_proj@reducedDims$Harmony@listData[['matDR']], num_k = 30)
  save(mutualknn30, file="HSPC_mutualknn30.rda") 
}


SE_data <- SE_gvar
SE_data_bg <- SE_gvar_bg

fileDir="./LDSC/data/GWAS-stats/fmSNP_PICS/forSCAVENGE/"
flst = list.files(fileDir)
flst = list.files(fileDir,pattern = "*.bed",include.dirs = FALSE)

for (k in 1:length(flst)) {

  
  # Pull trait info:
  f = flst[k]
  trait_file = paste0(fileDir,"/",f)
  traitName = gsub("\\..*","",f)
  # traitName = gsub("","",traitName)
  print(paste0(k,"/",length(flst),": ",traitName))
  
  # Read in SNP data and perform gchromvar Z-scores:
  trait_import <- importBedScore(rowRanges(SE_data), trait_file, colidx=5)
  SE_data_DEV <- computeWeightedDeviations(SE_data, trait_import, background_peaks =
                                             SE_data_bg)
  z_score_mat <- data.frame(colData(SE_data), z_score=t(assays(SE_data_DEV)[["z"]]) %>% c)
  
  # Seed cells:
  seed_idx <- seedindex(z_score_mat$z_score, 0.05)
  scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)
  
  # Network propagation from the seed cells
  np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)
  
  # remove cells not connected to other cells (singleton cells)
  omit_idx <- np_score==0
  mutualknn30.2 <- mutualknn30[!omit_idx, !omit_idx]
  np_score <- np_score[!omit_idx]
  
  # Calculate TRS:
  TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
  TRS <- TRS * scale_factor
  
  # Output matrix:
  trait_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)
  trait_mat = merge(data.frame(names=rownames(atac_proj@cellColData),
                               # cell=atac_proj@cellColData$cell,
                               Sample=atac_proj@cellColData$Sample,
                               group=atac_proj@cellColData$group,
                               Disease=atac_proj@cellColData$Disease,
                               group_celltype=atac_proj@cellColData$group_celltype,
                               disease_celltype=atac_proj@cellColData$disease_celltype,
                               FineClust=atac_proj@cellColData$FineClust,
                               FineClust_RNA=atac_proj@cellColData$FineClust_RNA,
                               # cell.types.final=atac_proj@cellColData$cell.types.final,
                               Clusters=atac_proj@cellColData$Clusters_Harmony,
                               UMAP_1=atac_proj@embeddings@listData$UMAP_Harmony[['df']][,1],
                               UMAP_2=atac_proj@embeddings@listData$UMAP_Harmony[['df']][,2]),
                    trait_mat,
                    by='names')
  
  nperm = 1000
  trait_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30.2, seed_idx=trait_mat$seed_idx,
                                    topseed_npscore=trait_mat$np_score, permutation_times=nperm,
                                    true_cell_significance=0.05, rda_output=F, mycores=numThreads, rw_gamma=0.05)
  colnames(trait_permu)[2] <- "enriched"
  trait_mat <- data.frame(trait_mat, trait_permu)
  # trait_mat[order(trait_mat$TRS,decreasing = T)[1:30],]
  
  # Save output matrix:
  f.out = paste0("./post-GWAS/scavenge/",traitName,".txt")
  trait_mat.save = trait_mat[,c("names","Sample","group","Disease","group_celltype"
                                ,"disease_celltype","FineClust","FineClust_RNA","Clusters","UMAP_1","UMAP_2","TRS","enriched")]
  fwrite(trait_mat.save,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
  # TRS U-MAP:
  p2 <- ggplot(data=trait_mat.save, aes(UMAP_1, UMAP_2, color=TRS)) + 
    geom_point(size=1, na.rm = TRUE, alpha = 0.6) +
    scale_color_gradientn(colors = viridis) + 
    scale_alpha()+
    theme_bw() + 
    theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
  f.plot = paste0("./post-GWAS/scavenge_plots/",traitName,".pdf")
  pdf(f.plot,width=7,height=7)
  print(p2)
  dev.off()
  
  # TRS boxplot:
  g=ggplot(trait_mat.save,aes(x=FineClust,y=TRS,fill=group)) +
    theme_bw() +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1") +
    theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
    labs(x="Cell type",y=paste0("TRS (",traitName,")"))
  f.plot = paste0("./post-GWAS/scavenge_boxplots/boxplot.",traitName,".pdf")
  pdf(f.plot,width = 10,height=6)
  print(g)
  dev.off()
  
}
for (traitName in c("Platelet_count")) { 
  f = paste0(traitName,".txt")
  df = fread(f,data.table = F,stringsAsFactors = F)
}

df %>% filter(FineClust %in% c("HSC.MPP")) %>% 
  group_by(Sample,group) %>% 
  summarise(TRS_ex = mean(enriched) * 100) %>% 
  group_by(group) %>% 
  summarise(mean = mean(TRS_ex) ) 
