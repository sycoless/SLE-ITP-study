{
  # library(RRHO)
  library(viridis)
  library(ggpubr)
  library(Seurat)
  library(patchwork)
  library(cowplot)
  library(tidyverse)
  
}

dbs <-c('GO_Biological_Process_2021','MSigDB_Hallmark_2020','UK_Biobank_GWAS_v1','GWAS_Catalog_2023','Reactome_Pathways_2024', 'WikiPathways_2024_Human', 'KEGG_2021_Human')

degs_ad <- read.csv(paste0('./SLE-ITP_ct/SLE-ITP_ct_MKP.csv')) %>% 
  filter(pct.1 > 0.1)
degs_adds <- read.csv(paste0('./SLE-NP_ct/SLE-NP_ct_MKP.csv')) %>% 
  filter(pct.1 > 0.1)


colnames(degs_ad)
# plot settings
logfc_thresh <- 0.2
combined_output <- data.frame()

genes.keep <- intersect(degs_ad$gene, degs_adds$gene)
cur_degs_ad <- subset(degs_ad, gene %in% genes.keep)
cur_degs_adds <- subset(degs_adds, gene %in% genes.keep)

# make sure they are in the same order:
rownames(cur_degs_ad) <- cur_degs_ad$gene
rownames(cur_degs_adds) <- cur_degs_adds$gene
cur_degs_adds <- cur_degs_adds[cur_degs_ad$gene,]

colnames(cur_degs_ad$p_val_adj)
# get gene sets for enrichr:
up_ad <- cur_degs_ad %>% subset(avg_log2FC >= logfc_thresh & p_val_adj < 0.05) %>% .$gene
down_ad <- cur_degs_ad %>% subset(avg_log2FC <= -logfc_thresh & p_val_adj < 0.05) %>% .$gene
up_adds <- cur_degs_adds %>% subset(avg_log2FC >= logfc_thresh & p_val_adj < 0.05) %>% .$gene
down_adds <- cur_degs_adds %>% subset(avg_log2FC <= -logfc_thresh & p_val_adj < 0.05) %>% .$gene

# list of inputs to enrichr
input_list <- list(
  up_ad = up_ad[!(up_ad %in% up_adds)],
  down_ad = down_ad[!(down_ad %in% down_adds)],
  up_adds = up_adds[!(up_adds %in% up_ad)],
  down_adds = down_adds[!(down_adds %in% down_ad)],
  up_both = intersect(up_ad, up_adds),
  down_both = intersect(down_ad, down_adds),
  up_ad_down_adds = intersect(up_ad, down_adds),
  up_adds_down_ad = intersect(up_adds, down_ad)
)

# size of lists
lapply(input_list, function(x){
  print(length(x))
})

library(enrichR)
enriched_df <- do.call(rbind, lapply(names(input_list), function(x){
  if(length(input_list[[x]]) > 0){
    cur_enrich <- enrichr(input_list[[x]], dbs)
  } else{return(data.frame())}
  cur_df <- do.call(rbind, lapply(dbs, function(cur_db){
    df <- cur_enrich[[cur_db]]
    if(nrow(df) > 1){df$degs <- x; 
    # df$group <- cur_group; 
    df$db <- cur_db}
    else{df <- data.frame()}
    df
  }))
}))

combined_output2 <- data.frame()
combined_output2 <- rbind(combined_output2, enriched_df) 
res_out <- combined_output2%>% 
  filter(Adjusted.P.value < 0.05) 

wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
unique(res_out$degs)


# ---------------------------------------
#       Figure 4B
# ---------------------------------------

{
  GO_result <- res_out
  GO_result$Adjusted.P.value <- as.numeric(GO_result$Adjusted.P.value)
  GO_result <- GO_result %>% 
    mutate(log10q = -log10(Adjusted.P.value)) %>% 
    filter(Term %in% c('Interferon Alpha Beta Signaling',
                       'Interferon Gamma Signaling',
                       'regulation of cell cycle (GO:0051726)',
                       'T cell receptor signaling pathway (GO:0050852)',
                       'Interleukin-12 Signaling',
                       'regulated exocytosis (GO:0045055)',
                       'MAP2K and MAPK Activation',
                       'cellular response to type I interferon (GO:0071357)',
                       'interleukin-27-mediated signaling pathway (GO:0070106)',
                       "negative regulation of programmed cell death (GO:0043069)",
                       'Response to Elevated Platelet Cytosolic Ca2+',
                       'Platelet Degranulation',
                       "regulation of cell migration (GO:0030334)",
                       "positive regulation of tumor necrosis factor production (GO:0032760)",
                       'Platelet Count',
                       'Factors Involved in Megakaryocyte Development and Platelet Production',
                       'regulation of megakaryocyte differentiation (GO:0045652)',
                       'Platelet Homeostasis',
                       "Myc Targets V1",
                       'MAPK signaling pathway'
                 
    )) %>% 
    filter(!(Term %in% c('Interferon Alpha Beta Signaling','Interferon Gamma Signaling') & degs %in% c('up_ad'))) %>%
    filter(degs != "down_adds")%>%
    filter(degs != "up_adds") %>%
    filter(degs != "up_adds_down_ad")
  unique(GO_result$degs)
  GO_result$Term <- str_replace(GO_result$Term, " \\s*\\([^\\)]+\\)", "")
  GO_result$log10q = ifelse(GO_result$degs %in% c('up_ad','up_both'), GO_result$log10q , GO_result$log10q* -1)
  GO_result$yax <- ifelse(GO_result$log10q > 0, -0.02, 0.02)
  GO_result$log10q <- ifelse(GO_result$log10q>15,15,GO_result$log10q)
  # Plot the GO results
  GO_result$wrap <- wrapText(GO_result$Term, 45)
  ggplot(GO_result) +
    geom_bar(aes(
      x = reorder(Term, log10q),
      y = log10q,
      fill = degs
    ),
    stat = "identity",
    color = "white") +
    geom_text(aes(
      x = reorder(Term, log10q),
      y = yax,
      label = Term,
      hjust = log10q > 0
    ),
    size = 7 * 0.35,
    angle = 0) +
    scale_fill_manual(
      name = "",
      values = c("up_ad" = "#FB8072", "down_ad" = "#BEBADA","up_both" ="#8DD3C7","down_both" ="gray")
    ) +
    coord_flip() +
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_text(size = 7),
      text = element_text(size = 7),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_blank(),
    ) +
    labs(x = "", y = "") 
}

# ---------------------------------------
#       Figure 4A
# ---------------------------------------

degs_ad <- read.csv(paste0('./SLE-ITP_ct/SLE-ITP_ct_MKP.csv')) %>% 
  filter(pct.1 > 0.1)
degs_adds <- read.csv(paste0('./SLE-NP_ct/SLE-NP_ct_MKP.csv')) %>% 
  filter(pct.1 > 0.1)

colnames(degs_ad)
# plot settings
logfc_thresh <- 0.2
combined_output <- data.frame()

genes.keep <- intersect(degs_ad$gene, degs_adds$gene)
cur_degs_ad <- subset(degs_ad, gene %in% genes.keep)
cur_degs_adds <- subset(degs_adds, gene %in% genes.keep)

# make sure they are in the same order:
rownames(cur_degs_ad) <- cur_degs_ad$gene
rownames(cur_degs_adds) <- cur_degs_adds$gene
cur_degs_adds <- cur_degs_adds[cur_degs_ad$gene,]

colnames(cur_degs_ad$p_val_adj)
# get gene sets for enrichr:
up_ad <- cur_degs_ad %>% subset(avg_log2FC >= logfc_thresh & p_val_adj < 0.05) %>% .$gene
down_ad <- cur_degs_ad %>% subset(avg_log2FC <= -logfc_thresh & p_val_adj < 0.05) %>% .$gene
up_adds <- cur_degs_adds %>% subset(avg_log2FC >= logfc_thresh & p_val_adj < 0.05) %>% .$gene
down_adds <- cur_degs_adds %>% subset(avg_log2FC <= -logfc_thresh & p_val_adj < 0.05) %>% .$gene

# list of inputs to enrichr
input_list <- list(
  up_ad = up_ad[!(up_ad %in% up_adds)],
  down_ad = down_ad[!(down_ad %in% down_adds)],
  up_adds = up_adds[!(up_adds %in% up_ad)],
  down_adds = down_adds[!(down_adds %in% down_ad)],
  up_both = intersect(up_ad, up_adds),
  down_both = intersect(down_ad, down_adds),
  up_ad_down_adds = intersect(up_ad, down_adds),
  up_adds_down_ad = intersect(up_adds, down_ad)
)


input_list_orig <- input_list

input_list <- input_list_orig[c(1,2,5,6)]
class_df <- data.frame("class"=rep(names(input_list),sapply(input_list,length)),"gene"=do.call(rbind,lapply(input_list,data.frame))[,1])
df1.full.mg <- degs_ad %>% left_join(.,class_df,by="gene")
colnames(df1.full.mg)
{
  tmp1 <- df1.full.mg
  tmp1$P <- -log10(tmp1$p_val_adj)
  tmp1 <- tmp1[,c("avg_log2FC","P","class",'gene')]
  colnames(tmp1) <- c("LFC","P","class",'names')
  tmp1$class <- factor(tmp1$class,levels = unique(class_df$class))
  tmp1$geneLabel <- NA
  genes_of_interest <- c("IGF1R","GP9","FCER1G","VWF",
    "PF4","ITGA2B",
                         "NFE2","STAT1","MX1","IFI27","OAS2","IFI6","IFI44L","XIST",
                         "RUNX1",
                         "HNRNPA2B1",
                         "PPIA",
                         "PSME2",
                         "CXCR4",
                         "PTCRA",
                         "HSPA9",
                         "TUBB1","IL2RG","STAT4","IRF9","IRF7","RGS18","MYL9"
                         )
  tmp1$geneLabel[tmp1$names %in% genes_of_interest] <- tmp1$names[tmp1$names %in% genes_of_interest]
  gene_set = read.csv("./resources/gene_set/MK_geneset.csv")
  colnames(gene_set)[2] <- "names"
  tmp2 <- tmp1 %>% 
    left_join(.,gene_set,by="names")
  library(scales)

  unique(tmp1$class)
  color_df <- data.frame("class" = c("up_ad","down_ad","up_both", "down_both"),"color"=c("#FB8072" ,"#BEBADA","#8DD3C7" ,"#FFFFB3"))
  cols <- tmp1 %>%
  distinct(class, color) %>%
  tibble::deframe()  

  tmp1 <- tmp1 %>%
    left_join(.,color_df,by="class")
  tmp1$color <- ifelse(tmp1$P < 1.30103, 'gray', tmp1$color)
  tmp1$P <- ifelse(tmp1$P > 20,20,tmp1$P)
  tmp1$LFC <- ifelse(tmp1$LFC  > 5,5,tmp1$LFC)

  tmp1 <- tmp1 %>%
   mutate(
    class_lab = recode(as.character(class),
      down_ad   = "Down-regulated in SLE-ITP",
      down_both = "Down-regulated in both SLE-ITP and SLE-NP",
      up_ad     = "Up-regulated in SLE-ITP",
      up_both   = "Up-regulated in both SLE-ITP and SLE-NP",
      "NA"      = "n.s.",          # class 是字符串 "NA" 的情况
      .default  = as.character(class)
    ),
    class_lab = if_else(is.na(class), "n.s.", class_lab)  # class 是真正 NA 的情况
  )

  cols <- tmp1 %>% distinct(class_lab, color) %>% deframe()
  tmp1 %>%
  ggplot(aes(LFC, P, color = class_lab)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = cols, name = "class") +
  theme_cowplot() 

  p <- tmp1  %>%
    ggplot(aes(x=LFC, y=P)) +
    geom_hline(yintercept=-log10(0.05), linetype='dashed')+ 
    ggrastr::rasterise(geom_point(alpha=0.5,color=tmp1 %>% .$color), dpi=500)

  p <- p +
    geom_point(
      inherit.aes=FALSE,
      data=subset(tmp1, !is.na(geneLabel)),
      aes(LFC, P),
      fill=subset(tmp1, !is.na(geneLabel)) %>% .$color,
      shape=21, size=3, color='black'
    ) +
    geom_text_repel(aes(label=geneLabel), color='black', fontface='italic',  min.segment.length=0, max.overlaps=Inf)
  p <- p + xlab(bquote("log"[2]~"(Fold Change)")) +
    ylab(bquote("-log"[10]~"(FDR)")) + 
    theme_cowplot() +
    theme(
      panel.border = element_rect(color='black', fill=NA, size=1),
      panel.grid.major = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5)
      # legend.position='bottom'
    ) 

  p 
}

# ---------------------------------------
#       Figure 4C
# ---------------------------------------
library(tidyverse)
library(data.table)
# library(rpgm)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
coul <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(30))

dat <- read.csv("./RNA/HSPC/tables/MKP-geneset-UCell-final.csv")

dat = dat[c('HC_1','HC_2','HC_3','HC_4','HC_5','HC_6','HC_7','ITP_1','ITP_2','ITP_3','ITP_4','SLE-ITP_1',
        'SLE-ITP_2','SLE-ITP_3','SLE-ITP_4','SLE-ITP_5','SLE-ITP_6','SLE-ITP_7',
        'SLE-ITP_8','SLE-ITP_9','SLE-ITP_10','SLE-ITP_11','SLE-ITP_12','SLE-ITP_13','SLE-ITP_14','SLE-NP_1','SLE-NP_2','SLE-NP_3'),]

df <- t(scale(dat))
col_fun = colorRamp2(c(min(df),0,max(df)), c("#4B98C5",  "white","#BB2B34")) 

ha = HeatmapAnnotation(Disease = factor(c(rep("HC",7),rep("ITP",4),rep("SLE-ITP",14),rep("SLE-NP",3)),levels = c("HC","ITP","SLE-ITP","SLE-NP")),
                       annotation_name_side = "left",
                       col = list( Disease=c("HC"="#B2CFE5","ITP"= "#9E6D7F" ,"SLE-ITP"="#8A2A46","SLE-NP"="#D97351"
                       )))

pdf("./RNA/HSPC/figures/MKP-geneset-heatmap.pdf",width=11,height=2.4)
Heatmap(df,name="Z-scores",
        col = col_fun,
        border_gp = gpar(col = "black", lty = 1),
        rect_gp = gpar(col = "white", lwd = 0.2),
        cluster_rows = F,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        cluster_columns = F,
        column_dend_side="bottom",
        column_dend_height = unit(8, "mm"),
        row_dend_width = unit(8, "mm"),
        column_title_gp = gpar(fontface = "bold",fontsize = 14),
        ht_opt(heatmap_column_names_gp=gpar(fontsize = 12), #fontface = "bold",
               heatmap_row_names_gp = gpar(fontsize = 14)), #fontface = "bold"
        row_names_side = "left", 
        top_annotation = ha

)
dev.off()

# ---------------------------------------
#       Figure 4D,E,F
# ---------------------------------------
{
  library(JASPAR2022)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(monaLisa)
  library(ComplexHeatmap)
  library(circlize)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(tidyverse)
  library(plyranges)
  library(data.table)
  library(ggpubr)
  library(cowplot)
}

{
  # pwms <- getMatrixSet(JASPAR2022,opts = list(matrixtype = "PWM",collection="CORE",species =9606))
  pwms <- getMatrixSet(JASPAR2022,opts = list(matrixtype = "PWM",tax_group = "vertebrates"))
  # names(pwms) <- name(pwms)
  filter_motif_data <- function(monalisa_out,log10_p_thresh = 4,type = 'negLog10Padj'){
    monalisa_filter <- apply(assay(monalisa_out, type), 1, function(x) max(abs(x), 0, na.rm = TRUE)) >log10_p_thresh 
    filtered_data <- monalisa_out[monalisa_filter,]}
  
  pull_motifs_monalisa_output <- function(monalisa_output,motif_list){
    selected_TF_indexes <- data.frame(motif = rowData(monalisa_output)$motif.name) %>%
      rownames_to_column('motif_id') %>%
      mutate(index = row_number()) %>%
      filter(toupper(motif) %in% motif_list | motif_id %in% motif_list) %>%
      pull(index)
    selected_TFs <- monalisa_output[selected_TF_indexes,]
    selected_TFs
    
    plotMotifHeatmaps2 <- function(x,
                                   which.plots = c("negLog10P", "pearsonResid", 
                                                   "negLog10Padj", "log2enr"),
                                   width = 4,
                                   col.enr = c("#053061", "#2166AC", "#4393C3",
                                               "#92C5DE", "#D1E5F0", "#F7F7F7",
                                               "#FDDBC7", "#F4A582", "#D6604D",
                                               "#B2182B", "#67001F"),
                                   col.sig = c("#F0F0F0", "#D9D9D9", "#BDBDBD",
                                               "#969696", "#737373", "#525252",
                                               "#252525", "#000000"),
                                   col.gc = c("#F7FCF5", "#E5F5E0", "#C7E9C0",
                                              "#A1D99B", "#74C476", "#41AB5D",
                                              "#238B45", "#006D2C", "#00441B"),
                                   maxEnr = NULL,
                                   maxSig = NULL,
                                   highlight = NULL,
                                   cluster = FALSE,
                                   show_dendrogram = FALSE,
                                   show_motif_GC = FALSE,
                                   show_seqlogo = FALSE,
                                   show_bin_legend = FALSE,
                                   width.seqlogo = 1.5,
                                   use_raster = FALSE,
                                   na_col = "white", 
                                   doPlot = TRUE,
                                   ...) {
      stopifnot(exprs = {
        is(x, "SummarizedExperiment")
        all(which.plots %in% assayNames(x))
        "bins" %in% names(metadata(x))
        (!show_motif_GC || "motif.percentGC" %in% colnames(rowData(x)))
      })
      b <- metadata(x)$bins
      stopifnot(exprs = {
        ncol(x) == nlevels(b)
        all(which.plots %in% c("negLog10P", "negLog10Padj", 
                               "pearsonResid", "log2enr"))
        is.null(highlight) || (is.logical(highlight) && 
                                 length(highlight) == nrow(x))
      })
      bincols <- attr(getColsByBin(b), "cols")
      if (identical(cluster, TRUE)) {
        clAssayName <- "pearsonResid"
        clAssay <- assay(x, clAssayName)
        allNA <- rowSums(is.na(clAssay)) == ncol(clAssay)
        if (any(allNA)) {
          warning("removing motifs without finite values in '",
                  clAssayName, "': ",
                  paste(rownames(clAssay)[allNA], collapse = ", "))
          x <- x[!allNA, ]
          clAssay <- clAssay[!allNA, ]
        }
        clres <- hclust(dist(clAssay))
      } else if (identical(cluster, FALSE)) {
        clres <- FALSE
      } else if (is(cluster, "hclust")) {
        clres <- cluster
      } else {
        stop("'cluster' must be either TRUE, FALSE or an hclust-object.")
      }
      hmBin <- HeatmapAnnotation(df = data.frame(bin = colnames(x)), name = "bin",
                                 col = list(bin = bincols),
                                 show_annotation_name = FALSE,
                                 which = "column", width = unit(width,"inch"),
                                 annotation_height = unit(width / 16, "inch"),
                                 show_legend = show_bin_legend)
      tmp <- matrix(if (!is.null(highlight)) {
        as.character(highlight) 
      } else {
        rep(NA, nrow(x))
      },
      ncol = 1, dimnames = list(unname(rowData(x)$motif.name), NULL))
      hmSeqlogo <- NULL
      if (show_seqlogo) {
        pfms <- rowData(x)$motif.pfm
        maxwidth <- max(vapply(TFBSTools::Matrix(pfms), ncol, 0L))
        grobL <- lapply(pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")
        hmSeqlogo <- HeatmapAnnotation(
          logo = annoSeqlogo(grobL = grobL, which = "row",
                             space = unit(0.5, "mm"),
                             width = unit(width.seqlogo, "inch")),
          show_legend = FALSE, show_annotation_name = FALSE, which = "row")
      }
      hmMotifs <- Heatmap(
        matrix = tmp, name = "names",
        width = unit(if (!is.null(highlight)) .2 else 0, "inch"),
        na_col = NA, col = c("TRUE" = "green3", "FALSE" = "white"),
        cluster_rows = clres, show_row_dend = show_dendrogram,
        cluster_columns = FALSE, show_row_names = TRUE,
        row_names_side = "left", show_column_names = FALSE,
        show_heatmap_legend = FALSE, left_annotation = hmSeqlogo,
        ...
      )
      
      assayNameMap1 <- c(negLog10P = "P value",
                         negLog10Padj = "adj. P value",
                         pearsonResid = "Pearson residual",
                         log2enr = "log2 enrichment")
      assayNameMap2 <- c(negLog10P = "P value (-log10)",
                         negLog10Padj = "adj. P value (-log10)",
                         pearsonResid = "Pearson residual (o-e)/sqrt(e)",
                         log2enr = "enrichment (log2)")
      L <- list(labels = hmMotifs)
      if (show_motif_GC) {
        tmp <- as.matrix(rowData(x)[, "motif.percentGC", drop = FALSE])
        hmPercentGC <- Heatmap(
          matrix = tmp, name = "Percent G+C",
          width = unit(0.2, "inch"), na_col = NA,
          col = colorRamp2(breaks = c(0, seq(20, 80, length.out = 254), 100),
                           colors = colorRampPalette(col.gc)(256)),
          cluster_rows = FALSE, cluster_columns = FALSE,
          show_row_names = FALSE, show_column_names = FALSE,
          show_heatmap_legend = TRUE,
          heatmap_legend_param = list(color_bar = "continuous"),
          use_raster = use_raster,
          ...
        )
        L <- c(L, list("percentGC" = hmPercentGC))
      }
      ret <- c(L, lapply(which.plots, function(w) {
        dat <- assay(x, w)
        if ((w == "pearsonResid") | (w == "log2enr")) {
          rng <- c(-1, 1) * if (is.null(maxEnr)) {
            quantile(abs(dat), .995, na.rm = TRUE) 
          } else {
            maxEnr
          }
          cols <- col.enr
        } else {
          rng <- c(0, 
                   if (is.null(maxSig)) {
                     quantile(dat, .995, na.rm = TRUE) 
                   } else {
                     maxSig
                   })
          cols <- col.sig
        }
        Heatmap(
          matrix = dat,
          name = assayNameMap1[w],
          width = unit(width,"inch"),
          column_title = assayNameMap2[w],
          col = colorRamp2(breaks = seq(rng[1], rng[2], length.out = 256),
                           colors = colorRampPalette(cols)(256)),
          cluster_rows = FALSE, cluster_columns = FALSE,
          show_row_names = FALSE, show_column_names = FALSE,
          # column_names_side = "bottom", 
          # column_names_max_height = unit(1.5,"inch"),
          top_annotation = hmBin, show_heatmap_legend = TRUE,
          heatmap_legend_param = list(color_bar = "continuous"),
          use_raster = use_raster,
          na_col = na_col, 
          ...
        )
      }))
      names(ret)[seq(length(ret) - length(which.plots) + 1L, length(ret))] <- 
        which.plots
      if (doPlot) {
        show(Reduce(ComplexHeatmap::add_heatmap, ret))
      }
      invisible(ret)
    }
    
  }
}


{
    load_dat <- function(up.bed,down.bed){
    lmr.down <- rtracklayer::import(con = down.bed, format = "bed")
    lmr.down$direction <- "down"
    lmr.up <- rtracklayer::import(con = up.bed, format = "bed")
    lmr.up$direction <- "up"
    lmr <- c(lmr.down,lmr.up)
    lmr
  }

  lmr <- load_dat(up.bed = "./ATAC-output/bed/marker.peak.pb_MKP.0.05.Up2·.bed",down.bed = "./ATAC-output/bed/marker.peak.pb_MKP.0.05.Down2.bed")
  lmr$direction <- factor(lmr$direction)
  input_bins <- lmr$direction
}

{
  peak_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38,lmr)
  monalisa_res <- calcBinnedMotifEnrR(seqs = peak_seqs,bins =as.factor(input_bins) ,pwmL =pwms,verbose = TRUE)
  monalisa_res
}

{
  # MKP
  selected_list <- c('NFE2','FOS::JUNB','KLF1','ETV6','MYC','FLI1','PBX1','GATA2','GATA1','NRF2','RUNX1') 
  # HSC.MPP
  # selected_list <- c('NFE2','FOS::JUNB','KLF1','FLI1','GATA2','GATA1::TAL1','NRF2','GFI1B','MYB')  #'RUNX1',
  selected_TF_indexes <- data.frame(motif = rowData(monalisa_res)$motif.name) %>%
    mutate(index = row_number()) %>%
    filter(toupper(motif) %in% selected_list) %>%
    pull(index)
  selected_TFs <- monalisa_res[selected_TF_indexes,]

  plotMotifHeatmaps2(x = selected_TFs, which.plots = c("log2enr"), 
                    width = 1, cluster = TRUE, maxEnr = 1, maxSig = 5,show_dendrogram = TRUE,show_bin_legend=T)
}

{
    pdf("./ATAC-output/figures/monalisa_motif_pb.MKP.pdf",height = 4,width = 5)
    plotMotifHeatmaps2(x = selected_TFs, which.plots = c("log2enr"), 
                    width = 1, cluster = TRUE, maxEnr = 1, maxSig = 5,show_dendrogram = TRUE,show_bin_legend=T)
    dev.off()
}

EPD_promoters <- fread("./LDSC/data/GWAS-stats/fmSNP_PICS/Hs_EPDnew_006_hg38.bed") %>% 
  separate(V4,into = c('geneSymbol','promoter_number')) %>%
  dplyr::select(-promoter_number) %>% 
  dplyr::rename('chr' =1,'start' =2,'end' = 3) %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns =  TRUE) 

EPD_ATAC_overlap_expression_annotated <- lmr %>%
  join_overlap_left(EPD_promoters) %>% 
  mutate(promoter = case_when(is.na(geneSymbol) ~ FALSE,TRUE~TRUE)) %>% 
  dplyr::select(-V5,-V6,-V7,-V8) %>%
  data.frame() %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns =  TRUE)
merge_res.df <- read.delim("./ATAC-output/bed/marker.peak.pb_MKP.all.txt")

{
    df <- EPD_ATAC_overlap_expression_annotated %>% data.frame() 
    # runs go enrichment analysis and binds results on a column
    run_GO_enrich <- function(input_df,gene_column = 'geneSymbol'){
    require(enrichR)
    dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021',
    'MSigDB_Hallmark_2020','Reactome_Pathways_2024', 'WikiPathway_2021_Human', 'KEGG_2021_Human')
    enrich_res <- enrichr(input_df %>% data.frame()  %>% pull({gene_column}),dbs) %>%
        bind_rows(.id = 'column_label')
    enrich_res
    }

    enrich_res.up <- run_GO_enrich(EPD_ATAC_overlap_expression_annotated %>% filter(direction=="up"))
    enrich_res.down <- run_GO_enrich(EPD_ATAC_overlap_expression_annotated %>% filter(direction=="down"))

    df <- enrich_res.up %>%
    mutate(log2_OR = log2(Odds.Ratio)) %>%  # 计算 log2(OR)
    filter(Term %in% c("G2 M Transition",
                        "Myc Targets V2",
                        "Interferon Alpha Response",
                        "PI3K/AKT/mTOR  Signaling")) %>% 
    arrange(desc(log2_OR))             # 按 log2(OR) 降序排序
    df$log10q =log2(df$P.value)* -1
    df$yax <- ifelse(df$log10q > 0, -0.02, 0.02)
    # 绘制条形图
    ggplot(df, aes(x = log2_OR, y = reorder(Term, log2_OR))) +
    geom_bar(stat = "identity",fill="grey") +   # 绘制水平条形图, fill = "gray"
    labs(
        title = "",
        x = expression(log[2]~"(OR) Pathway Enrichment"),
        y = ""
    ) +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 10,color="black"),     
        axis.text.x = element_text(size = 8,color="black"), 
        axis.title.x = element_text(size = 8,color="black") 
    )
    ggsave("./ATAC-output/figures/monalisa_motif_pb.MKP.enrich_res.up.pdf",width = 3.8,height = 2)


    df <- enrich_res.down %>%
    mutate(log2_OR = log2(Odds.Ratio)) %>% 
    filter(Term %in% c("VEGFA-VEGFR2 Signaling Pathway WP3888",
                        "Factors Involved in Megakaryocyte Development and Platelet Production",
                        "PKR-mediated Signaling",
                        "TNF-alpha signaling pathway WP231")) %>% 
    arrange(desc(log2_OR))             
    df$Term2 <- c("PKR-mediated Signaling","Factors Involved in Megakaryocyte \nDevelopment and Platelet Production",
                "TNF-alpha signaling pathway","VEGFA-VEGFR2 Signaling Pathway")
    df$log10q =log2(df$P.value)* -1
    df$yax <- ifelse(df$log10q > 0, -0.02, 0.02)
    # 绘制条形图
    df <- df %>% mutate(Term_wrapped = str_wrap(Term, width = 35))
    ggplot(df, aes(x = log2_OR, y = reorder(Term2, log2_OR))) +
    geom_bar(stat = "identity",fill="grey") +   # 绘制水平条形图, fill = "gray"
    labs(
        title = "",
        x = "Log2(OR) Pathway Enrichment",
        y = ""
    ) +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 10,color="black"),      # 调整Y轴字体大小
        axis.text.x = element_text(size = 8,color="black"),
        axis.title.x = element_text(size = 8,color="black"),# 调整X轴字体大小
        plot.title = element_text(hjust = 0.5, size = 14)  # 图标题居中
    )
    ggsave("./ATAC-output/figures/monalisa_motif_pb.MKP.enrich_res.down.pdf",width = 4.1,height = 2)

}

{
    dat <- df %>% left_join(.,merge_res.df,by="Name") 
    dat$direction <- ""
    dat <- dat %>% mutate(direction = ifelse(P.Value <= 0.05 & logFC > 0 ,"Up",
        ifelse(P.Value <= 0.05 & logFC < 0,"Down","")))
    color_df <- data.frame("direction" = c("Up","Down","up_both", "down_both"),"color"=c("#FB8072" ,"#BEBADA","#8DD3C7" ,"#FFFFB3"))
    dat <- dat %>%
        left_join(.,color_df,by="direction")

    dat$color <- ifelse(dat$P.Value > 0.05, 'gray', dat$color)
    p <- dat  %>%
        ggplot(aes(x=logFC, y=-log10(P.Value))) +
        geom_hline(yintercept=-log10(0.05), linetype='dashed')+ 
        ggrastr::rasterise(geom_point(alpha=0.5,color=dat %>% .$color), dpi=500)
    p <- p + xlab(bquote("log"[10]~"(Fold Change)")) +
    ylab(bquote("-log"[10]~"(P.value)")) + #Adj. P-value
    theme_cowplot() +
    theme(
      panel.border = element_rect(color='black', fill=NA, size=1),
      panel.grid.major = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position='bottom'
    ) 
    p
}

# ---------------------------------------
#       Figure 4G
# ---------------------------------------
meta <- read.csv("./matrix/WGS_metadata_v5.csv") 
prs <- read.table("./PRS-score/PRS_platelet_more-0.1-250.profile",header=T)

prs_ <- prs %>% 
    mutate(WGS_id= sapply(strsplit(IID, "\\_"), `[`, 3)) %>%
    left_join(.,meta,by="WGS_id")

prs_ <- prs %>% 
    mutate(WGS_id= sapply(strsplit(IID, "\\_"), `[`, 3)) %>%
    inner_join(.,meta,by="WGS_id") %>%
    mutate(group = ifelse(group %in% c("SLE-TP","un-SLE-TP","nr-SLE-TP"),"SLE-ITP",group)) %>% 
    mutate(group = ifelse(group %in% c("SLE-ITP"),"SLE-ITP",
              ifelse(group %in% c("HD"),"HC",group))) %>%
    mutate(group = ifelse(!group %in% c("HC","SLE-ITP"),"Other SLE",group)) %>%
    arrange(desc(SCORESUM))

my_comparisons <- list(c("HC", "Other SLE"),c("Other SLE","SLE-ITP"),c("HC", "SLE-ITP"))
library(ggpubr)
prs_$SCORESUM_scaled <- scale(prs_$SCORESUM, center = T, scale = T)
ggplot(prs_, aes(x=group, y=SCORESUM_scaled, fill=group)) +
  geom_boxplot()+
  labs(x="", y = "PLT Polygenic Score")+
  scale_x_discrete(limits=c("HC", "Other SLE","SLE-ITP"))+
  scale_fill_manual(values=c("HC"="#B2CFE5","ITP"= "#9E6D7F" ,"SLE"="gray87","Other SLE"="gray87","SLE-ITP"="#8A2A46","SLE-NP"="#D97351"))+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "none") +
  stat_compare_means(
      aes(group = group),
      comparisons=my_comparisons,
      method = "t.test", 
    #   label = "p.signif",
      label = "p.format",
      # method.args = list(alternative = "less")) +
    method.args = list(alternative = "greater")) +
   theme(legend.title = element_blank(),
          legend.text = element_text(color = 'black',size = 10, face = 'plain'),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black',size = 12,  face = 'plain'),
          axis.title = element_text(color = 'black',size = 12, face = 'plain'),
          axis.ticks = element_line(color = 'black'))

# ---------------------------------------
#       Figure 4I
# ---------------------------------------
meta <- read.csv("./matrix/WGS_metadata_v5.csv")
dat <- read.table("./prs/CXCR4/snp.summary.txt",sep="\t",header=T) %>% select(-group)
data_ <- dat %>% 
   filter(WGS_id %in% prs_$WGS_id) %>%
    inner_join(.,meta,by="WGS_id") %>%
    mutate(group = ifelse(group %in% c("SLE-TP"),"GNR-SLE-ITP",group)) %>%  #,"nr-SLE-TP"
    mutate(group = ifelse(group %in% c("SLE-ITP"),"SLE-ITP",
              ifelse(group %in% c("HD"),"HC",ifelse(group %in% c("nr-SLE-TP","un-SLE-TP"),"GR-SLE-ITP",group)))) %>%
    mutate(group = ifelse(!group %in% c("HC","GR-SLE-ITP","GNR-SLE-ITP"),"Other SLE",group))

prs2 <- prs_ %>%
  left_join(.,data_,by="WGS_id")  
prs2 <- prs2 %>% filter(WGS_id != "S549")
prs2$SCORESUM_scaled <- scale(prs2$SCORESUM, center = T, scale = T)
prs2 <- prs2 %>% filter(group.y %in% c("GNR-SLE-ITP","GR-SLE-ITP"))
my_comparisons <- list(c("GR-SLE-ITP","GNR-SLE-ITP"))
library(ggpubr)
ggplot(prs2, aes(x=group.y, y=SCORESUM_scaled, fill=group.y)) +
  geom_boxplot()+
  labs(x="", y = "PLT Polygenic Score")+
# labs(x="", y = "SLE Polygenic Risk Score")+
  scale_x_discrete(limits=c("GR-SLE-ITP","GNR-SLE-ITP"))+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "none") +
  stat_compare_means(
      aes(group = group),
      comparisons=my_comparisons,
      method = "t.test", 
    #   label = "p.signif",
      label = "p.format",
      # method.args = list(alternative = "less")) +
    method.args = list(alternative = "greater")) +
   theme(legend.title = element_blank(),
          legend.text = element_text(color = 'black',size = 10, face = 'plain'),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black',size = 12,  face = 'plain'),
          axis.title = element_text(color = 'black',size = 12, face = 'plain'),
          axis.ticks = element_line(color = 'black'))
ggsave("./plt_prs_scaled-GR-GNR.pdf",width=3.5,height=4)
