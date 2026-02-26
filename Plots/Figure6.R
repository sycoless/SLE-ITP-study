library(tidyverse)
options(stringsAsFactors = F)
library(rstatix)
# ---------------------------------------
#       Figure 6D
# ---------------------------------------

filter_all_df <- read.csv("./BCR/changeO-Final.csv")
{
  filter_all_df_ <- filter_all_df %>% 
    filter(cell_type_r1 != "Doublet") %>% 
    mutate(id_clone_id = paste0(Sample,clone_id)) 
  clone_sizes <- countClones(dplyr::filter(filter_all_df_, locus == "IGH"),
                             groups = "Sample",clone="clone_id")
  clone_sizes <- clone_sizes %>%  mutate(id_clone_id = paste0(Sample,clone_id))
  clone <- filter_all_df_ %>% 
    select(barcode_id,Disease,cell_type_r1,id_clone_id) %>% 
    left_join(.,clone_sizes,by="id_clone_id")
  
}

{
    clone2 <- clone %>% 
        filter(Sample %in% c(
        "SLE-ITP_5P"
        , "SLE-ITP_6P","SLE-ITP_8P",
        "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P",
        "SLE-ITP_5"
        , "SLE-ITP_6","SLE-ITP_8",
        "SLE-ITP_9","SLE-ITP_10","SLE-ITP_11"
    )) %>%
    mutate(Clone_state=ifelse((seq_count==1),"No Clonal",
                                ifelse((  seq_count==2 ),"Clonal = 2",# @seq_count>=2 &
                                        "Clonal >= 3")))
    clone2[which(is.na(clone2$Clone_state)),"Clone_state"] <- "Non-VDJ"

    df <- clone2 %>% 
        group_by(Sample,cell_type_r1,Clone_state,.drop=FALSE) %>% 
        summarise(n=n())
    df <- df %>% 
    group_by(Sample, cell_type_r1) %>%
    mutate(Proportion = n / sum(n) * 100) %>%
    ungroup() %>% 
    mutate(Disease=as.character(lapply(strsplit(df$Sample,"_",fixed=TRUE),function(x)x[1])))
    df <- df %>%
        mutate(Time = ifelse(str_ends(Sample, "P"), "Post", "Pre")) %>%
        mutate(id = str_remove(Sample, "P$")) %>%
        mutate(Response = ifelse(id %in% c("SLE-ITP_5","SLE-ITP_9","SLE-ITP_10"),"CR","PR"))
    df$Response <- factor(df$Response,levels = c("CR","PR"))
    df$Time <- factor(df$Time,levels = c("Pre","Post"))
    ct_levels <- c("Pre ProB","ProB","Large PreB","Small PreB","Immature B","Transitional B","Naive B",'ACB',"USM B",'SM B',"AtM","ASC")
    df$cell_type_r1 <- factor(df$cell_type_r1,levels = ct_levels)
    df$Clone_state <- factor(df$Clone_state ,levels = c("Non-VDJ", "No Clonal","Clonal = 2", "Clonal >= 3" ))
    df$id <- factor(df$id,levels = c("SLE-ITP_5"
        , "SLE-ITP_6","SLE-ITP_8",
        "SLE-ITP_9","SLE-ITP_10","SLE-ITP_11"))
    df <- df %>% group_by(Time,Response,Clone_state,cell_type_r1) %>% summarise(n=mean(Proportion))

    df %>% 
        filter(cell_type_r1 %in% c("USM B",'SM B',"AtM","ASC")) %>%
        ggplot( aes(x =Time, y = n, fill = Clone_state)) + 
        geom_bar(position = "fill",stat = "identity") + #position="stack" gives numbers
        scale_fill_manual("Clone state", values = c("#E1E1E1","#A0CBE8","#00366C")) +# "#F9F9F9""black",
        theme_classic()+ 
        facet_grid(~cell_type_r1+Response) + 
        labs(y="Distribution of clone status (%)",x="") + 
        theme(
            legend.text= element_text(color="black", size=10), #family = "Arial", 
            axis.text.x = element_text(angle=90, vjust=0.6,size=12,hjust = 1, color="black"),
            axis.title = element_text(size = 12), #family = "Arial",
            # aspect.ratio = 3/1,
            axis.text.y = element_text(size = 12, color="black"), #family = "Arial",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background= element_blank(),
            axis.line = element_line(colour=NA),
            plot.background = element_rect(colour = NA),
            strip.background = element_blank(),
            strip.text  = element_text(size = 12,color="black") #family = "Arial", 
            # panel.background = element_rect(colour = NA)
    ) 
     ggsave("./RNA/BMMC/figures/cart/ASC-clone-state-prepost-group.png",width = 7,height = 4)
}

# ---------------------------------------
#       Figure 6E
# ---------------------------------------
{
for (ct in unique(filter_all_df_$cell_type_r1)){
  print(ct)
#   ct = "Naive B"
  tmp <- filter_all_df_ %>% 
    filter(cell_type_r1==ct) %>% 
    group_by(Sample) %>%  
    mutate(sample_total=n())
  tmp <- tmp %>% group_by(Sample,c_call) %>% 
    mutate(cell_total=n())
  
  df<- tmp %>% distinct(Sample,Disease, c_call,
                        sample_total, cell_total) %>% mutate(freq = cell_total/sample_total) #group,
  
  df <- df %>% filter(Sample %in% c(
      "SLE-ITP_5P"
      , "SLE-ITP_6P","SLE-ITP_8P",
      "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P",
      "SLE-ITP_5"
      , "SLE-ITP_6","SLE-ITP_8",
      "SLE-ITP_9","SLE-ITP_10","SLE-ITP_11")) %>%
      mutate(Time = ifelse(str_ends(Sample, "P"), "Post", "Pre")) %>%
      mutate(id = str_remove(Sample, "P$")) %>%
      mutate(Response = ifelse(id %in% c("SLE-ITP_5","SLE-ITP_9","SLE-ITP_10"),"CR","PR"))
  # head(df)
  unique(df$Sample)
  df$Response <- factor(df$Response,levels = c("CR","PR"))
  df$Time <- factor(df$Time,levels = c("Pre","Post"))
  p <- ggplot(df, aes(x=Time, y = freq, fill=c_call))+ 
    geom_bar(position="fill", stat="identity", width = 0.8)+
    labs(y="Cell proportions (%)",x="",title=ct) +
    # labs(y="Distribution of clone status (%)",x="") + 
    # facet_grid(~group) + 
    # theme_foundation() +
    facet_grid(~Response) + 
    scale_fill_manual(values=isotype.colors)  +
    guides(fill=guide_legend(title = "Isotype",ncol = 1)) +
    theme(
      legend.text= element_text( color="black", size=10),
      axis.text.x = element_text(angle=90, vjust=0.6,size=12,hjust = 1,color="black"),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12,hjust = 0.5),
      # aspect.ratio = 3/1,
      axis.text.y = element_text(size = 12,color="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background= element_blank(),
      # axis.line = element_line(colour="black"),
      axis.line = element_line(colour=NA),
      # panel.border = element_rect(colour = "black", fill=NA, size=0.8),
      # panel.border = element_blank(),
      plot.background = element_rect(colour = NA),
      strip.background = element_blank(),
      strip.text  = element_text(size = 12,color="black"),
      # panel.background = element_rect(colour = NA)
    )  
  prefix = gsub(" ","",ct)
  p
  ggsave(sprintf("./RNA/BMMC/figures/cart/%s-Isotype-group.pdf",prefix),
         width = 4, height = 4)
}
}


# ---------------------------------------
#       Figure 6H
# ---------------------------------------
library(tidyverse)
library(reshape2)
library(scales)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

df.myeloid.score <- read.csv("./RNA/BMMC/subclustering/Myeloid/myeloid.signature.csv") %>%
  filter(cell_type_r1 == "Mono_C1Q") %>%
  mutate(M1=M1.blood_gene_signature_ucell_score) %>%
  mutate(M2=M2.blood_gene_signature_ucell_score) %>%
  select(Sample,M1,M2) 
df.myeloid.score.summ <-  df.myeloid.score %>% group_by(Sample) %>%
  summarise(across(starts_with("M"), mean, na.rm = TRUE)) %>%
  mutate("M1.M2"=M1/M2) %>%melt()
head(df.myeloid.score.summ)

df<- df.myeloid.score.summ %>%
  filter(Sample %in% c(
    "HC_5","HC_6","HC_7",
      "SLE-ITP_5P",
       "SLE-ITP_6P","SLE-ITP_8P",
      "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P",
      "SLE-ITP_5"
      , "SLE-ITP_6","SLE-ITP_8",
      "SLE-ITP_9","SLE-ITP_10","SLE-ITP_11"
    )) %>%
    mutate(Time = ifelse(str_ends(Sample, "P"), "Post",
     ifelse(str_starts(Sample, "HC"),"HC","Pre"))) %>%
    #  "Pre")) %>%
    mutate(id = str_remove(Sample, "P$")) %>%
    mutate(Response = ifelse(id %in% c("SLE-ITP_5","SLE-ITP_9","SLE-ITP_10"),"CR",
    ifelse(id %in% c( "SLE-ITP_6","SLE-ITP_8","SLE-ITP_11"),"PR","HC") )) %>%
    group_by(Time,Response,variable) %>%
   summarize(across(where(is.numeric),mean)) %>%
   mutate(group=paste(Time,Response,sep=".")) %>%
   as.data.frame() 
plot.df <- df %>%
  as.data.frame() %>%
  select(group,value,variable) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  as.data.frame() 
  # %>%
  # column_to_rownames("group")
dd <- plot.df %>% filter(group %in% c("HC.HC","Post.CR","Pre.CR")) %>%column_to_rownames("group")
dd <- plot.df %>% filter(group %in% c("HC.HC","Post.PR","Pre.PR")) %>%column_to_rownames("group")


dd <- t(scale(dd))
dd[dd<-1] = -0.5
col_fun = colorRamp2(c(-1,0,0.8), c("#4B98C5",  "white","#BB2B34")) 

dd = dd[,c("HC.HC","Pre.CR" ,"Post.CR" )]
colnames(dd) <- c("HC","CR Pre","CR Post")
rownames(dd) <- c("M1","M2","M1/M2")

dd = dd[,c("HC.HC","Pre.PR" ,"Post.PR" )]
colnames(dd) <- c("HC","PR Pre","PR Post")
rownames(dd) <- c("M1","M2","M1/M2")
p <- Heatmap(dd,name=" ",
          column_title = "Mono_C1Q",
          # col = col_fun,
          col = col_fun,
          # col = coul,
          
          border_gp = gpar(col = "black", lty = 1),
          rect_gp = gpar(col = "white", lwd = 0.2),
          cluster_rows = F,
          # cluster_rows = as.dendrogram(row_clust),
          show_row_dend = FALSE,
          # show_column_dend = FALSE,
          cluster_columns = F,
          # column_title = ,
          column_dend_side="bottom",
          column_dend_height = unit(8, "mm"),
          row_dend_width = unit(8, "mm"),
          show_column_names = T,
          # cluster_columns = as.dendrogram(col_clust),
          # cluster_columns = OR.hclust.col$branch,
          column_title_gp = gpar(fontface = "bold",fontsize = 14),
          #ht_opt(heatmap_border = TRUE),
          ht_opt(heatmap_column_names_gp=gpar(fontsize = 12), #fontface = "bold",
                 heatmap_row_names_gp = gpar(fontsize = 14)), #fontface = "bold"
          row_names_side = "left",
          # top_annotation = ha2,
          heatmap_legend_param = (list(at = c(-1,0, 0.8),
                                       labels = c("Min","","Max")
          )))
draw(p)
getwd()