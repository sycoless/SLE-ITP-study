library(tidyverse)
# setwd("")
options(stringsAsFactors = F)
library(rstatix)

# ---------------------------------------
#       Figure 2C
# ---------------------------------------

{   
  dat <- read.csv("./B_pbmc_cellprop.csv") %>%
    column_to_rownames("X") %>% 
    mutate(Sample=rownames(.)) %>% 
    # mutate(Disease=sub("_\\d+$", "", rownames(.))) %>% 
    mutate(
      Disease = str_extract(Sample, "^[^_]+"),
      Source  = str_extract(Sample, "(PBMC|BMMC)"),   # 提取 PBMC 或 BMMC
      Disease = paste(Disease, Source,sep=".")                # 拼接成新列
    ) %>%
      as.data.frame() %>%
      mutate(Disease = str_replace(Disease, "^(SLE-TP|Nr-SLE-TP)", "SLE-ITP")) %>%
      mutate(Disease = str_replace(Disease, "^(SLE-nonTP)", "SLE-NP")) %>%
      mutate(Disease = str_replace(Disease, "(BMMC)$", "BM")) %>%
      mutate(Disease = str_replace(Disease, "(PBMC)$", "PB")) %>%
      filter(!Sample %in% c("SLE-TP_8_BMMC-Post-B","SLE-TP_11_BMMC-Post",
      "SLE-TP_10_BMMC-Post","SLE-TP_9_BMMC-Post","SLE-TP_8_BMMC-Post",
      "SLE-TP_6_BMMC-Post","SLE-TP_5_BMMC-Post"
      ))
    

  df_long <- dat %>%
    pivot_longer(
      cols = -c(Sample, Disease,Source),  
      names_to = "CellType",     
      values_to = "Value"        
    )
  df_long$CellType <- gsub("\\.", " ", df_long$CellType)

  unique(df_long$CellType)
  unique(df_long$Sample)
  unique(df_long$Disease)
  my_comparisons <- list(c("HC.PB", "HC.BM"),c("SLE-ITP.PB", "SLE-ITP.BM"),c("HC.BM", "SLE-NP.BM"),c("HC.PB", "SLE-ITP.PB"),
  c("HC.BM", "SLE-ITP.BM"),c("SLE-NP.BM", "SLE-ITP.BM"))
  # df_long$Disease <- factor(df_long$Disease,levels=c("HC.PBMC","HC.BMMC","SLE-NP.BMMC","SLE-ITP.PBMC" ,"SLE-ITP.BMMC"))
  df_long$Disease <- factor(df_long$Disease,levels=c("HC.PB","HC.BM","SLE-NP.BM","SLE-ITP.PB" ,"SLE-ITP.BM"))
  library(ggpubr)
  
  for(i in unique(df_long$CellType)){
    print(i)
    ggplot(df_long %>% filter(CellType==i), aes(x = Disease, y = Value, fill = Disease)) +
      geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(0.8)) +
      geom_point(position = position_dodge(0.8), size = 2, shape = 21) +
      stat_compare_means(
        aes(group = Disease),
        comparisons=my_comparisons,
        method = "wilcox.test", label = "p.signif",method.args = list(alternative = "two.sided")) +
      theme_classic()+
      labs(title = i, y = "% of B cells", x = "",fill="") +
      # scale_fill_manual(values = c("CR" = "#4B98C5", "PR" = "#BB2B34")) +
      scale_fill_manual(values = c("HC.PB"="#31599B",
      "HC.BM"="#B2CFE5",
      "SLE-NP.BM"="#D97351",
      "SLE-ITP.PB"="#CD7F97", 
      "SLE-ITP.BM"="#8A2A46" )) +
      # scale_y_continuous(expand=c(0,3))+
      scale_y_continuous(expand = expansion(mult = c(0, .1)))+
      # guide_legend(fill=())
      theme(
        # axis.text.x = element_text(angle = 45, hjust = 1)
        axis.text.x = element_text(angle=45,hjust=1,size=12,color = "black"),
        plot.title = element_text( hjust = 0.5,size=14),
        # axis.text.x = element_text(size=10,color = "black"),
        axis.text.y = element_text(size=10,color = "black"),
        legend.position="none"
    )
  ggsave(sprintf("./RNA/BMMC/figures/B_%s_proportion_pbmc.pdf",i),width = 3.5,height = 4)
}
}

# ---------------------------------------
#       Figure 2D
# ---------------------------------------
setwd("./Bcells/BCR/data_aa0.1")
library(tidyverse)
library(alakazam)
library(shazam)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(stats)
library(stringr)
library(ggrepel)

files <- dir('./')
HC_df <- data.frame()
for (file in c("HC_5_BMMC-BCR","HC_6_BMMC-BCR","HC_7_BMMC-BCR",
               'HC_5_PBMC-BCR', 'HC_6_PBMC-BCR','HC_7_PBMC-BCR'
)){ 

  print(file)
  name = str_split(file,'-')[[1]][1]
  print(name)
  tmp <- read.csv(paste0('./',file,'/10X_clone-pass_germ-pass.tsv'),header = T,sep = '\t') %>% 
    mutate(id=name) %>% 
    mutate(group="HC")
  if(dim(tmp)[1]==0){
    HC_df <- tmp
  }else {
    HC_df <- rbind(HC_df,tmp)
  }
}

SLE_df <- data.frame()
for (file in c("SLE-TP_1_BMMC-BCR","SLE-TP_2_BMMC-BCR",
               "SLE-TP_1_PBMC-BCR","SLE-TP_2_PBMC-BCR",
               "SLE-TP_3_BMMC-BCR",
               "SLE-TP_4_BMMC-BCR","SLE-TP_5_BMMC-BCR","SLE-TP_6_BMMC-BCR",
               "SLE-TP_7_BMMC-BCR","SLE-TP_8_BMMC-BCR","SLE-TP_8_BMMC-Post-B-BCR",
               'SLE-TP_9_BMMC-BCR',
               'SLE-TP_10_BMMC-BCR',
               'SLE-TP_11_BMMC-BCR',
               'SLE-TP_8_BMMC-Post-B-BCR',
               'SLE-TP_5_BMMC-Post-BCR',
               'SLE-TP_6_BMMC-Post-BCR',
               'SLE-TP_8_BMMC-Post-BCR',
               'SLE-TP_9_BMMC-Post-BCR',
               'SLE-TP_10_BMMC-Post-BCR',
               'SLE-TP_11_BMMC-Post-BCR',
               'Nr-SLE-TP_1_BMMC-BCR',
               'Nr-SLE-TP_2_BMMC-BCR',
               'Nr-SLE-TP_3_BMMC-BCR'
)){
  print(file)
  name = str_split(file,'-BCR')[[1]][1]
  print(name)
  tmp <-  read.csv(paste0('./',file,'/10X_clone-pass_germ-pass.tsv'),header = T,sep = '\t') %>% 
    mutate(id=name) %>% 
    mutate(group="SLE-ITP")
  if(dim(tmp)[1]==0){
    SLE_df <- tmp
  }else {
    SLE_df <- rbind(SLE_df,tmp)
  }
}

SLEnonTP_df <- data.frame()
for (file in c("SLE-nonTP_1_BMMC-BCR","SLE-nonTP_2_BMMC-BCR","SLE-nonTP_3_BMMC-BCR")){
  print(file)
  name = str_split(file,'-B')[[1]][1]
  tmp <-  read.csv(paste0('./',file,'/10X_clone-pass_germ-pass.tsv'),header = T,sep = '\t') %>% 
    mutate(id=name) %>% 
    mutate(group="SLE-nonTP")
  if(dim(tmp)[1]==0){
    SLEnoTPn_df <- tmp
  }else {
    SLEnonTP_df <- rbind(SLEnonTP_df,tmp)
  }
}

all_df <- rbind(HC_df,SLEnonTP_df,SLE_df) 
all_df <- all_df %>% 
  mutate(barcode = as.character(lapply(
    strsplit(cell_id,"-",fixed=T), function(x) x[1]
  ))) %>% 
  mutate(barcode_id = paste(barcode,id,sep="_"))

b_cell_rna <- read.csv("./RNA/BMMC/tables/BMMC_B.metadata.csv") %>%
  mutate(barcode_id = paste(as.character(lapply(strsplit(barcode_id,"_",fixed=T),function(x) x[1])),id,sep="_")) 
filter_all_df <- inner_join(all_df,b_cell_rna,by="barcode_id")

{
  filter_all_df_ <- filter_all_df %>% 
    mutate(id_clone_id = paste0(Sample,clone_id)) 
  clone_sizes <- countClones(dplyr::filter(filter_all_df_, locus == "IGH"),
                             groups = "Sample",clone="clone_id")
  clone_sizes <- clone_sizes %>%  mutate(id_clone_id = paste0(Sample,clone_id))
  #id_
  clone <- filter_all_df_ %>% 
    select(barcode_id,Disease,cell_type_r1,id_clone_id) %>% 
    left_join(.,clone_sizes,by="id_clone_id")
  
}

{
  clone2 <- clone %>% 
    filter(!Sample %in% c("SLE-ITP_8P","SLE-ITP_8P2","SLE-ITP_5P","SLE-ITP_6P",
     "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P")) %>%
    mutate(Clone_state=ifelse((seq_count==1),"No Clonal",
                              ifelse(( seq_count==2 ),"Clonal = 2",# @seq_count>=2 &
                                     "Clonal >= 3")))
  clone2[which(is.na(clone2$Clone_state)),"Clone_state"] <- "Non-VDJ"
  table(clone2$Clone_state)
  
  df <- clone2 %>% 
    group_by(Sample,cell_type.PBMC,Clone_state,.drop=FALSE) %>% 
    summarise(n=n())
  df <- df %>% 
    group_by(Sample, cell_type.PBMC) %>%
    mutate(Proportion = n / sum(n) * 100) %>%
    ungroup() %>% 
    mutate(Disease=as.character(lapply(strsplit(df$Sample,"_",fixed=TRUE),function(x)x[1])))
  plot.df <- df %>% group_by(Disease,Clone_state,cell_type.PBMC) %>% summarise(n=mean(Proportion))
  plot.df$Disease <- factor(plot.df$Disease,levels = c("HC","SLE-NP","SLE-ITP"))
  ct_levels <- c("Pre ProB","ProB","Large PreB","Small PreB","Immature B","Naive B",'ACB',"USM B",'SM B',"AtM","ASC")
  plot.df$cell_type.PBMC <- factor(plot.df$cell_type.PBMC,levels = ct_levels)
  plot.df$Clone_state <- factor(plot.df$Clone_state ,levels = c("Non-VDJ", "No Clonal","Clonal = 2", "Clonal >= 3"  ))
  plot.df %>% 
  filter(cell_type.PBMC %in% c("Naive B","USM B",'SM B',"AtM","ASC")) %>%
   ggplot( aes(x =Disease, y = n, fill = Clone_state)) + 
    geom_bar(position = "fill",stat = "identity") + #position="stack" gives numbers
    scale_fill_manual("Clone state", values = c("#E1E1E1","#A0CBE8","#00366C")) +# "#F9F9F9"
    theme_classic()+ 
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(~cell_type.PBMC) + 
    labs(y="Distribution of clone status (%)",x="") + 
    theme(
      legend.text= element_text(family = "Arial", color="black", size=10),
      axis.text.x = element_text(angle=90, vjust=0.6,size=12,hjust = 1,family = "Arial", color="black"),
      axis.title = element_text(family = "Arial",size = 12),
      # aspect.ratio = 3/1,
      axis.text.y = element_text(size = 12,family = "Arial", color="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background= element_blank(),
      axis.line = element_line(colour=NA),
      plot.background = element_rect(colour = NA),
      strip.background = element_blank(),
      strip.text  = element_text(size = 12,family = "Arial", color="black"),
    ) 
  
}

# ---------------------------------------
#       Figure 2E
# ---------------------------------------
# load("./Input_data/celltype_2.colors.rda")
filter_all_df_ <- filter_all_df
isotype.colors <- c(IGHA1="#fec44f",IGHA2="#fe9929",IGHG1="#c6dbef",IGHG2="#6baed6",IGHG3="#2171b5",IGHG4="#08306b",
                    IGHD="#807dba",IGHE="#D51F26",noBCR="#d9d9d9",IGHM="#fa9fb5")
outdir = "./RNA/BMMC/subclustering/Bcells/figures/"
for (ct in unique(filter_all_df_$cell_type_r1)){
  print(ct)
  tmp <- filter_all_df_ %>% 
    filter(cell_type_r1==ct) %>% 
    group_by(Disease) %>% 
    mutate(sample_total=n())
  
  tmp <- tmp %>% group_by(Disease,c_call) %>% 
    mutate(cell_total=n())
  
  df<- tmp %>% distinct(Disease, c_call,
                        sample_total, cell_total) %>% mutate(freq = cell_total/sample_total) #group,
  
  df<-df %>% group_by(Disease, c_call) %>% 
    summarise(average_props = mean(freq))
  
  df$Disease <- factor(df$Disease,levels = c("HC","SLE-NP","SLE-ITP"))
  p <- ggplot(df, aes(x=Disease, y = average_props, fill=c_call))+ 
    geom_bar(position="fill", stat="identity", width = 0.8)+
    labs(y="Cell proportions (%)",x="",title=ct) +
    # labs(y="Distribution of clone status (%)",x="") + 
    # facet_grid(~group) + 
    # theme_foundation() +
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
  ggsave(paste0(outdir,prefix,"-BMMC_changeO-isotype-Disease.png"),width = 3.5, height = 4)
}


# ---------------------------------------
#       Figure 2F
# ---------------------------------------
filter_all_df_ <- filter_all_df
db_obs_all <- observedMutations(filter_all_df_, sequenceColumn="sequence_alignment",
                                germlineColumn="germline_alignment_d_mask",
                                regionDefinition=NULL,
                                frequency=T,
                                combine = T, nproc=5)
colnames(filter_all_df_)
unique(filter_all_df_$Sample)
table(filter_all_df_$cell_type_r1)
library(plyr)

db_obs_all %<>% filter(c_call != '')

dat <- db_obs_all %>% 
  filter(cell_type_r1 %in% c("ASC")) %>%
  filter(!Sample %in% c("SLE-ITP_8P2","SLE-ITP_8P","SLE-ITP_5P","SLE-ITP_6P",
  "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P")) %>%
  group_by(Sample) %>% 
  mutate(mean_mu_freq=mean(mu_freq)) %>% 
  distinct(Sample,.keep_all = T) %>% 
  select(c("Sample","Disease","mean_mu_freq"))
cp <- c( "HC"="#B2CFE5","ITP"= "#9E6D7F" ,"SLE-ITP"="#8A2A46","SLE-NP"="#D97351")
library(ggpubr)
dat$Disease <- factor(dat$Disease,levels=c("HC","SLE-NP","SLE-ITP"))
ggplot(dat,aes(x=Disease,y=mean_mu_freq*100,fill=Disease))+geom_boxplot()+
  theme_classic()+
  geom_jitter(width = 0.2, alpha = 0.8) +
  stat_compare_means(aes(label=..p.signif..),comparisons = list(c("HC","SLE-ITP"),
                                                                c("HC","SLE-NP"),
                                                                c("SLE-NP","SLE-ITP")),
                     method="wilcox.test")+
  ylab("Mutation frquency (%)")+xlab('')+
  labs(title = "ASC") +
  scale_fill_manual(values=cp)+
  scale_y_continuous(expand = c(0,0.5)) +
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle=90,hjust=1,size=10,color = "black"),
    axis.text.y = element_text(size=10,color = "black"),
    legend.position = "None"
  )

# ---------------------------------------
#       Figure 2G
# ---------------------------------------
library(ggpubr)
library(rstatix)
library(DESeq2)

runVDJ.diffUsage<-function(df,v.call="germline_v_call",
                           sample.label="samples",group.label="group",
                           GroupA,GroupB,reference,
                           min.cells=10){
  sample.size=table(df[,sample.label] %>% as.vector()) %>% as.data.frame()
  rownames(sample.size)<-sample.size$Var1
  sample.size<-sample.size[c(GroupA,GroupB,reference),]
  keep.samples<-sample.size[sample.size$Freq>=min.cells,]
  keep.samples$group<-NA
  keep.samples[GroupA,]$group<-"GroupA"
  keep.samples[GroupB,]$group<-"GroupB"
  keep.samples[reference,]$group<-"reference"
  
  
  new<-df
  new$IGHV<-NA
  i=1
  for(i in 1:dim(new)[1]){
    new[,v.call][i] %>%as.vector() %>%strsplit(split='\\*') ->tmp
    new[i,"IGHV"]<-tmp[[1]][1]
  }
  Usage<-prop.table(table(new$IGHV,new[,sample.label]))
  for(j in 1:dim(Usage)[2]){
    Usage[,j]<-Usage[,j]/sum(Usage[,j])
  }
  Usage <- as.data.frame.array(t(Usage[,rownames(keep.samples)]))
  Usage$group<-keep.samples$group
  Usage[,-dim(Usage)[2]]<-Usage[,-dim(Usage)[2]]*1000 +1
  
  condition <- factor(Usage$group)
  colData <- data.frame(row.names=rownames(Usage), condition)
  Usage<-floor(Usage[,-dim(Usage)[2]])
  Usage<-t(Usage) %>% as.data.frame()
  
  dds <- DESeqDataSetFromMatrix(countData = Usage,
                                colData = colData, 
                                design= ~condition)
  dds2 <- DESeq(dds)
  res1 <- results(dds2, contrast=c("condition","GroupA","reference")) %>% as.data.frame()
  res2 <- results(dds2, contrast=c("condition","GroupB","reference"))%>% as.data.frame()
  
  
  data1<-na.omit(res1)
  data2<-na.omit(res2)
  
  
  m<-Usage
  z<-m %>% apply(2,sum)
  for(i in 1:dim(m)[2]){
    m[,i]<-m[,i] / as.numeric(z[i])
  }
  n<-m %>% apply(1,median) %>% as.data.frame()
  y<-intersect(data1 %>% rownames(),data2%>% rownames())
  k<-n[y,]
  
  data<-data.frame(X1=data1[y,]$log2FoldChange,X2=data2[y,]$log2FoldChange,row.names = y,Size=k)
  data$ID<-rownames(data)
  library(ggrepel)
  data$color<-"0"
  data[data$X1>0 &data$X2>0,]$color<-"1"
  data[data$X1<0 &data$X2>0,]$color<-"2"
  data[data$X1<0 &data$X2<0,]$color<-"3"
  data[data$X1>0 &data$X2<0,]$color<-"4"
  P1<-ggplot(data,mapping = aes(x=X1,y=X2,size=Size,color=color))+geom_point()+
    # DimPlot_theme+
    geom_text_repel(data=data[data$X1 >0 & data$X2 < 0,],aes(x=X1,y=X2,label=ID))+
    geom_text_repel(data=data[data$X1 >0 & data$X2 > 0,],aes(x=X1,y=X2,label=ID))+
    geom_hline(aes(yintercept=0))+
    geom_vline(aes(xintercept=0))+scale_color_manual(values = c("black","red","yellow","blue","green"))
  return(list("data"=data,plot=P1,"data1"=data1,"data2"=data2))
}

runVDJ.diffUsage<-function(df,v.call="germline_v_call",
                           sample.label="samples",group.label="group",
                           GroupA,GroupB,reference,
                           min.cells=10){
  sample.size=table(df[,sample.label] %>% as.vector()) %>% as.data.frame()
  rownames(sample.size)<-sample.size$Var1
  sample.size<-sample.size[c(GroupA,GroupB,reference),]
  keep.samples<-sample.size[sample.size$Freq>=min.cells,]
  keep.samples$group<-NA
  keep.samples[GroupA,]$group<-"GroupA"
  keep.samples[GroupB,]$group<-"GroupB"
  keep.samples[reference,]$group<-"reference"
  
  
  new<-df
  new$IGHV<-NA
  i=1
  for(i in 1:dim(new)[1]){
    new[,v.call][i] %>%as.vector() %>%strsplit(split='\\*') ->tmp
    new[i,"IGHV"]<-tmp[[1]][1]
  }
  Usage<-prop.table(table(new$IGHV,new[,sample.label]))
  for(j in 1:dim(Usage)[2]){
    Usage[,j]<-Usage[,j]/sum(Usage[,j])
  }
  Usage <- as.data.frame.array(t(Usage[,rownames(keep.samples)]))
  Usage$group<-keep.samples$group
  Usage[,-dim(Usage)[2]]<-Usage[,-dim(Usage)[2]]*1000 +1
  
  condition <- factor(Usage$group)
  colData <- data.frame(row.names=rownames(Usage), condition)
  Usage<-floor(Usage[,-dim(Usage)[2]])
  Usage<-t(Usage) %>% as.data.frame()
  
  dds <- DESeqDataSetFromMatrix(countData = Usage,
                                colData = colData, 
                                design= ~condition)
  dds2 <- DESeq(dds)
  res1 <- results(dds2, contrast=c("condition","GroupA","reference")) %>% as.data.frame()
  res2 <- results(dds2, contrast=c("condition","GroupA","GroupB"))%>% as.data.frame()
  
  
  data1<-na.omit(res1)
  data2<-na.omit(res2)
  
  
  m<-Usage
  z<-m %>% apply(2,sum)
  for(i in 1:dim(m)[2]){
    m[,i]<-m[,i] / as.numeric(z[i])
  }
  n<-m %>% apply(1,median) %>% as.data.frame()
  y<-intersect(data1 %>% rownames(),data2%>% rownames())
  k<-n[y,]
  
  data<-data.frame(X1=data1[y,]$log2FoldChange,X2=data2[y,]$log2FoldChange,row.names = y,Size=k)
  data$ID<-rownames(data)
  library(ggrepel)
  data$color<-"0"
  data[data$X1>0 &data$X2>0,]$color<-"1"
  data[data$X1<0 &data$X2>0,]$color<-"2"
  data[data$X1<0 &data$X2<0,]$color<-"3"
  data[data$X1>0 &data$X2<0,]$color<-"4"
  P1<-ggplot(data,mapping = aes(x=X1,y=X2,size=Size,color=color))+geom_point()+
    # DimPlot_theme+
    geom_text_repel(data=data[data$X1 >0 & data$X2 < 0,],aes(x=X1,y=X2,label=ID))+
    geom_text_repel(data=data[data$X1 >0 & data$X2 > 0,],aes(x=X1,y=X2,label=ID))+
    geom_hline(aes(yintercept=0))+
    geom_vline(aes(xintercept=0))+scale_color_manual(values = c("black","red","yellow","blue","green"))
  return(list("data"=data,plot=P1,"data1"=data1,"data2"=data2))
}

plot.df <- filter_all_df %>% 
     filter(cell_type_r1 %in% c("ASC")) %>%  filter(!Sample %in% c("SLE-ITP_8P","SLE-ITP_8P2",
                                                 "SLE-ITP_5P","SLE-ITP_6P",
                                                 "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P"))
res <- runVDJ.diffUsage(plot.df,v.call="v_call_10x",
                        sample.label="Sample",group.label="Disease",
                        GroupA=c(
                          "SLE-ITP_1",
                          "SLE-ITP_2",
                          "SLE-ITP_3",
                          "SLE-ITP_4","SLE-ITP_5", 
                          "SLE-ITP_6" ,
                          "SLE-ITP_7","SLE-ITP_8",   "SLE-ITP_9",
                          "SLE-ITP_10",
                          "SLE-ITP_11",
                          "SLE-ITP_12" , "SLE-ITP_13" ,
                          "SLE-ITP_14"
                        ),GroupB=c("SLE-NP_1","SLE-NP_2","SLE-NP_3"),reference=c("HC_5" ,"HC_6" ,"HC_7"),min.cells=10)
df = res$data
data = res$data
p <- ggplot(df,aes(x=X1,y=X2,size=Size*100,color=color))+geom_point()+
  # DimPlot_theme+
  geom_text_repel(data=data[data$X1 >0 & data$X2 < 0,],aes(x=X1,y=X2,label=ID),size=2,min.segment.length=0.2)+
  geom_text_repel(data=data[data$X1 >0 & data$X2 > 0,],aes(x=X1,y=X2,label=ID),size=2,min.segment.length=0.2)+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept=0))
p + theme_bw() + 
  xlab(bquote("Log"[2]~"(Fold Change)")) +
  ylab(bquote("Log"[2]~"(Fold Change)")) +
  theme(
    axis.text.x = element_text(size = 12,color = "black"),
    axis.text.y = element_text(size = 12,color = "black"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    panel.border = element_blank(),
    # legend.position = "none"
  ) +
  theme(
    panel.border = element_rect(color='black', fill=NA, size=1),
    panel.grid.major = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position='bottom'
  ) +
  scale_color_manual(values = c("#a05528"  ,"#D97351" ,"#bac4d0" ,"#916ba6")) + 
  guides(size=guide_legend(title="Usage (%)"),color="none") +
  scale_y_continuous(limits=c(-5,5))


# ---------------------------------------
#       Figure 2H
# ---------------------------------------
plot.df <- filter_all_df %>% 
     filter(cell_type_r1 %in% c("ASC")) 
table(plot.df$Sample)

p <- runVDJ.Usage(plot.df, v.call="germline_v_call",sample.label="Sample",group.label="Disease",min.cells=10)
df2 <- p$df 
df2 <- df2 %>% filter(!Samples %in% c(
    "SLE-ITP_8P2","SLE-ITP_8P",
     "SLE-ITP_5P","SLE-ITP_6P",
     "SLE-ITP_9P","SLE-ITP_10P","SLE-ITP_11P"
     ))
clinial <- read.csv("./data/clinical/Platelet_meta.csv")
colnames(clinial)[1] <- "Samples"
df_ <- df2  %>% 
    left_join(.,clinial,by="Samples") %>% 
    filter(group=="SLE-ITP") 

# df_ <- df2  %>% 
#     left_join(.,clinial,by="Samples") %>% 
#     na.omit()
#     filter(group %in% c("SLE-ITP","HC")) 

result <- df_ %>%
  group_by(IGHV) %>%
  summarize(
    cor = cor(Usage, log(platelet_count), method = "spearman"),
    p_value = cor.test(Usage, platelet_count, method = "spearman")$p.value
  )

result <- df_ %>%
  group_by(IGHV) %>%
  summarize(
    cor = cor(Usage, log(platelet_count), method = "spearman"),
    p_value = cor.test(Usage, platelet_count, method = "spearman")$p.value
  )
# print(result  %>% arrange(desc(cor)),n=20)
print(result  %>% arrange(cor),n=20)

p <- ggplot(df_ %>% filter(IGHV=="IGHV3-9"), aes(x = Usage, y = platelet_count)) +
  geom_point(shape=21,aes(fill=group),size=3, color="white")+
  scale_fill_manual(name="group",values = c("#374E55FF","#DF8F44FF","#00A1D5FF","#B24745FF"))+
  geom_smooth(se = F,method = "lm") + 
  stat_cor(method = "spearman") +
  # stat_cor(method = "pearson") +
  theme_bw() 
p

# ---------------------------------------
#       Figure 2I
# ---------------------------------------
library(Startrac)
library(tidyverse)
library(reshape2)
library(tictoc)
library(tidyverse)
options(stringsAsFactors = F)
setwd("./RNA/BMMC/subclustering/Bcells/BCR")

BCR.data <- read.delim("./RNA/BMMC/subclustering/Bcells/pbmc/BCR/changeO-Final.csv")
# head(BCR.data)
BCR.data$majorCluster <- as.character(BCR.data$cell_type_r1)
BCR.data$clone.id <- BCR.data$clone_id #changeo_
ncell.patient.cluster <- sort(unclass(table(sprintf("%s.%s",BCR.data$Sample,BCR.data$majorCluster))))
BCR.data <- BCR.data[ncell.patient.cluster[sprintf("%s.%s",BCR.data$Sample,BCR.data$majorCluster)]>=5,]
table(BCR.data$Sample)
dim(BCR.data)

BCR.data$Cell_Name=BCR.data$barcode_id
BCR.data$clone.id=BCR.data$clone_id
BCR.data$patient=BCR.data$id.x
BCR.data$majorCluster=BCR.data$cell_type_r1
BCR.data$loc="BMMC"

dir.startrac <- "./Additional_data"
resultsfile <- "B.out.startrac_changO.rds"
##all type
if(!file.exists(sprintf("%s/%s",dir.startrac,resultsfile))){
  tic("Startrac.run")
  out <- Startrac.run(BCR.data, proj="panC",verbose=T,cores=3,n.perm=50)
  toc()
  out2 <- list()
  out2[['proj']] <- out@proj
  out2[['cluster.data']] <- out@cluster.data
  out2[['cluster.sig.data']] <- out@cluster.sig.data
  out2[['pIndex.migr']] <- out@pIndex.migr
  out2[['pIndex.tran']] <- out@pIndex.tran
  out2[['pIndex.sig.migr']] <- out@pIndex.sig.migr
  out2[['pIndex.sig.tran']] <- out@pIndex.sig.tran
  saveRDS(out2,sprintf("%s/%s",dir.startrac,resultsfile))
}

patient.results <- "B.out.startrac_changO_BMMC.patient.rds"
patient.vec <- unique(BCR.data$Sample)
if(!file.exists(sprintf("%s/%s",dir.startrac,patient.results))){
  res.bypatient <- lapply(patient.vec,function(x){
    dat <- BCR.data %>%
      filter(Sample == x)
    Startrac.run(dat,
                 proj=x,verbose=F,cores=4,n.perm=NULL)
  })
  names(res.bypatient) <- patient.vec
  saveRDS(res.bypatient,sprintf("%s/%s",dir.startrac,patient.results))
}

tese <- lapply(names(out),function(x) a=out[[x]]@pIndex.sig.tran)
pIndex.sig.tran <- do.call("rbind",tese) %>% 
  filter(majorCluster %in% c("ASC")) %>% 
  select(c("aid","index","value")) %>% 
  pivot_wider(names_from = "aid",values_from = "value") %>% 
  column_to_rownames("index") 

pIndex.sig.tran[is.na(pIndex.sig.tran)]=0
z.max <- 1.5
dat.plot.mtx <- t(scale(t(pIndex.sig.tran)))
dat.plot.mtx[ dat.plot.mtx > z.max ] <- z.max
dat.plot.mtx[ dat.plot.mtx < -z.max ] <- -z.max


col_fun = colorRamp2(c(-1,0,1), c("#4B98C5",  "white","#BB2B34")) 
dat.plot.mtx <- dat.plot.mtx[c("Immature B","Naive B","ACB","AtM","USM B","SM B"),]

colnames(dat.plot.mtx) <- gsub("SLE-ITP_8P2","SLE-ITP_8P",colnames(dat.plot.mtx))
x = Heatmap(dat.plot.mtx,name="pTrans",col = col_fun,
        border_gp = gpar(col = "black", lty = 1),
        rect_gp = gpar(col = "white", lwd = 0.2),
        cluster_rows = F,
        # cluster_rows = as.dendrogram(row_clust),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        cluster_columns = T,
        # column_title = ,
        column_dend_side="bottom",
        column_dend_height = unit(8, "mm"),
        row_dend_width = unit(8, "mm"),
        # cluster_columns = as.dendrogram(col_clust),
        # cluster_columns = OR.hclust.col$branch,
        column_title_gp = gpar(fontface = "bold",fontsize = 14),
        #ht_opt(heatmap_border = TRUE),
        ht_opt(heatmap_column_names_gp=gpar(fontsize = 10), #fontface = "bold",
               heatmap_row_names_gp = gpar(fontsize = 10)), #fontface = "bold"
        row_names_side = "left")


library(RColorBrewer)
coul <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(30))
