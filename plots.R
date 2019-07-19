library(ggplot2)
library(scales)
library(stringr)
library(gridExtra)
library(gtable)
library(grid)

# Composition of clusters/tissues

clus_tis <- read.csv("clus_tis_df.csv")

tissue_cols <- c('#d62728', # blood
                 '#279e68', # brain
                 '#964B00', # colon
                 '#ff00ef', # decidua
                 '#aa40fc', # fetal liver
                 '#E8E848', # fetal skin
                 '#aec7e8', # kidney
                 '#1406D6', # lung
                 '#b5bd61', # microglia
                 '#65069C', #mLN            
                 '#ff7f0e', # pancreas
                 '#ff9896', # placenta
                 '#00fff3', # prostate
                 '#00ff89', # testis
                 '#F7C8D8') # upper airway

Clusters <- c()
Tissues <- c()
Percents <- c()
for(i  in 0:length(unique(clus_tis$clusters))){
  tmp <- subset(clus_tis, clus_tis$clusters==i)
  p <- c()
  for(j in unique(tmp$tissues)){
    p <- c(p, (length(tmp$tissues[tmp$tissues==j])/length(tmp$tissues))*100)}
  Clusters <- c(Clusters, rep(paste("Cluster", i), length(p)))
  Tissues <- c(Tissues,as.character(unique(tmp$tissues) ))
  Percents <- c(Percents, p)
}


df <- data.frame(Clusters, Tissues, Percents)
df$Clusters <- factor(df$Clusters, levels = unique(df$Clusters))

png("cluster_comp.png", width = 1200, height = 650)
ggplot(df, aes(x=Clusters, y = Percents, fill=Tissues)) +
  geom_col(width=0.9, col='black', size = 0.02) +
  scale_y_continuous(limits = c(0,101), breaks = seq(25, 100,25), labels = paste(seq(25, 100,25), "%"), expand = c(0,0)) +
  labs(x = "", y="") +
  scale_fill_manual(values= tissue_cols, name= NULL) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme_classic() +
  theme(axis.text.x =element_text(size = 28, angle = 90, hjust = 1, color='black'),
        axis.text.y = element_text(size = 24, color='black'),
        axis.ticks.x = element_line(),
        axis.ticks=element_blank(), 
        legend.position = "top",
        legend.text = element_text(size=26.5, color = 'black'),
        legend.key.size = unit(1.5,"line")) 
dev.off()


cluster_cols <- c('#FF230A', #1
                  '#BA2223', #2
                  '#993636', #3
                  '#7D4D43', #4
                  '#DC5BFC', #5
                  '#DBD767', #6
                  '#98df8a', #7
                  '#279e68', #8
                  '#17becf', #9
                  '#aec7e8', #10
                  '#005EFF', #11
                  '#FF70DE', #12
                  '#ffbb78', #13
                  '#c5b0d5', #14
                  '#65069C', #15
                  '#f7b6d2') #16 

Clusters <- c()
Tissues <- c()
Percents <- c()

for(i  in unique(clus_tis$tissues)){
  tmp <- subset(clus_tis, clus_tis$tissues==i)
  p <- c()
  for(j in unique(tmp$clusters)){
    p <- c(p, (length(tmp$clusters[tmp$clusters==j])/length(tmp$clusters))*100)}
  Tissues <- c(Tissues, rep(i, length(p)))
  Clusters <- c(Clusters, unique(tmp$clusters) )
  Percents <- c(Percents, p)
}


df <- data.frame(Clusters, Tissues, Percents)
df$Tissues <- factor(df$Tissues, levels = sort(unique(df$Tissues)))
df$Clusters <- factor(df$Clusters, levels = seq(1,19))

png("tissue_comp.png", width = 1200, height = 650)
ggplot(df, aes(x=Tissues, y = Percents, fill=Clusters)) +
  geom_col(width=0.9, col='black', size = 0.02) +
  scale_y_continuous(limits = c(0,101), breaks = seq(25, 100,25), labels = paste(seq(25, 100,25), "%"), expand = c(0,0)) +
  labs(x = "", y="") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  scale_fill_manual(values= cluster_cols, name = NULL, labels = paste("Cluster", seq(1,16))) +
  theme_classic() +
  theme(axis.text.x =element_text(size = 28, angle = 90, hjust=1, color='black'),
        axis.text.y = element_text(size = 24, color='black'),
        axis.ticks.x = element_line(),
        axis.ticks=element_blank(), 
        legend.position = "top",
        legend.text = element_text(size=24, color = 'black'),
        legend.key.size = unit(1.5,"line")) 
dev.off()

-----

# Functional enrichment plots

  # Manhattan plot of tissue-specific GO terms (Fig 7A)
  
# read in list of selected terms that I want to include
sig_BPs_tiss <- read.csv("~/Documents/Project/GO-BPs.csv")
sig_BPs <- unlist(sig_BPs_tiss, use.names = F)
sig_BPs <- as.character(unique(sig_BPs_tiss))

# read in tissue DEGs and filter non-significant results
tissue_degs <- read.csv('de_genes_tissues.csv', sep = '\t')
tissue_degs <- subset(tissue_degs, abs(tissue_degs$avg_logFC >= 1 & tissue_degs$p_val_adj <= 0.05))

all_gos <- c()
tissue <- c()
p_vals <- c()
recall <- c()
for(i in tissue_degs$cluster){
  filename <- paste( i, "_enrichment_BP.csv", sep = "")
  csv <- read.csv(filename, header = T, sep = ",")
  all_gos <- c(all_gos,as.character(csv[,7]) )
  all_gos <- trimws(all_gos)
  p_vals <- c(p_vals, csv[,2])
  tissue <- c(tissue, rep(i, nrow(csv)))
  recall <- c(recall, csv[,6])
}
counts <- c()
for(i in all_gos){
  counts <- c(counts, as.integer(length(all_gos[all_gos==i])))
}

df1 <- data.frame(tissue,p_vals, all_gos, recall,counts, stringsAsFactors = F)
df1 <- subset(df1, df1$all_gos %in% sig_BPs)
df1$all_gos <- factor(df1$all_gos, levels = rev(unique(df1$all_gos)))
df1$all_gos2 <- str_wrap(df1$all_gos, width = 35)
df1$cat <- df1$tissue
df1$t_c <- "t"

all_gos <- c()
cluster <- c()
p_vals <- c()
recall <- c()
for(i in 0:15){
  filename <- paste("GO_BP_clusters/", i, "_enrichment_BP2.csv", sep = "")
  csv <- read.csv(filename, header = T, sep = ",")
  all_gos <- c(all_gos,as.character(csv[,7]) )
  all_gos <- trimws(all_gos)
  p_vals <- c(p_vals, csv[,2])
  cluster <- c(cluster, rep(i+1, nrow(csv)))
  recall <- c(recall, csv[,6])
}
counts <- c()
for(i in all_gos){
  counts <- c(counts, as.integer(length(all_gos[all_gos==i])))
}
df2 <- data.frame(cluster,p_vals, all_gos, recall,counts, stringsAsFactors = F)
df2 <- subset(df2, df2$all_gos %in% sig_BPs)
df2$cluster <- factor(df2$cluster, levels = unique(df2$cluster))
df2$all_gos <- factor(df2$all_gos, levels = rev(unique(df2$all_gos)))
df2$all_gos2 <- str_wrap(df2$all_gos, width = 35)
df2$cat <- paste("Cluster", df2$cluster)
df2$t_c <- "c"


df3 <- merge(df, df2, by = c("cat", "all_gos", "p_vals", "counts", "recall", "t_c", "all_gos2"), all = T)
df3$cat <- factor(df3$cat, levels = c(unique(df2$cat), unique(df$cat))) 

p1 <- ggplot(df3, aes(x = df3$all_gos2, y = -log10(df3$p_vals))) +
  geom_point(aes(fill = df3$cat, size = df3$recall, alpha=0.9975, shape = df3$t_c)) +
  scale_size(range = c(2,17)) +
  scale_shape_manual(values = c(22, 21)) +
  scale_fill_manual(values = c(cluster_cols, tissue_cols)) +
  guides(alpha = F, shape=F, fill = guide_legend(nrow=2,byrow=TRUE, override.aes = list(size=7, shape = c(rep(22,16), rep(21,15))), title = NULL), size = F) +
  labs(title = "", x = "", y = "-log10(p value)") +
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5, size = 20, color = 'black'),
        axis.text.y = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size=18),
        legend.text=element_text(size=18, color='black'),
        legend.position = "top")
leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box")

p2 <-  ggplot(df3, aes(x = df3$all_gos2, y = -log10(df3$p_vals))) +
  geom_point(aes(fill = df3$cat, size = df3$recall, alpha=0.9975, shape = df3$t_c)) +
  scale_size(range = c(2,17)) +
  scale_shape_manual(values = c(22, 21)) +
  scale_fill_manual(values = c(cluster_cols, tissue_cols)) +
  guides(alpha = F, fill = F, shape = F, size = guide_legend(title = "Recall", title.hjust = 0.5, override.aes = list(shape = 21))) +
  labs(title = "", x = "", y = "-log10(p value)") +
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5, size = 20, color = 'black'),
        axis.text.y = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size=18),
        legend.text=element_text(size=18, color='black'),
        legend.title = element_text(size = 18),
        legend.justification = c(1,0.92))
leg2 <- gtable_filter(ggplot_gtable(ggplot_build(p2)), "guide-box") 

plot <- ggplot(df3, aes(x = df3$all_gos2, y = -log10(df3$p_vals))) +
  geom_point(aes(fill = df3$cat, size = df3$recall, alpha=0.9975, shape = df3$t_c)) +
  scale_size(range = c(2,17)) +
  scale_shape_manual(values = c(22, 21)) +
  scale_fill_manual(values = c(cluster_cols, tissue_cols)) +
  guides(alpha = F, fill = F, size = F, shape = F) +
  labs(title = "", x = "", y = "-log10(p value)") +
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5, size = 20, color = 'black'),
        axis.text.y = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size=18),
        legend.text=element_text(size=18, color='black'),
        legend.position = "top")

plotNew <- arrangeGrob(leg1, plot, 
                       heights = unit.c(leg1$height, unit(1, "npc") - leg1$height), ncol = 1)
plotNew <- arrangeGrob(plotNew, leg2,
                       widths = unit.c(unit(1, "npc") - leg2$width, leg2$width), nrow = 1)


png("BP_enrichment.png", width = 1900, height = 800)
grid.draw(plotNew)
dev.off()

# Cytokine signaling matrixplot (Fig 7B)

# read in list of genes involved in cytokine signaling
cyto_genes <- read.csv("cyto_genes.csv", sep='\t', header=F)

cyto_tissues <- c('Blood', 'Brain', 'Colon', 'Decidua', 'Microglia', 'mLN', 'Prostate', 'Testis', 'Upper airway')
cyto_df <- subset(sig_tissue_degs, sig_tissue_degs$cluster %in% cyto_tissues)
cyto_df <- subset(cyto_df, cyto_df$gene %in% cyto_genes)
cyto_df$gene <- factor(cyto_df$gene, levels = rev(unique(cyto_df$gene)))

colfunc <- colorRampPalette(c("#2D0059", "#4C0099", "#FF6A4D",  "#FF7221", "#FFFF99"))
my_pal <- colfunc(100)

cyto_df$cluster <- factor(cyto_df$cluster, levels = rev(levels(cyto_df$cluster)))
cyto_df$gene <- factor(cyto_df$gene, levels = rev(levels(cyto_df$gene)))

png("cyto_signaling_matrixplot.png", height= 350, width=1600)
ggplot(cyto_df, aes(x=cluster, y=gene))+
  geom_tile(aes(fill = avg_logFC), width=1) +
  labs(fill = "logFC", x ="", y="", title = "Cytokine Signaling in Immune system") +
  scale_fill_gradientn(colours = my_pal) +
  guides( barheight = 10000) +
  #theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size = 14, color = 'black'),
        axis.text.y = element_text(size = 16, colour = 'black'),
        plot.title = element_text(size = 20, hjust = 0.5, color = 'black' ),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        legend.text=element_text(size=12, color='black'),
        legend.title = element_text(size = 12, color='black'),
        legend.key.height = unit(2,"line"),
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "black",
                                        colour = "black")) +
  coord_flip()
dev.off()


