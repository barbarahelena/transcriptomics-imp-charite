## DESeq2 plot

## Libraries
library(ComplexHeatmap)
library(tidyverse)
library(ggrepel)
library(dplyr)

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.0), hjust = 0.5),
            text = element_text(family = 'Helvetica'),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            # legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing  = unit(0, "cm"),
            # legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

} 

heatmap_plot <- function(df, group = "unequal", number = 10) {
    if (group == "unequal")
        df3 <- df %>% group_by(group) %>% 
            arrange(padj) %>% 
            dplyr::slice(1:number) %>% 
            ungroup(.) %>% 
            arrange(log2FoldChange) %>% 
            dplyr::select(Gene.name, padj, contains("ImP"), contains("Ctrl"))
    else if (group == "equal")
        df3 <- df %>%
            arrange(padj) %>% 
            dplyr::slice(1:number) %>% 
            arrange(log2FoldChange) %>% 
            dplyr::select(Gene.name, padj, contains("ImP"), contains("Ctrl"))
    else if (group == "sig")
        df3 <- df %>% 
            arrange(padj) %>% 
            dplyr::filter(padj<0.058) %>% 
            arrange(log2FoldChange) %>% 
            dplyr::select(Gene.name, padj, contains("ImP"), contains("Ctrl"))
    
    df4 <- df3 %>% 
        pivot_longer(., cols = c(3:8), 
                     names_to = "condition",
                     names_prefix = "X12h.",
                     values_to = "count") %>% 
        mutate(condition = fct_recode(condition,
                                      'Control (1)' = 'Ctrl.1', 'Control (2)' = 'Ctrl.2','Control (3)' =  'Ctrl.3',
                                      'ImP (1)' = '100.nM.ImP.1','ImP (2)' = '100.nM.ImP.2','ImP (3)' = '100.nM.ImP.3')
        ) %>% 
        pivot_wider(., id_cols = c(1,2), names_from = "condition", values_from = "count",
                    names_sort = TRUE)
    
    rownames(df4) <- df4$Gene.name
    df3_counts <- as.matrix(df4[,3:8])
    rownames(df3_counts) <- rownames(df4)
    df3_counts <- t(scale(t(df3_counts)))
    
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick"))
    
    Heatmap(df3_counts, name = "Row Z-score", col = col_fun, rect_gp = gpar(col = "white", lwd = 2),
            column_split = c(rep("ImP 100 nM", 3), rep("Control",3)),  column_order = colnames(df3_counts)[1:6],
            clustering_method_rows = "ward.D2",
            width = unit(6, "cm"), row_names_side = "left", column_gap = unit(2, "mm"),
            row_dend_side = "right",
            row_names_max_width = max_text_width(
                rownames(df3_counts),
                gp = gpar(fontsize = 10)))
}


## Data
df <- read.csv2("data/220905_deseq2results.csv")
df_old <- read.csv("DEG/deseq2/Ctrl-vs-Treatment/Differential_expression_analysis_table.csv")
names(df)
df$X <- NULL

df$group <- case_when(
    df$log2FoldChange < 0 ~ paste0("lower in ImP"),
    df$log2FoldChange > 0 ~ paste0("higher in ImP")
)
df$sig <- case_when(df$pvalue < 0.05 ~ paste0("sig"), 
                  df$pvalue >0.05 ~ paste0("insig"))
df$Gene.name <- df$symbol

# df$missing <- case_when(df$X12h.100.nM.ImP.3 == 0 | df$X12h.100.nM.ImP.2 == 0 | df$X12h.100.nM.ImP.1 == 0 |
#                             df$X12h.Ctrl.1 == 0 | df$X12h.Ctrl.2 == 0 | df$X12h.Ctrl.1 == 0 ~ TRUE)

# dfsig <- df %>% filter(pvalue < 0.05)
df$padj2 <- p.adjust(df$pvalue, method = "BH", n = length(df$pvalue))
summary(is.na(df$padj))

df2 <- df %>% group_by(group) %>% 
    arrange(padj) %>% 
    dplyr::filter(padj<0.058) %>% 
    ungroup(.) %>% 
    arrange(log2FoldChange) %>% 
    mutate(Gene.name = forcats::fct_inorder(symbol))
nrow(df2)

ggplot(df2, aes(x = Gene.name, y = log2FoldChange, fill = forcats::fct_rev(group))) +
        theme_Publication() +
        geom_bar(stat = "identity", alpha = 0.8) +
        scale_y_continuous(breaks = c(-0.50, -0.25, 0, 0.25, 0.50), limits = c(-0.60, 0.60)) +
        ggsci::scale_fill_lancet() +
        coord_flip() +
        labs(y = 'Log 2 fold change with ImP', x='', 
             title = 'Differential gene expression',
             fill = '') +
        theme(axis.text.x = element_text(size=11)) + 
        theme(axis.text.y = element_text(size=10)) +
        theme(legend.key.size= unit(0.5, "cm")) +
        theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/pdf/220905_diff_exp_sig.pdf", height = 6, width = 6, device = "pdf")
ggsave("results/tiff/220905_diff_exp_sig.tiff", height = 6, width = 6, device = "tiff")
ggsave("results/svg/220905_diff_exp_sig.svg", height = 6, width = 6, device = "svg")
ggsave("results/png/220905_diff_exp_sig.png", height = 6, width = 6, device = "png")

df5 <- df %>% group_by(group) %>% 
    arrange(padj) %>% 
    dplyr::slice(1:10) %>% 
    ungroup(.) %>% 
    arrange(log2FoldChange) %>% 
    mutate(Gene.name = forcats::fct_inorder(symbol))

ggplot(df5, aes(x = Gene.name, y = log2FoldChange, fill = forcats::fct_rev(group))) +
    theme_Publication() +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_y_continuous(breaks = c(-0.50, -0.25, 0, 0.25, 0.50), limits = c(-0.60, 0.60)) +
    ggsci::scale_fill_lancet() +
    coord_flip() +
    labs(y = 'Log 2 fold change with ImP', x='', 
         title = 'Differential gene expression',
         fill = '') +
    theme(axis.text.x = element_text(size=11)) + 
    theme(axis.text.y = element_text(size=10)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/pdf/220905_diff_exp_top20.pdf", height = 6, width = 5, device = "pdf")
ggsave("results/tiff/220905_diff_exp_top20.tiff", height = 6, width = 5, device = "tiff")
ggsave("results/svg/220905_diff_exp_top20.svg", height = 6, width = 5, device = "svg")
ggsave("results/png/220905_diff_exp_top20.png", height = 6, width = 5, device = "png")

pdf("results/pdf/220905_heatmap_significant.pdf", width = 7, height = 7)
heatmap_plot(df, group = "sig", number = 30)
dev.off()
tiff("results/tiff/220905_heatmap_significant.tiff")
heatmap_plot(df, group = "sig", number = 30)
dev.off()
svg("results/svg/220905_heatmap_significant.svg", width = 7, height = 7)
heatmap_plot(df, group = "sig", number = 30)
dev.off()
png("results/png/220905_heatmap_significant.png", width = 600, height = 600)
heatmap_plot(df, group = "sig", number = 30)
dev.off()

pdf("heatmap_top15_uneven.pdf", width = 7, height = 5)
    heatmap_plot(df, group = "unequal", number = 15)
dev.off()

pdf("heatmap_top5_uneven.pdf", width = 5, height = 4)
heatmap_plot(df, group = "unequal", number = 5)
dev.off()


## Boxplot differences
res <- list()
for(a in 1:nrow(df2)){
    gene <- df2[a,]
    genename <- df2$symbol[a]
    genetab <- gene %>% pivot_longer(., cols = c(1:6)) %>% 
        mutate(group2 = case_when(
            str_detect(name, "ImP") ~ paste0("ImP 100 nM"),
            str_detect(name, "Ctrl") ~ paste0("Control")
        ))
    means_genetab <- genetab %>% group_by(group2) %>% 
        summarise(mean_count = mean(value), sd_count = sd(value))
    (pl <- ggplot(data = genetab, aes(x = group2, y = value)) + 
        ggpubr::stat_compare_means(method = "t.test", label = "p.signif", ref.group = "Control",
                                   hide.ns = TRUE) +
        geom_bar(data = means_genetab, aes(x = group2, y = mean_count, 
                                           fill = group2), stat = "identity", color = "black") +
        geom_errorbar(data = means_genetab, 
                      aes(x=group2, y=mean_count, ymin=mean_count, 
                          ymax=mean_count+(sd_count/sqrt(3))), width = 0.5) +
        geom_jitter(color = "black", width = 0.2, height = 0) +
        scale_fill_manual(guide = "none", values = c("white",ggsci::pal_lancet()(2)[2])) +
        labs(title=paste0(genename), y="Normalized counts", x = "") +
        theme_Publication())
    ggsave(str_c("results/pdf/", genename, ".pdf"), width = 3, height = 4, device = "pdf")
    ggsave(str_c("results/tiff/", genename, ".tiff"), width = 3, height = 4, device = "tiff")
    ggsave(str_c("results/svg/", genename, ".svg"), width = 3, height = 4, device = "svg")
    ggsave(str_c("results/png/", genename, ".png"), width = 3, height = 4, device = "png")
    res[[a]] <- pl
}

pdf("results/pdf/boxplots_counts.pdf", width = 12, height = 15)
grid.arrange(grobs=res, ncol=6)
dev.off()

tiff("results/tiff/boxplots_counts.tiff")
grid.arrange(grobs=res, ncol=6)
dev.off()

svg("results/svg/boxplots_counts.svg", width = 12, height = 15)
grid.arrange(grobs=res, ncol=6)
dev.off()

png("results/png/boxplots_counts.png", width = 1500*0.8, height = 1500)
grid.arrange(grobs=res, ncol=6)
dev.off()


df3 <- df %>% group_by(group) %>% 
    arrange(padj) %>% 
    dplyr::filter(str_detect(Gene.name,"PIK3")) %>% 
    ungroup(.) %>% 
    arrange(log2FoldChange) %>% 
    mutate(Gene.name = forcats::fct_inorder(symbol))

res2 <- list()
for(a in 1:nrow(df3)){
    gene <- df3[a,]
    genename <- df3$symbol[a]
    genetab <- gene %>% pivot_longer(., cols = c(1:6)) %>% 
        mutate(group2 = case_when(
            str_detect(name, "ImP") ~ paste0("ImP"),
            str_detect(name, "Ctrl") ~ paste0("Control")
        ))
    means_genetab <- genetab %>% group_by(group2) %>% 
        summarise(mean_count = mean(value), sd_count = sd(value))
    (pl <- ggplot(data = genetab, aes(x = group2, y = value)) + 
            ggpubr::stat_compare_means(method = "t.test", label = "p.signif", ref.group = "Control",
                                       hide.ns = TRUE) +
            geom_bar(data = means_genetab, aes(x = group2, y = mean_count, 
                                               fill = group2), stat = "identity", color = "black") +
            geom_errorbar(data = means_genetab, 
                          aes(x=group2, y=mean_count, ymin=mean_count, 
                              ymax=mean_count+(sd_count/sqrt(3))), width = 0.5) +
            geom_jitter(color = "black", width = 0.2, height = 0) +
            scale_fill_manual(guide = "none", values = c("white",ggsci::pal_lancet()(2)[2])) +
            labs(title=paste0(genename), y="Normalized counts", x = "") +
            theme_Publication())
    ggsave(str_c("results/pdf/", genename, ".pdf"), width = 3, height = 4, device = "pdf")
    ggsave(str_c("results/tiff/", genename, ".tiff"), width = 3, height = 4, device = "tiff")
    ggsave(str_c("results/svg/", genename, ".svg"), width = 3, height = 4, device = "svg")
    ggsave(str_c("results/png/", genename, ".png"), width = 3, height = 4, device = "png")
    res2[[a]] <- pl
}

pdf("results/pdf/pik3_isoforms.pdf", width = 12, height = 12)
grid.arrange(grobs=res2, ncol=5)
dev.off()

tiff("results/tiff/pik3_isoforms.tiff")
grid.arrange(grobs=res2, ncol=5)
dev.off()

svg("results/svg/pik3_isoforms.svg", width = 12, height = 12)
grid.arrange(grobs=res2, ncol=5)
dev.off()

png("results/png/pik3_isoforms.png", width = 700, height = 700)
grid.arrange(grobs=res2, ncol=5)
dev.off()

## Volcano plot

df$sigdir <- case_when(
    df$sig == "sig" & df$group == "higher in ImP" ~ paste0("up"),
    df$sig == "sig" & df$group == "lower in ImP" ~ paste0("down"),
    df$sig == "insig" ~ paste0("no")
)
df$sigdir <- as.factor(df$sigdir)
levels(df$sigdir)

df <- df %>% 
    mutate(
        sigdir = case_when(
            padj < 0.058 & group == "higher in ImP" ~ paste0("up"),
            #padj < 0.01 & group == "higher in ImP" ~ paste0("up"),
            padj < 0.058 & group == "lower in ImP" ~ paste0("down"),
            #padj < 0.01 & group == "lower in ImP" ~ paste0("down"),
            padj >= 0.058 ~ paste0("no")
        ),
        sigdir = as.factor(sigdir),
        sigdir = fct_infreq(sigdir),
        delabel = case_when(
            Gene.name == "PIK3C2A" ~ paste0(Gene.name),
            Gene.name == "ADGRL2" ~ paste0(Gene.name),
            Gene.name == "ARHGAP22" ~ paste0(Gene.name),
            #Gene.name == "ROMO1" ~ paste0(Gene.name),
            Gene.name == "CALCRL" ~ paste0(Gene.name),
            #Gene.name == "STAT1" ~ paste0(Gene.name),
            Gene.name == "ANGPTL4" ~ paste0(Gene.name)
        )
    )

seed(1234)
ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = group, label = delabel)) +
    theme_Publication() +
    geom_hline(aes(yintercept = -log10(0.058)), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 1.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey50", 
                    force = 1, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_lancet()(2)), guide = "none") +
    #geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
    scale_x_continuous(limits = c(-0.6, 0.6)) +
    labs(x = "Log2 fold change with ImP",
         y = "-log10(p-value)") 
ggsave("results/pdf/230320_volcanoplot_adjpval.pdf", width = 5, height = 5, device = "pdf")    
ggsave("results/tiff/230320_volcanoplot_adjpval.tiff", width = 5, height = 5, device = "tiff")  
ggsave("results/svg/230320_volcanoplot_adjpval.svg", width = 5, height = 5, device = "svg")
ggsave("results/png/230320_volcanoplot_adjpval.png", width = 5, height = 5, device = "png") 

seed(1234)
ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = group, label = delabel)) +
    theme_Publication() +
    #geom_hline(aes(yintercept = -log10(0.058)), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 1.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey50", 
                    force = 1, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_lancet()(2)), guide = "none") +
    #geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
    scale_x_continuous(limits = c(-0.6, 0.6)) +
    labs(x = "Log2 fold change with ImP",
         y = "-log10(p-value)") 
ggsave("results/pdf/230320_volcanoplot_noline.pdf", width = 5, height = 5, device = "pdf")    
ggsave("results/tiff/230320_volcanoplot_noline.tiff", width = 5, height = 5, device = "tiff")  
ggsave("results/svg/230320_volcanoplot_noline.svg", width = 5, height = 5, device = "svg")
ggsave("results/png/230320_volcanoplot_noline.png", width = 5, height = 5, device = "png") 

seed(1234)
ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = group, label = delabel)) +
    theme_Publication() +
    #geom_hline(aes(yintercept = -log10(0.058)), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 1.5, min.segment.length = 0,
                    point.padding = 0, color = "black", fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey50", 
                    force = 1, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_lancet()(2)), guide = "none") +
    #geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
    scale_x_continuous(limits = c(-0.6, 0.6)) +
    labs(x = "Log2 fold change with ImP",
         y = "-log10(p-value)")
ggsave("results/pdf/230320_volcanoplot_noline.pdf", width = 5, height = 5, device = "pdf")    
ggsave("results/tiff/230320_volcanoplot_noline.tiff", width = 5, height = 5, device = "tiff")  
ggsave("results/svg/230320_volcanoplot_noline.svg", width = 5, height = 5, device = "svg")
ggsave("results/png/230320_volcanoplot_noline.png", width = 5, height = 5, device = "png") 
