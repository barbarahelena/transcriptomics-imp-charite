# Pathway analysis

# libraries
library(gage)
library(DESeq2)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(pathview)

# Open file with DESeq2 results
deg <- read.csv2("data/220905_deseq2results.csv")
rownames(deg) <- deg$ENSEMBL

# GO
# set up kegg database
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP] # biological process
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF] # molecular function
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC] # cellular component

# annotate the deseq2 results with additional gene identifiers
deg$symbol <- mapIds(EnsDb.Hsapiens.v86, keys=deseq2Results$ENSEMBL,
                     column="SYMBOL", keytype="GENEID", multiVals="first")
deg$entrez <- mapIds(EnsDb.Hsapiens.v86, keys=deseq2Results$ENSEMBL,
                     column="ENTREZID", keytype="GENEID", multiVals="first")
deg$name <- mapIds(EnsDb.Hsapiens.v86, keys=deseq2Results$ENSEMBL,
                   column="GENENAME", keytype="GENEID", multiVals="first")

# grab the log fold changes for everything
foldchanges <- deg$log2FoldChange
names(foldchanges) <- deg$entrez

# Run enrichment analysis on all log fc
fc.kegg.sigmet.p <- gage(foldchanges, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(foldchanges, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(foldchanges, gsets = go.bp.gs)
fc.go.mf.p <- gage(foldchanges, gsets = go.mf.gs)
fc.go.cc.p <- gage(foldchanges, gsets = go.cc.gs)

# covert the kegg results to data frames
fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

# convert the go results to data frames
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)

# Check if any significant pathways
any(fc.go.bp.p.up$q.val < 0.05)
any(fc.go.bp.p.down$q.val < 0.05)
any(fc.go.cc.p.down$q.val < 0.05)
any(fc.go.cc.p.up$q.val < 0.05)
any(fc.go.mf.p.down$q.val < 0.05)
any(fc.go.mf.p.up$q.val < 0.05)

kegg_a <- fc.kegg.sigmet.p.up %>% 
    arrange(p.val) %>% 
    dplyr::slice(1:10) %>% 
    arrange(-p.val) %>% 
    mutate(group = "up", pathway = rownames(.),
           pathway = fct_inorder(pathway))

kegg_b <- fc.kegg.sigmet.p.down %>% 
    arrange(p.val) %>% 
    dplyr::slice(1:10) %>% 
    arrange(-p.val) %>% 
    mutate(group = "down", pathway = rownames(.),
           pathway = fct_inorder(pathway))

# fc.kegg.sigmet.ab <- rbind(a,b)
# fc.kegg.sigmet.ab$pathway <- rownames(fc.kegg.sigmet.ab)

kegg_c <- fc.kegg.dise.p.down %>% 
    arrange(p.val) %>% 
    dplyr::slice(1:10) %>% 
    arrange(-p.val) %>% 
    mutate(group = "down", pathway = rownames(.),
           pathway = fct_inorder(pathway))

kegg_d <- fc.kegg.dise.p.up %>% 
    arrange(p.val) %>% 
    dplyr::slice(1:10) %>% 
    arrange(-p.val) %>% 
    mutate(group = "up", pathway = rownames(.),
           pathway = fct_inorder(pathway))


down <- ggsci::pal_lancet()(2)[1]
up <- ggsci::pal_lancet()(2)[2]

plc <- ggplot(kegg_c, aes(x = pathway, y = p.val)) +
    theme_Publication() +
    geom_bar(stat = "identity", fill = down) +
    coord_flip() +
    labs(y = 'p-value', x='', 
         title = 'KEGG disease pathways: downregulated',
         fill = '') +
    theme(axis.text.x = element_text(size=10)) + 
    theme(axis.text.y = element_text(size=8)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')

pld <- ggplot(kegg_d, aes(x = pathway, y = p.val)) +
    theme_Publication() +
    geom_bar(stat = "identity", fill = up) +
    coord_flip() +
    labs(y = 'p-value', x='', 
         title = 'KEGG disease pathways: upregulated',
         fill = '') +
    theme(axis.text.x = element_text(size=10)) + 
    theme(axis.text.y = element_text(size=8)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')

ggpubr::ggarrange(plc, pld)
ggsave("kegg_diseasepathways.pdf", width = 13, height = 5)

pla <- ggplot(kegg_b, aes(x = pathway, y = p.val)) +
    theme_Publication() +
    geom_bar(stat = "identity", fill = down) +
    coord_flip() +
    labs(y = 'p-value', x='', 
         title = 'KEGG sig met: downregulated',
         fill = '') +
    theme(axis.text.x = element_text(size=10)) + 
    theme(axis.text.y = element_text(size=8)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')

plb <- ggplot(kegg_a, aes(x = pathway, y = p.val)) +
    theme_Publication() +
    geom_bar(stat = "identity", fill = up) +
    coord_flip() +
    labs(y = 'p-value', x='', 
         title = 'KEGG sig met: upregulated',
         fill = '') +
    theme(axis.text.x = element_text(size=10)) + 
    theme(axis.text.y = element_text(size=8)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')

ggpubr::ggarrange(pla, plb)
ggsave("kegg_sigmet.pdf", width = 13, height = 5)

# Overlay the expression data onto this pathway
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa04613")
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa04930")
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa05034")
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa05322")

# Overlay the expression data onto this pathway
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa04613", kegg.native=FALSE)
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa04930", kegg.native=FALSE)
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa05034", kegg.native=FALSE)
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa05322", kegg.native=FALSE)
