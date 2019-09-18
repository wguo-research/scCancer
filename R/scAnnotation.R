filterCell <- function(cell.manifest, filter.thres){
    rownames(filter.thres) <- filter.thres$Index
    ix.cr <- which(cell.manifest$droplet.type == "cell")
    ix.thres <- getCellix(cell.manifest, filter.thres,
                          arg = c("nUMI", "nGene", "mito.percent", "ribo.percent", "diss.percent"))
    cell.manifest <- cell.manifest[intersect(ix.cr, ix.thres), ]
    rownames(cell.manifest) <- cell.manifest$barcodes

    return(cell.manifest)
}




filterGene <- function(gene.manifest,
                       anno.filter = c("mitochondrial", "ribosome", "dissociation"),
                       nCell.min = 3,
                       bgPercent.max = 0.1,
                       savePath = NULL){
    if(!is.null(anno.filter)){
        gene.manifest.left <- subset(gene.manifest, !(Annotation %in% anno.filter))
    }else{
        gene.manifest.left <- gene.manifest
    }
    gene.manifest.left <- subset(gene.manifest.left, nCell >= nCell.min & bg.percent <= bgPercent.max)
    gene.manifest.filter <- subset(gene.manifest, !(EnsemblID %in% gene.manifest.left$EnsemblID))

    if(!is.null(savePath)){
        write.table(gene.manifest.filter, file = file.path(savePath, "gene.manifest.filter.txt"),
                    quote = F, sep = "\t", row.names = F)
    }
    return(gene.manifest.left)
}




#' getFilterData
#'
#' According to the QC results of scStatistics, filter cells and genes.
#'
#' @inheritParams runScAnnotation
#'
#' @return A list containing expr.data, cell.manifest, gene.manifest, filter.thres.
#' @export
#'
#' @examples
getFilterData <- function(dataPath, statPath, savePath = NULL,
                          bool.filter.cell = T, bool.filter.gene = T,
                          anno.filter = c("mitochondrial", "ribosome", "dissociation"),
                          nCell.min = 3, bgPercent.max = 0.1,
                          hg.mm.mix = F){

    expr.data <- Read10Xdata(data.dir = get10Xpath(dataPath, raw.data = T))
    gene.manifest <- read.table(file.path(statPath, 'geneManifest.txt'), header = T, sep = "\t")
    cell.manifest <- read.table(file.path(statPath, 'cellManifest-all.txt'), header = T)
    filter.thres <- read.table(file.path(statPath, 'cell.QC.thres.txt'), header = T)

    if(hg.mm.mix){
        rownames(expr.data) <- substr(rownames(expr.data), 6, 60)
    }

    if(bool.filter.cell){
        cell.manifest <- filterCell(cell.manifest, filter.thres = filter.thres)
    }
    if(bool.filter.gene){
        gene.manifest <- filterGene(gene.manifest,
                                    anno.filter = anno.filter,
                                    nCell.min = nCell.min,
                                    bgPercent.max = bgPercent.max,
                                    savePath = savePath)
    }

    expr.data <- expr.data[gene.manifest$Symbol, rownames(cell.manifest)]

    exprList <- list(expr.data = expr.data,
                     cell.manifest = cell.manifest,
                     gene.manifest = gene.manifest,
                     filter.thres = filter.thres)
    return(exprList)
}



prepareSeurat <- function(expr.data, cell.manifest, label = "sc",
                          vars.to.regress = c("nUMI", "mito.percent", "ribo.percent")){
    expr = CreateSeuratObject(counts = expr.data,
                              min.cells = 0,
                              min.features = 0,
                              project = label)

    metaVar <- c("mito.percent", "ribo.percent", "diss.percent")
    for(mv in metaVar){
        tmp <- cell.manifest[[mv]]
        names(tmp) <- rownames(cell.manifest)
        expr[[mv]] <- tmp
    }

    expr <- NormalizeData(object = expr,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)

    message("[", Sys.time(), "] -----: scale data")
    expr <- ScaleData(object = expr,
                      vars.to.regress = vars.to.regress)

    return(expr)
}




#' runSeurat
#'
#' Perform usual Seurat step and cell type prediction.
#'
#' @param exprList A list containing expr.data, cell.manifest, gene.manifest and filter.thres.
#' It should be the return value of function 'getFilterData' generally.
#' @inheritParams runScAnnotation
#'
#' @return A list containing a Seurat object, differential expressed genes and annotation information for cells.
#' @export
#'
#' @examples
runSeurat <- function(exprList, savePath, sampleName = "sc",
                      vars.to.regress = c("nUMI", "mito.percent", "ribo.percent"),
                      pc.use = 30, resolution = 0.8,
                      clusterStashName = "default",
                      bool.runDiffExpr = T){

    if(!dir.exists(file.path(savePath))){
        dir.create(file.path(savePath), recursive = T)
    }

    expr.data <- exprList$expr.data
    cell.manifest <- exprList$cell.manifest

    message("[", Sys.time(), "] -----: Seurat object preparation")
    expr <- prepareSeurat(expr.data, cell.manifest, label = sampleName,
                          vars.to.regress = vars.to.regress)

    message("[", Sys.time(), "] -----: highly variable genes")
    expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000)

    message("[", Sys.time(), "] -----: PCA")
    expr <- RunPCA(expr, features = VariableFeatures(expr), verbose = F)

    message("[", Sys.time(), "] -----: Clustering")
    expr <- FindNeighbors(expr, dims = 1:pc.use, verbose = F)
    expr <- FindClusters(expr, resolution = resolution, verbose = F)
    expr[[clusterStashName]] <- as.numeric(Idents(object = expr))

    message("[", Sys.time(), "] -----: tSNE")
    expr <- RunTSNE(object = expr, dims = 1:pc.use)

    if(bool.runDiffExpr){
        message("[", Sys.time(), "] -----: differential expression analysis")
        if(!dir.exists(file.path(savePath, "diff.expr.genes"))){
            dir.create(file.path(savePath, "diff.expr.genes"), recursive = T)
        }
        diff.expr.genes <- FindAllMarkers(expr, only.pos = TRUE,
                                          min.pct = 0.25,
                                          logfc.threshold = 0.25,
                                          verbose = F)
        # write.table(diff.expr.genes[, c("cluster", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")],
        #             file = file.path(savePath, "diff.expr.genes.txt"), quote = F, sep = "\t", row.names = F)

        diff.expr.genes <- diff.expr.genes[, c("cluster", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")]
        nCluster <- length(unique(diff.expr.genes$cluster))
        for(ci in 1:nCluster){
            # write.table(subset(diff.expr.genes, cluster == ci),
            #             file = file.path(savePath, "diff.expr.genes", paste0("cluster", ci ,".txt")),
            #             quote = F, sep = "\t", row.names = F)
            write.csv(subset(diff.expr.genes, cluster == ci),
                      file = file.path(savePath, "diff.expr.genes", paste0("cluster", ci ,".csv")),
                      quote = F, row.names = F)
        }
        # write.xlsx(subset(diff.expr.genes, cluster == 0),
        #            file = file.path(savePath, "diff.expr.genes.xlsx"),
        #            sheetName = "cluster0", row.names = F)
        # nCluster <- length(unique(diff.expr.genes$cluster))
        # for(ci in 1:(nCluster - 1)){
        #     write.xlsx(subset(diff.expr.genes, cluster == ci),
        #                file = file.path(savePath, "diff.expr.genes.xlsx"),
        #                sheetName = paste0("cluster", ci), append = T, row.names = F)
        # }
    }else{
        diff.expr.genes <- NULL
    }

    cell.annotation <- data.frame(barcodes = colnames(x = expr), stringsAsFactors = F)
    cell.annotation <- cbind(cell.annotation, expr@reductions$tsne@cell.embeddings)
    # cell.annotation$Cluster <- factor(expr@meta.data[[clusterStashName]],
    #                                   levels = 0:(length(unique(expr@meta.data[[clusterStashName]]))-1))
    cell.annotation$Cluster <- factor(expr@meta.data[[clusterStashName]])


    message("[", Sys.time(), "] -----: cell cycle score estimation")
    cellCycle.genes <- read.table(system.file("txt", "cellCycle-genes.txt", package = "scCancer"),
                                  header = F, stringsAsFactors = F)$V1
    expr <- AddModuleScore(expr, features = list(cellCycle.genes), name = "cellCycle")
    cell.annotation$CellCycle.score <- expr[["cellCycle1"]]$cellCycle1


    return(list(expr = expr,
                diff.expr.genes = diff.expr.genes,
                cell.annotation = cell.annotation))
}




singleGenePlot <- function(expr.data, gene,
                           coor.df, coor.names = c("tSNE_1", "tSNE_2"),
                           color = "blue", font.size = 8,
                           legend = F, axis = F){
    minx <- min(coor.df[, coor.names[1]])
    miny <- min(coor.df[, coor.names[2]])

    ju.zero <- F
    if(!(gene %in% rownames(expr.data))){
        ju.zero <- T
    }else if(sum(expr.data[gene, ] > 0) == 0){
        ju.zero <- T
    }

    if(ju.zero){
        warning(gene, " isn't detected in the expression matrix.\n")
        new.label <- paste0(gene, " (0/", dim(expr.data)[2], ")")
        p <- ggplot(coor.df, aes(x = coor.df[, coor.names[1]], y = coor.df[, coor.names[2]])) +
            geom_point(shape = 21, size = 1, stroke = 0.25, color = "grey", fill = "white") +
            annotate("text", x = minx, y = miny, label = new.label, hjust = 0, vjust = 0,
                     size = font.size, fontface = 'italic') +
            theme_classic()
    }else{
        coor.df$cur.value <- expr.data[gene, ]
        new.label <- paste0(gene, " (", sum(expr.data[gene, ] > 0), "/", dim(expr.data)[2], ")")
        p <- ggplot(coor.df, aes(x = coor.df[, coor.names[1]], y = coor.df[, coor.names[2]],
                                 fill = cur.value)) +
            geom_point(shape = 21, size = 1, stroke = 0.25, color = "grey") +
            scale_fill_gradientn(colors = c("white", color),
                                 guide = guide_colorbar(title = gene)) +
            annotate("text", x = minx, y = miny, label = new.label, hjust = 0, vjust = 0,
                     size = font.size, fontface = 'italic') +
            theme_classic()
    }

    if(!legend){
        p <- p + theme(legend.position = 'none')
    }
    if(!axis){
        p <- p + theme(
            panel.border = element_rect(fill = NA),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank()
        )
    }
    return(p)
}



markerPlot <- function(expr.data, coor.df, coor.names = c("tSNE_1", "tSNE_2"),
                       features = NULL, add = T,
                       species = "human",
                       font.size = 8, color = "blue"){
    feature.def <- getDefaultMarkers(species = species)
    feature.def <- unlist(feature.def)

    if(is.null(features)){
        features <- feature.def
    }else if(add){
        features <- c(feature.def, unlist(features))
    }

    ps <- list()
    for(g in features){
        ps[[g]] <- singleGenePlot(expr.data = expr.data, gene = g,
                                  coor.df = coor.df, coor.names = coor.names,
                                  font.size = font.size, color = "blue")
    }

    return(ps)
}



pointDRPlot <- function(cell.annotation, value,
                        colors, discrete = T,
                        limit.quantile = 0,
                        legend.position = "right",
                        legend.title = NULL){
    if(is.null(legend.title)){
        legend.title <- value
    }

    ratio <- with(cell.annotation, diff(range(tSNE_1))/diff(range(tSNE_2)))

    fill.value <- cell.annotation[[value]]
    if(!discrete){
        low.thres <- quantile(cell.annotation[[value]], limit.quantile)
        up.thres <- quantile(cell.annotation[[value]], 1 - limit.quantile)
        fill.value <- ifelse(fill.value < low.thres, low.thres,
                             ifelse(fill.value > up.thres, up.thres, fill.value))
    }

    p <- ggplot(cell.annotation, aes(x = tSNE_1, y = tSNE_2, fill = fill.value)) +
        geom_point(shape = 21, size = 1, stroke = 0.25, color = "lightgrey") +
        coord_fixed(ratio = ratio) +
        ggplot_config(base.size = 6) +
        labs(x = "t-SNE 1", y = "t-SNE 2") +
        guides(fill = guide_legend(override.aes = list(size = 3),
                                   keywidth = 0.1,
                                   keyheight = 0.15,
                                   default.unit = "inch",
                                   title = legend.title)) +
        theme(legend.position = legend.position)

    if(discrete){
        p <- p + scale_fill_manual(values = colors)
    }else{
        p <- p + scale_fill_gradientn(colors = colors)
    }

    return(p)
}



clusterBarPlot <- function(cell.annotation, cell.colors, sel.col = "Cell.Type", legend.title = NULL){
    if(is.null(legend.title)){
        legend.title <- sel.col
    }

    bar.df <- melt(table(cell.annotation[c(sel.col, "Cluster")]))
    # bar.df$Cluster = factor(bar.df$Cluster,
    #                         levels = 0:(length(unique(bar.df$Cluster))-1))
    bar.df$Cluster = factor(bar.df$Cluster)

    p <- ggplot(bar.df, aes(x = Cluster, y = value, fill = bar.df[[sel.col]])) +
        geom_bar(stat = "identity") +
        ggplot_config(base.size = 6) +
        scale_fill_manual(values = cell.colors) +
        labs(y = "Number of cells") +
        guides(fill = guide_legend(override.aes = list(size=1), title = legend.title))
    return(p)
}



preDEheatmap <- function(expr, cell.annotation, genes = NULL, cells = NULL,
                         sel.col.anno = c("Cluster", "Cell.Type"),
                         slot = "scale.data",
                         min.value = -2.5, max.value = 2.5){

    expr.data <- GetAssayData(object = expr, slot = slot)
    if(!is.null(genes)){
        genes <- intersect(genes, rownames(expr.data))
        expr.data <- expr.data[genes, ]
    }
    if(!is.null(cells)){
        cells <- intersect(cells, colnames(expr.data))
        expr.data <- expr.data[, cells]
    }

    rownames(cell.annotation) <- cell.annotation$barcodes
    cell.annotation <- cell.annotation[colnames(expr.data), ]
    cell.cluster <- cell.annotation[order(cell.annotation$Cluster), sel.col.anno, drop = FALSE]
    expr.data <- expr.data[, rownames(cell.cluster)]

    ## limitData
    expr.data <- limitData(expr.data, min = min.value, max = max.value)

    ## gaps_col
    num.cluster <- table(cell.cluster$Cluster)
    num.cluster <- num.cluster[as.character(1 : length(num.cluster))]
    gaps_col <- cumsum(num.cluster)

    return(list(expr.data = expr.data,
                cell.cluster = cell.cluster,
                gaps_col = gaps_col))
}





#' plotSeurat
#'
#' Construct and save plots of Seurat analysis.
#'
#' @param expr A Seurat object return by function 'runSeurat'.
#' @param cell.annotation A data.frame of cells' annotation return by function 'runSeurat'.
#' @param diff.expr.genes A data.frame of differential expressed genes return by function 'runSeurat'.
#' @inheritParams runScAnnotation
#'
#' @return A list of all plots generated by Seurat analyses.
#' @export
#'
#' @examples
plotSeurat <- function(expr,
                       cell.annotation = cell.annotation,
                       show.features = NULL, bool.add.features = T,
                       coor.names = c("tSNE_1", "tSNE_2"),
                       bool.runDiffExpr = T,
                       diff.expr.genes = NULL, n.markers = 5,
                       species = "human",
                       savePath){

    message("[", Sys.time(), "] -----: Seurat plotting and saving")

    if(!dir.exists(file.path(savePath, "figures/singleMarkerPlot/"))){
        dir.create(file.path(savePath, "figures/singleMarkerPlot/"), recursive = T)
    }

    p.results <- list()

    # message(sprintf('------p.hvg------'))
    top10 <- head(VariableFeatures(expr), 10)
    suppressWarnings(
        plot1 <- VariableFeaturePlot(expr, cols = c("grey", "#ec7d89"))
    )
    suppressWarnings(
        p.results[["p.hvg"]] <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + NoLegend()
    )


    # message(sprintf('------p.markers.all------'))
    p.results[["ps.markers"]] <- markerPlot(expr.data = GetAssayData(expr),
                                            coor.df = cell.annotation,
                                            coor.names = coor.names,
                                            features = show.features, add = bool.add.features,
                                            species = species,
                                            font.size = 4, color = "blue")
    if(length(p.results[["ps.markers"]]) > 0){
        p.results[["p.markers.all"]] <- plot_grid(plotlist = p.results[["ps.markers"]], ncol = 4)
    }else{
        p.results[["p.markers.all"]] <- NULL
    }


    # message(sprintf('------p.cluster------'))
    def.colors <- getDefaultColors()
    clusters <- unique(cell.annotation$Cluster)
    clusters <- sort(clusters)
    cluster.colors <- c()
    for(i in 1:length(clusters)){
        cluster.colors[as.character(clusters[i])] = def.colors[i]
    }

    p.results[["p.cluster"]] <- pointDRPlot(cell.annotation, value = "Cluster",
                                            colors = cluster.colors,
                                            legend.title = "Cluster")


    if(bool.runDiffExpr && !(is.null(diff.expr.genes))){
        # message(sprintf('------p.DE.heatmap------'))
        top.genes <- diff.expr.genes %>% group_by(cluster) %>% top_n(n = n.markers, wt = avg_logFC)
        de.pre <- preDEheatmap(expr = expr,
                               cell.annotation = cell.annotation,
                               genes = top.genes$gene,
                               sel.col.anno = c("Cluster"),
                               slot = "scale.data",
                               min.value = -2.5, max.value = 2.5)


        ann_colors = list(Cluster = cluster.colors)
        p.results[["p.de.heatmap"]] <-
            pheatmap(de.pre$expr.data,
                     color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100),
                     annotation_col = de.pre$cell.cluster,
                     annotation_colors = ann_colors,
                     fontsize = 7,
                     gaps_col = de.pre$gaps_col,
                     cluster_rows = F, cluster_cols = F,
                     show_colnames = F, legend = F,
                     silent = T)
    }

    # message(sprintf('------p.cellCycle------'))
    p.results[["p.cellCycle"]] <- pointDRPlot(cell.annotation, value = "CellCycle.score",
                                              colors = c("white", "#009b45"),
                                              discrete = F,
                                              legend.position = "bottom",
                                              legend.title = "Cell cycle score")

    # message(sprintf('------save images------'))
    ggsave(filename = file.path(savePath, "figures/hvg.png"), p.results[["p.hvg"]],
           width = 8, height = 4, dpi = 800)

    for(gene in names(p.results[["ps.markers"]])){
        ggsave(filename = paste0(savePath, "/figures/singleMarkerPlot/", gene, ".png"),
               p.results[["ps.markers"]][[gene]], width = 5, height = 5, dpi = 800)
    }

    if(length(p.results[["ps.markers"]]) > 0){
        markersPlot.height <- 2 * ceiling(length(p.results[["ps.markers"]]) / 4)
        ggsave(filename = file.path(savePath, "figures/markers-all.png"),
               p.results[["p.markers.all"]], width = 8, height = markersPlot.height, dpi = 800)
    }

    ggsave(filename = file.path(savePath, "figures/cluster-point.png"),
           p.results[["p.cluster"]], width = 6, height = 5, dpi = 800)


    if(bool.runDiffExpr && !(is.null(diff.expr.genes))){
        DEplot.height <- 0.5 + 0.1 * n.markers * length(unique(cell.annotation$Cluster))
        ggsave(filename = file.path(savePath, "figures/DE-heatmap.png"),
               p.results[["p.de.heatmap"]], width = 8, height = DEplot.height, dpi = 800)
    }

    ggsave(filename = file.path(savePath, "figures/cellCycle-point.png"),
           p.results[["p.cellCycle"]], width = 6, height = 5, dpi = 800)

    # saveRDS(cell.annotation, file = file.path(savePath, "cell.annotation.RDS"))

    return(p.results)
}




#' predCellType
#'
#' @param X.test A cells expression matrix (row for genes, column for cells).
#' @param ct.templates A list of gene weight vectors for each cell type.
#'
#' @return A list of predicted cell types and the relative correlations.
#' @export
#'
predCellType <- function(X.test, ct.templates = NULL){

    if(is.null(ct.templates)){
        ct.templates <- readRDS(system.file("rds", "cellTypeTemplates.RDS", package = "scCancer"))
    }

    cor.df <- list()
    for(i in 1:length(ct.templates)){
        type <- names(ct.templates)[i]
        common.gene <- intersect(rownames(X.test), names(ct.templates[[type]]))
        s.cor <- cor(as.matrix(X.test[common.gene, ]),
                     ct.templates[[type]][common.gene],
                     method = "spearman")
        s.cor[is.na(s.cor)] <- 0
        cor.df[[type]] <- s.cor
    }
    cor.df <- as.data.frame(cor.df)

    type.pred <- colnames(cor.df)[unlist(apply(cor.df, 1, which.max))]
    # type.pred[rowSums(cor.df > 0.25) == 0] <- "Unknown"

    thres <- apply(cor.df, 2, quantile, probs = 0.1)
    type.pred[rowSums(t(t(cor.df) > thres)) == 0] <- "Unknown"
    type.pred[rowSums(cor.df > 0.1) == 0] <- "Unknown"

    names(type.pred) <- colnames(X.test)

    colnames(cor.df) <- paste0(colnames(cor.df), ".corr")
    rownames(cor.df) <- colnames(X.test)

    return(list(type.pred = type.pred,
                cor.df = cor.df))
}





#' runCellClassify
#'
#' Use a one-class logistic regression (OCLR) model to predict cancer microenvironment cell types.
#'
#' @param expr A Seurat object return by function 'runSeurat'.
#' @param cell.annotation A data.frame of cells' annotation return by function 'runSeurat'.
#' @inheritParams runScAnnotation
#'
#' @return A list of updated Seurat object, cell.annotation, and the plots for cell type annotation.
#' @export
#'
#' @examples
runCellClassify <- function(expr, cell.annotation, savePath, ct.templates = NULL){
    t.results <- predCellType(X.test = GetAssayData(expr), ct.templates = ct.templates)

    expr[["Cell.Type"]] <- t.results$type.pred

    cell.annotation$Cell.Type <- t.results$type.pred
    cell.annotation <- cbind(cell.annotation, t.results$cor.df)

    cell.colors <- c(
        "T.cells.CD4" = "#07a2a4",
        "T.cells.CD8" = "#9a7fd1",
        "B.cells" = "#588dd5",
        "NK.cells" = "#f5994e",
        "Myeloid.cells" = "#c05050",
        "Endothelial" = "#59678c",
        "Fibroblast" = "#c9ab00",
        "Epithelial" = "#7eb00a",
        "Unknown" = "grey")
    cti = 1
    for(ct in setdiff(cell.annotation$Cell.Type, names(cell.colors))){
        cell.colors[ct] <- getDefaultColors()[cti]
        cti = cti + 1
    }

    # message(sprintf('------p.type------'))
    p.type <- pointDRPlot(cell.annotation, value = "Cell.Type",
                          colors = cell.colors,
                          legend.title = "Cell type")

    # message(sprintf('------p.bar------'))
    p.bar <- clusterBarPlot(cell.annotation = cell.annotation,
                            cell.colors = cell.colors,
                            sel.col = "Cell.Type",
                            legend.title = "Cell type")

    ggsave(filename = file.path(savePath, "figures/cellType-point.png"),
           p.type, width = 7, height = 4, dpi = 800)
    ggsave(filename = file.path(savePath, "figures/cellType-bar.png"),
           p.bar, width = 6, height = 4, dpi = 800)

    return(list(expr = expr,
                cell.annotation = cell.annotation,
                p.results = list(p.type = p.type,
                                 p.bar = p.bar)))
}



#' runGeneSets
#'
#' Calculate gene set signature scores for cells.
#'
#' @param expr A Seurat object.
#' @param method The method to be used in calculate gene set scores. Currently, only "average" and "GSVA" are allowed.
#' @inheritParams runScAnnotation
#'
#' @return A data.frame of calculated gene set signature scores.
#' @export
#'
#' @importFrom GSVA gsva
#'
#' @examples
runGeneSets <- function(expr, geneSets = NULL, method = "average"){
    message("[", Sys.time(), "] -----: calculate gene sets scores")
    if(is.null(geneSets)){
        geneSets <- readLines(system.file("txt", "hallmark-pathways.txt", package = "scCancer"))
        geneSets <- strsplit(geneSets, "\t")
        geneSets <- as.data.frame(geneSets, stringsAsFactors = F)
        colnames(geneSets) <- geneSets[1, ]
        geneSets <- as.list(geneSets[2, ])
        geneSets <- sapply(geneSets, function(x) strsplit(x, ", "))
    }else{
        if(class(geneSets) != "list"){
            warning("The 'geneSets' should be a list of several gene sets.\n")
            return(NULL)
        }
    }
    if(method == "average"){
        expr <- AddModuleScore(expr, features = geneSets, name = "geneSets")
        t.scores <- expr[[paste0("geneSets", 1:length(geneSets))]]
        t.scores <- scale(t.scores)
    }else if(method == "GSVA"){
        tmp.data <- as.matrix(GetAssayData(object = expr, slot = "scale.data"))
        tmp.data <- tmp.data[VariableFeatures(expr), ]
        t.scores <- gsva(tmp.data, geneSets)
        t.scores <- t(t.scores)
    }else{
        warning("The 'method' ", method, " is not allowed.\n")
        return(NULL)
    }
    colnames(t.scores) <- paste0("GS__", names(geneSets))
    return(t.scores)
}




#' plotGeneSet
#'
#' @param cell.annotation A data.frame of cells' annotation containing gene set signature scores.
#' @param prefix A prefix string of column names for gene sets.
#' @param bool.limit A logical value indicating whether to set upper and lower limit when plot heatmap.
#' @inheritParams runScAnnotation
#'
#' @return A heatmap for gene set signature scores.
#' @export
#'
#' @examples
plotGeneSet <- function(cell.annotation, prefix = "GS__", bool.limit = T, savePath = NULL){
    gs.name <- colnames(cell.annotation)
    gs.name <- grep(paste0("^", prefix), gs.name, value = TRUE)
    data <- cell.annotation[, gs.name]
    colnames(data) <- substr(gs.name, 1+nchar(prefix), 200)

    if(bool.limit){
        low.bound <- quantile(as.matrix(data), 0.01)
        up.bound <- quantile(as.matrix(data), 0.99)
        data <- limitData(data, low.bound, up.bound)
    }

    tmp.results <- getClusterInfo(cell.annotation)
    cluster.info <- tmp.results$cluster.info
    cluster.colors <- tmp.results$cluster.colors
    cluster.pos <- tmp.results$cluster.pos

    p <- pheatmap(t(data)[, rownames(cluster.info)],
                  show_colnames = F,
                  cluster_cols = F,
                  fontsize = 7,
                  annotation_col = cluster.info,
                  annotation_colors = cluster.colors,
                  gaps_col = cluster.pos,
                  legend = F,
                  silent = T)

    if(!is.null(savePath)){
        ggsave(filename = file.path(savePath, "figures/geneSet-heatmap.png"),
               p, width = 10, height = 6, dpi = 1000)
    }

    return(p)
}




#' runExprProgram
#'
#' Perform non-negative matrix factorization (NMF) to identify expression programs.
#'
#' @param expr A Seurat object.
#' @param rank An integer of decomposition rank used in NMF.
#' @inheritParams runScAnnotation
#'
#' @return A list of decomposed matrixes (W and H), and the relative genes of each programs.
#' @export
#'
#' @importFrom methods as
#'
#' @examples
runExprProgram <- function(expr, rank = 50, savePath = NULL){
    message("[", Sys.time(), "] -----: identify expression programs")

    data <- as(object = expr[["RNA"]]@data, Class = "dgTMatrix")
    ave.data <-  Matrix::rowSums(data) / Matrix::rowSums(data > 0)

    data@x <- data@x - ave.data[data@i + 1]
    data@x[data@x < 0] <- 0

    nmf.results <- nnmf(as.matrix(data), k = rank)

    W <- nmf.results$W
    colnames(W) <- paste0("p", 1:dim(W)[2])
    H <- nmf.results$H
    rownames(H) <- paste0("p", 1:dim(H)[1])

    all.genes <- rownames(W)
    sel.W <- (W > (1 - 50/dim(W)[1]))
    for(pi in 1:dim(W)[2]){
        tmp <- data.frame(program = colnames(W)[pi], gene = all.genes[sel.W[, pi]], value = W[sel.W[, pi], pi])
        tmp <- tmp[order(tmp$value, decreasing = T), ]
        if(pi == 1){
            program.gene.value <- tmp
        }else{
            program.gene.value <- rbind(program.gene.value, tmp)
        }
    }

    # programs.geneList <- apply(sel.W, 2, FUN = function(x){ return(all.genes[x])})

    if(!is.null(savePath)){
        if(!dir.exists(file.path(savePath, "expr.programs/"))){
            dir.create(file.path(savePath, "expr.programs/"), recursive = T)
        }
        write.table(W, file = file.path(savePath, "expr.programs/W-gene-program.txt"),
                    quote = F, sep = "\t")
        write.table(H, file = file.path(savePath, "expr.programs/H-program-cell.txt"),
                    quote = F, sep = "\t")
        write.table(program.gene.value, file = file.path(savePath, "expr.programs/program.gene.value.txt"),
                    quote = F, sep = "\t", row.names = F)

        # cat("", file = file.path(savePath, "expr.programs/programs.geneList.txt"))
        # for(p in names(programs.geneList)){
        #     cat(p, "\t", str_c(programs.geneList[[p]], collapse = ", "), "\n", append = T,
        #         file = file.path(savePath, "expr.programs/programs.geneList.txt"))
        # }
    }
    return(list(W = W, H = H,
                program.gene.value = program.gene.value))
}





#' plotExprProgram
#'
#' @param H The decomposed right matrix H.
#' @param cell.annotation A data.frame of cells' annotation containing cluster information.
#' @param bool.limit A logical value indicating whether to set upper and lower limit when plot heatmap.
#' @inheritParams runScAnnotation
#'
#' @return A heatmap for cells' expression programs.
#' @export
#' @importFrom NNLM nnmf
#'
#' @examples
plotExprProgram <- function(H, cell.annotation, bool.limit = F, savePath = NULL){
    if(bool.limit){
        up.bound <- quantile(as.matrix(H), 0.995)
        H <- limitData(H, low.bound, up.bound)
    }

    tmp.results <- getClusterInfo(cell.annotation)
    cluster.info <- tmp.results$cluster.info
    cluster.colors <- tmp.results$cluster.colors
    cluster.pos <- tmp.results$cluster.pos

    p <- pheatmap(H[, rownames(cluster.info)],
                  show_colnames = F,
                  cluster_cols = F, fontsize = 7,
                  annotation_col = cluster.info,
                  annotation_colors = cluster.colors,
                  gaps_col = cluster.pos,
                  color = colorRampPalette(colors = c("#f9fcfb","#009b45"))(100),
                  silent = T)

    if(!is.null(savePath)){
        ggsave(filename = file.path(savePath, "figures/exprProgram-heatmap.png"),
               p, width = 10, height = 5.5, dpi = 1000)
    }
    return(p)
}



#' runScAnnotation
#'
#' Use Seurat package to perform data normalization,
#' highly variable genes identification, dimensionality reduction, clustering and differential expression analysis.
#' Predict cell type. Infer CNV.
#'
#' @param dataPath A path containing the cell ranger processed data.
#' Under this path, folders 'filtered_feature_bc_matrix' and 'raw_feature_bc_matrix' exist generally.
#' @param statPath A path containing the results files of step 'runScStatistics'.
#' @param savePath A path to save the results files. If NULL, the 'statPath' will be used instead.
#' @param authorName A character string for authors name and will be shown in the report.
#' @param sampleName A character string giving a label for this sample.
#' @param bool.filter.cell A logical value indicating whether to filter the cells
#' according to the QC of 'scStatistics'.
#' @param bool.filter.gene A logical value indicating whether to filter the genes
#' according to the QC of 'scStatistics'.
#' @param anno.filter A vector indicating the types of genes to be filtered.
#' Must be some of c("mitochondrial", "ribosome", "dissociation")(default).
#' @param nCell.min An integer number used to filter gene. The default is 3.
#' Genes with the number of expressed cells less than this threshold will be filtered.
#' @param bgPercent.max A float number used to filter gene. The default is 0.001.
#' Genes with the background percentage larger than this threshold will be filtered.
#' @param vars.to.regress A vector indicating the variables to regress out in R package Seurat.
#' The default is c("nUMI", "mito.percent", "ribo.percent").
#' @param pc.use An integer number indicating the number of PCs to use as input features. The default is 30.
#' @param resolution A float number used in function 'FindClusters' in Seurat. The default is 0.8.
#' @param clusterStashName A character string used as the name of cluster identies. The default is "default".
#' @param show.features A list or vector for genes to be plotted in 'markerPlot'.
#' @param bool.add.features A logical value indicating whether to add default features to 'show.features' or not.
#' @param bool.runDiffExpr A logical value indicating whether to perform differential expressed analysis.
#' @param n.markers A integer indicating the number of differential expressed genes showed in the plot. The defalut is 5.
#' @param species A character string indicating what species the sample belong to.
#' Must be one of "human"(default) and "mouse".
#' @param hg.mm.mix  A logical value indicating whether the sample is a mix of
#' human cells and mouse cells(such as PDX sample).
#' If TRUE, the arguments 'hg.mm.thres' and 'mix.anno' should be set to corresponding values.
#' @param bool.runCellClassify A logical value indicating whether to predict the usual cell type. The default is TRUE.
#' @param ct.templates A list of vectors of several cell type templates.
#' The default is NULL and the templates prepared in this package will be used.
#' @param coor.names A vector indicating the names of two-dimension coordinate used in visualization.
#' @param bool.runMalignancy A logical value indicating whether to estimate malignancy.
#' @param cutoff A threshold used in the CNV inference.
#' @param p.value.cutoff A threshold to decide weather the bimodality distribution of malignancy score is significant.
#' @param bool.runGeneSets A logical value indicating whether to estimate gene sets signature scores.
#' @param geneSets A list of gene sets to be analyzed. The default is NULL and 50 hallmark gene sets from MSigDB will be used.
#' @param geneSet.method The method to be used in calculate gene set scores. Currently, only "average" and "GSVA" are allowed.
#' @param bool.runExprProgram A logical value indicating whether to run non-negative matrix factorization (NMF) to identify expression programs.
#' @param nmf.rank 	An integer of decomposition rank used in NMF.
#' @param genReport A logical value indicating whether to generate a .html/.md report (suggest to set TRUE).
#'
#' @return A results list with all useful objects used in the function.
#' @export
#'
#' @import Matrix knitr ggplot2 Seurat
#' @importFrom markdown markdownToHTML
#' @importFrom pheatmap pheatmap
#'
#' @examples
runScAnnotation <- function(dataPath, statPath, savePath = NULL,
                            authorName = NULL,
                            sampleName = "sc",
                            bool.filter.cell = T, bool.filter.gene = T,
                            anno.filter = c("mitochondrial", "ribosome", "dissociation"),
                            nCell.min = 3, bgPercent.max = 0.001,
                            vars.to.regress = c("nUMI", "mito.percent", "ribo.percent"),
                            pc.use = 30,
                            resolution = 0.8,
                            clusterStashName = "default",
                            show.features = NULL, bool.add.features = T,
                            bool.runDiffExpr = T,
                            n.markers = 5,
                            species = "human",
                            hg.mm.mix = F,
                            bool.runCellClassify = T,
                            ct.templates = NULL,
                            coor.names = c("tSNE_1", "tSNE_2"),
                            bool.runMalignancy = T,
                            cutoff = 0.1,
                            p.value.cutoff = 0.5,
                            bool.runGeneSets = T,
                            geneSets = NULL,
                            geneSet.method = "average",
                            bool.runExprProgram = T,
                            nmf.rank = 50,
                            genReport = T){

    message("[", Sys.time(), "] START: RUN scAnnotation")
    results <- as.list(environment())

    if(is.null(savePath)){
        savePath <- sataPath
    }
    suppressWarnings( dataPath <- normalizePath(dataPath, "/") )
    suppressWarnings( statPath <- normalizePath(statPath, "/") )
    suppressWarnings( savePath <- normalizePath(savePath, "/") )
    results[["dataPath"]] <- dataPath
    results[["statPath"]] <- statPath
    results[["savePath"]] <- savePath

    if(!dir.exists(file.path(savePath, "figures/"))){
        dir.create(file.path(savePath, "figures/"), recursive = T)
    }

    ## --------- filter data ---------
    message("[", Sys.time(), "] -----: cells and genes filtering")
    exprList <- getFilterData(dataPath = dataPath,
                              statPath = statPath,
                              savePath = savePath,
                              bool.filter.cell = bool.filter.cell,
                              bool.filter.gene = bool.filter.gene,
                              anno.filter = anno.filter,
                              nCell.min = nCell.min,
                              bgPercent.max = bgPercent.max,
                              hg.mm.mix = hg.mm.mix)
    results[["cell.manifest"]] = exprList$cell.manifest
    results[["gene.manifest"]] = exprList$gene.manifest
    results[["filter.thres"]] = exprList$filter.thres

    ## --------- seurat ---------
    t.results <- runSeurat(
        exprList = exprList,
        savePath = savePath,
        sampleName = sampleName,
        vars.to.regress = vars.to.regress,
        pc.use = pc.use,
        resolution = resolution,
        clusterStashName = clusterStashName,
        bool.runDiffExpr = bool.runDiffExpr
    )
    expr = t.results$expr
    cell.annotation = t.results$cell.annotation
    results[["diff.expr.genes"]] = t.results$diff.expr.genes
    rm(t.results)


    results[["seurat.plots"]] <- plotSeurat(
        expr = expr,
        cell.annotation = cell.annotation,
        show.features = show.features,
        bool.add.features = bool.add.features,
        coor.names = coor.names,

        bool.runDiffExpr = bool.runDiffExpr,
        diff.expr.genes = results[["diff.expr.genes"]],
        n.markers = n.markers,

        species = species,
        savePath = savePath
    )

    results[["DEplot.height"]] <- 0.5 + 0.1 * n.markers * length(unique(cell.annotation$Cluster))
    results[["markersPlot.height"]] <- 2 * ceiling(length(results[["seurat.plots"]]$ps.markers) / 4)



    ## --------- cell type ---------
    if(bool.runCellClassify){
        message("[", Sys.time(), "] -----: TME cell types annotation")
        t.results <- runCellClassify(expr, cell.annotation,
                                     savePath = savePath,
                                     ct.templates = ct.templates)

        expr <- t.results$expr
        cell.annotation <- t.results$cell.annotation
        results[["cellType.plot"]] <- t.results$p.results
        rm(t.results)
    }



    ## --------- malignancy ---------
    if(bool.runMalignancy){
        message("[", Sys.time(), "] -----: cells malignancy annotation")
        if(species != "human"){
            warning("To perform 'runMalignancy', the argument 'species' needs to be 'human'.\n")
            bool.runMalignancy = FALSE
        }else{
            t.results <- runMalignancy(dataPath, statPath, savePath,
                                       cell.annotation = cell.annotation,
                                       expr = expr,
                                       cutoff = cutoff, minCell = 3,
                                       p.value.cutoff = p.value.cutoff,
                                       hg.mm.mix = hg.mm.mix)
            expr <- t.results$expr
            cell.annotation <- t.results$cell.annotation
            results[["cnvList"]] <- t.results$cnvList
            results[["referScore"]] <- t.results$referScore
            results[["ju.exist.malign"]] <- t.results$ju.exist.malign
            results[["malign.thres"]] <- t.results$malign.thres
            results[["bimodal.pvalue"]] <- t.results$bimodal.pvalue
            results[["malign.plot"]] <- t.results$p.results
            rm(t.results)
        }
    }


    ## --------- gene sets ----------
    if(bool.runGeneSets){
        t.scores <- runGeneSets(expr = expr, geneSets = geneSets, method = geneSet.method)
        if(!is.null(t.scores)){
            cell.annotation <- cbind(cell.annotation, t.scores)
            rm(t.scores)

            bool.limit <- T
            if(geneSet.method == "GSVA"){
                bool.limit <- F
            }
            results[["geneSet.plot"]] <-
                plotGeneSet(cell.annotation, prefix = "GS__",
                            bool.limit = bool.limit, savePath = savePath)
        }else{
            bool.runGeneSets = FALSE
        }
    }


    ## ---------- expression programs ----------
    if(bool.runExprProgram){
        results[["exprProgram.results"]] <- runExprProgram(expr, rank = nmf.rank, savePath = savePath)
        results[["exprProgram.plot"]] <- plotExprProgram(H = results[["exprProgram.results"]]$H,
                                                         cell.annotation, savePath = savePath)
    }


    results[["expr"]] <- expr
    results[["cell.annotation"]] <- cell.annotation


    ## -------- save ---------
    saveRDS(expr, file = file.path(savePath, "expr.RDS"))
    write.table(cell.annotation, file = file.path(savePath, "cellAnnotation.txt"),
                quote = F, sep = "\t", row.names = F)

    if(genReport){
        message("[", Sys.time(), "] -----: report generating")
        if(!dir.exists(file.path(savePath, 'report-figures/'))){
            dir.create(file.path(savePath, 'report-figures/'), recursive = T)
        }
        suppressWarnings(
            knit(system.file("rmd", "main-scAnno.Rmd", package = "scCancer"),
                 file.path(savePath,'report-scAnno.md'), quiet = T)
        )
        markdownToHTML(file.path(savePath,'report-scAnno.md'),
                       file.path(savePath, 'report-scAnno.html'))
    }

    return(results)
}




#' genAnnoReport
#'
#' @param results A list generated by 'runScAnnotation'
#' @inheritParams runScAnnotation
#'
#' @return NULL
#' @export
#'
#' @examples
genAnnoReport <- function(results, savePath){
    message("[", Sys.time(), "] -----: report generating")
    results[['savePath']] <- normalizePath(savePath, "/")
    savePath <- normalizePath(savePath, "/")
    if(!dir.exists(file.path(savePath, 'report-figures/'))){
        dir.create(file.path(savePath, 'report-figures/'), recursive = T)
    }
    suppressWarnings(
        knit(system.file("rmd", "main-scAnno.Rmd", package = "scCancer"),
             file.path(savePath,'report-scAnno.md'), quiet = T)
    )
    markdownToHTML(file.path(savePath,'report-scAnno.md'),
                   file.path(savePath, 'report-scAnno.html'))
}

