filterCell <- function(cell.manifest, filter.thres){
    rownames(filter.thres) <- filter.thres$Index
    ix.cr <- which(cell.manifest$droplet.type == "cell")
    ix.thres <- getCellix(
        cell.manifest, filter.thres,
        arg = c("nUMI", "nGene", "mito.percent", "ribo.percent", "diss.percent")
    )
    return(cell.manifest$barcodes[intersect(ix.cr, ix.thres)])
}



filterGene <- function(gene.manifest,
                       anno.filter = c("mitochondrial", "ribosome", "dissociation"),
                       nCell.min = 3,
                       bgPercent.max = 1){
    ix.type <- which(!(gene.manifest$Annotation %in% anno.filter))
    ix.nCell <- which(gene.manifest$nCell >= nCell.min)
    ix <- intersect(ix.type, ix.nCell)

    if("bg.percent" %in% colnames(gene.manifest)){
        ix.bg <- which(gene.manifest$bg.percent <= bgPercent.max)
        ix <- intersect(ix, ix.bg)
    }else{
        if(bgPercent.max < 1){
            cat("- Warning in 'filterGene': Can not filter gene by 'bg.percent' due to the lack of background distribution.\n")
        }
    }
    return(gene.manifest$Symbol[ix])
}



# rmContamination <- function(expr.data, cell.manifest, contamination.fraction){
#     bg.bars <- subset(cell.manifest, nUMI >= 1 & nUMI <= 10)$barcodes
#     bg.sum <- Matrix::rowSums(expr.data[, bg.bars])
#     bg.percent <- bg.sum / sum(bg.sum)
#     soupProfile <- data.frame(row.names = rownames(expr.data),
#                          est = bg.percent,
#                          counts = bg.sum)
#
#     cell.manifest <- subset(cell.manifest, droplet.type == "cell")
#     sc <- list(toc = expr.data[, cell.manifest$barcodes],
#                metaData = data.frame(row.names = cell.manifest$barcodes,
#                                      nUMIs = cell.manifest$nUMI,
#                                      rho = contamination.fraction),
#                soupProfile = soupProfile)
#     class(sc) = c("list", "SoupChannel")
#     out = adjustCounts(sc, roundToInt = TRUE)
#
#     return(out)
# }


rmContamination <- function(statPath){
    if(file.exists(file.path(statPath, "soupx-object.RDS"))){
        soupx.obj = readRDS(file.path(statPath, "soupx-object.RDS"))
        out <- list(expr.data = adjustCounts(soupx.obj, roundToInt = FALSE, verbose = F),
                    contamination.frac = soupx.obj$fit$rhoEst)
    }else{
        out = NULL
    }
    return(out)
}



#' prepareSeurat
#'
#' According to the QC results of scStatistics, filter cells and genes.
#' Prepare a Seurat object.
#'
#' @inheritParams runScAnnotation
#'
#' @return A list of Seurat object and gene.manifest.
#' The Seurat object is after log-normalization, highly variable genes identification, scaling data.
#' @export
#'
prepareSeurat <- function(dataPath, statPath, savePath,
                          sampleName = "sc",
                          bool.filter.cell = T,
                          bool.filter.gene = T,
                          anno.filter = c("mitochondrial", "ribosome", "dissociation"),
                          nCell.min = 3, bgPercent.max = 1,
                          hg.mm.mix = F,
                          bool.rmContamination = T,
                          # contamination.fraction = NULL,
                          vars.add.meta = c("mito.percent", "ribo.percent", "diss.percent"),
                          vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent")){
    raw.data = T
    data.path <- get10Xpath(dataPath, raw.data = raw.data)
    if(is.null(data.path)){
        raw.data = F
        data.path <- get10Xpath(dataPath, raw.data = raw.data)
        if(is.null(data.path)){
            stop("Cannot find the raw data or filtered data.\n")
        }else{
            # bool.rmContamination <- F
            cat("- Warning in 'prepareSeurat': Cannot find the raw data, and use the filtered data instead.\n")
        }
    }

    cat("[", paste0(Sys.time()), "] -----: data preparation\n")
    expr.data <- Read10Xdata(data.dir = data.path, only.expr = T)
    rownames(expr.data) <- gsub("_", "-", rownames(expr.data))

    gene.manifest <- read.table(file.path(statPath, 'geneManifest.txt'),
                                header = T, sep = "\t", stringsAsFactors = F)
    rownames(gene.manifest) <- gene.manifest$Symbol

    cell.manifest <- read.table(file.path(statPath, 'cellManifest-all.txt'),
                                header = T, stringsAsFactors = F)
    rownames(cell.manifest) <- cell.manifest$barcodes

    filter.thres <- read.table(file.path(statPath, 'cell.QC.thres.txt'),
                               header = T, stringsAsFactors = F)

    if(bool.rmContamination){
        cat("[", paste0(Sys.time()), "] -----: contamination removing\n")
        soupx.out <- rmContamination(statPath = statPath)
        if(!is.null(soupx.out)){
            expr.data <- soupx.out$expr.data
            rownames(expr.data) <- gsub("_", "-", rownames(expr.data))

            cell.names <- colnames(expr.data)
            if (all(grepl(pattern = "\\-1$", x = cell.names))) {
                cell.names <- as.vector(x = as.character(
                    x = sapply(cell.names, FUN = ExtractField, field = 1, delim = "-")))
            }
            colnames(expr.data) <- cell.names

            contamination.frac <- soupx.out$contamination.frac
        }else{
            bool.rmContamination <- F
            contamination.frac <- NULL
            cat("- Warning in removing ambient RNA contamination: Cannot find the SoupX result file 'soupx-object.RDS'.\n")
        }
    }else{
        contamination.frac <- NULL
    }

    if(hg.mm.mix){
        rownames(expr.data) <- substr(rownames(expr.data), 6, 60)
    }

    if(bool.filter.cell){
        cells.select <- filterCell(cell.manifest,
                                   filter.thres = filter.thres)
    }else{
        cells.select <- subset(cell.manifest, droplet.type == "cell")$barcodes
    }

    # manifest.sub <- cell.manifest[cells.select,]
    # cells.select <- manifest.sub$barcodes[order(manifest.sub$nUMI, decreasing=T)[1:10000]]
    # print(length(cells.select))

    if(bool.filter.gene){
        genes.select <- filterGene(gene.manifest,
                                   anno.filter = anno.filter,
                                   nCell.min = nCell.min,
                                   bgPercent.max = bgPercent.max)
    }else{
        genes.select <- gene.manifest$Symbol
    }
    expr.data <- expr.data[genes.select, cells.select]


    cat("[", paste0(Sys.time()), "] -----: Seurat object creation\n")
    expr = CreateSeuratObject(counts = expr.data,
                              min.cells = 0,
                              min.features = 0,
                              project = sampleName)

    for(mv in vars.add.meta){
        if(!(mv %in% colnames(cell.manifest))){
            cat("- Warning in 'prepareSeurat': the", mv, " column is not in cell.manifest.\n")
        }else if(mv %in% colnames(expr@meta.data)){
            cat("- Warning in 'prepareSeurat': the", mv, " column has been in the Seurat object meta.data.\n")
        }else{
            tmp <- cell.manifest[cells.select, mv]
            names(tmp) <- cells.select
            expr[[mv]] <- tmp
        }
    }

    expr <- NormalizeData(object = expr,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000,
                          verbose = F)

    cat("[", paste0(Sys.time()), "] -----: highly variable genes\n")
    # expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000, verbose = F)
    # expr <- FindVariableFeatures(expr, selection.method = "vst",
    #                              nfeatures = min(10000, length(rownames(expr))), verbose = F)
    expr <- FindVariableFeatures(expr, selection.method = "vst",
                                 nfeatures = min(5000, length(rownames(expr))), verbose = F)

    cat("[", paste0(Sys.time()), "] -----: data scaling\n")
    # print(length(rownames(expr)))
    expr <- ScaleData(object = expr,
                      # features =  rownames(expr),
                      vars.to.regress = vars.to.regress,
                      verbose = F)

    gene.select.type <- rep("filter", dim(gene.manifest)[1])
    names(gene.select.type) <- gene.manifest$Symbol
    gene.select.type[genes.select] <- "keep"
    names(gene.select.type) <- gsub("_", "-", names(gene.select.type))
    gene.select.type[VariableFeatures(expr)] <- "hvg"
    names(gene.select.type) <- gene.manifest$Symbol
    gene.manifest$Select <- gene.select.type


    write.table(gene.manifest, file = file.path(savePath, "geneManifest.txt"),
                quote = F, sep = "\t", row.names = F)

    return(list(expr = expr,
                gene.manifest = gene.manifest,
                bool.rmContamination = bool.rmContamination,
                contamination.frac = contamination.frac))
}




#' runSeurat
#'
#' Perform usual Seurat step and cell type prediction.
#'
#' @param expr A Seurat object return by prepareSeurat.
#' @param comb.method The method used in combining samples. It worked only for multi-sample analysis.
#' @inheritParams runScAnnotation
#'
#' @return A list containing a Seurat object, differential expressed genes and annotation information for cells.
#' @export
#'
runSeurat <- function(expr,
                      savePath,
                      npcs = 50,
                      pc.use = 30, resolution = 0.8,
                      clusterStashName = "default",
                      bool.runDiffExpr = T,
                      comb.method = NULL){

    if(!dir.exists(file.path(savePath))){
        dir.create(file.path(savePath), recursive = T)
    }

    if(!is.null(comb.method)){
        if(comb.method == "Harmony"){
            reduction.type = "harmony"
        }else if(comb.method == "LIGER"){
            reduction.type = "inmf"
            pc.use <- min(pc.use, ncol(expr@reductions$inmf@cell.embeddings))
        }else{
            cat("[", paste0(Sys.time()), "] -----: PCA\n")
            expr <- RunPCA(expr, verbose = F)
            reduction.type <- "pca"
        }
    }else{
        cat("[", paste0(Sys.time()), "] -----: PCA\n")
        expr <- RunPCA(expr, npcs = npcs, verbose = F)
        reduction.type <- "pca"
    }

    cat("[", paste0(Sys.time()), "] -----: clustering\n")
    expr <- FindNeighbors(expr, reduction = reduction.type, dims = 1:pc.use, verbose = F)
    expr <- FindClusters(expr, resolution = resolution, verbose = F)
    expr[[clusterStashName]] <- as.numeric(Idents(object = expr))

    if(is.null(comb.method)){
        cat("[", paste0(Sys.time()), "] -----: tSNE\n")
        expr <- RunTSNE(object = expr, dims = 1:pc.use, reduction = reduction.type)
    }else{
        if(comb.method != "LIGER"){
            cat("[", paste0(Sys.time()), "] -----: tSNE\n")
            expr <- RunTSNE(object = expr, dims = 1:pc.use, reduction = reduction.type)
        }
    }

    cat("[", paste0(Sys.time()), "] -----: UMAP\n")
    suppressWarnings(
        tryCatch(expr <- RunUMAP(expr, dims = 1:pc.use, reduction = reduction.type, verbose = F),
                 error = function(err) {
                     cat("Error in 'RunUMAP': Please use 'pip install umap-learn' to install UMAP firstly.\n")}
        )
    )

    if(bool.runDiffExpr){
        cat("[", paste0(Sys.time()), "] -----: differential expression analysis\n")
        if(!dir.exists(file.path(savePath, "diff.expr.genes"))){
            dir.create(file.path(savePath, "diff.expr.genes"), recursive = T)
        }
        diff.expr.genes <- FindAllMarkers(expr, only.pos = TRUE,
                                          min.pct = 0.25,
                                          logfc.threshold = 0.25,
                                          verbose = F)
        saveRDS(diff.expr.genes, file.path(savePath, "DE-Genes.RDS"))
        # write.table(diff.expr.genes[, c("cluster", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")],
        #             file = file.path(savePath, "diff.expr.genes.txt"), quote = F, sep = "\t", row.names = F)
        if("avg_logFC" %in% colnames(diff.expr.genes)){
            de.FC <- "avg_logFC"   # packageVersion("Seurat") < 4
        }else if("avg_log2FC" %in% colnames(diff.expr.genes)){
            de.FC <- "avg_log2FC"   # packageVersion("Seurat") >= 4
        }else{
            stop("Couldn't find the differential expression results.\n")
        }
        diff.expr.genes <- diff.expr.genes[, c("cluster", "gene", "p_val", de.FC, "pct.1", "pct.2", "p_val_adj")]
        diff.expr.genes$cluster <- as.numeric(diff.expr.genes$cluster)
        nCluster <- length(unique(diff.expr.genes$cluster))
        for(ci in 1:nCluster){
            cur.diff.genes <- subset(diff.expr.genes, cluster == ci)
            cur.diff.genes <- cur.diff.genes[order(cur.diff.genes[[de.FC]], decreasing = T), ]
            write.csv(cur.diff.genes,
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
    if("umap" %in% names(expr@reductions)){
        cell.annotation <- cbind(cell.annotation, expr@reductions$umap@cell.embeddings)
    }
    # cell.annotation$Cluster <- factor(expr@meta.data[[clusterStashName]],
    #                                   levels = 0:(length(unique(expr@meta.data[[clusterStashName]]))-1))
    cell.annotation$Cluster <- factor(expr@meta.data[[clusterStashName]])

    return(list(expr = expr,
                diff.expr.genes = diff.expr.genes,
                cell.annotation = cell.annotation))
}



singleGenePlot <- function(expr.data, gene,
                           coor.df, coor.names = c("UMAP_1", "UMAP_2"),
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
        new.label <- paste0(gene, " (0/", dim(expr.data)[2], ")")
        p <- ggplot(coor.df, aes(x = .data[[coor.names[1]]], y = .data[[coor.names[2]]])) +
            geom_point(shape = 21, size = 0.8, stroke = 0.2, color = "grey", fill = "white") +
            annotate("text", x = minx, y = miny, label = new.label, hjust = 0, vjust = 0,
                     size = font.size, fontface = 'italic') +
            theme_classic()
    }else{
        coor.df$cur.value <- expr.data[gene, ]
        new.label <- paste0(gene, " (", sum(expr.data[gene, ] > 0), "/", dim(expr.data)[2], ")")
        p <- ggplot(coor.df, aes(x = .data[[coor.names[1]]], y = .data[[coor.names[2]]],
                                 fill = cur.value)) +
            geom_point(shape = 21, size = 0.8, stroke = 0.2, color = "grey") +
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



#' markerPlot
#'
#' Generate plots of interested genes' expression profile.
#'
#' @param expr.data A matrix of expression (gene by cell)
#' @param coor.df A data.frame which contains cells' 2D coordinates.
#' @param features A vector of genes to plot.
#' @param add A logical value indicating whether to present the default markers.
#' @param font.size The size of labels.
#' @param color The color of point.
#' @inheritParams runScAnnotation
#'
#' @return A list of ggplot obejects for each maker genes.
#' @export
#'
markerPlot <- function(expr.data, coor.df, coor.names = c("UMAP_1", "UMAP_2"),
                       features = NULL, add = T,
                       species = "human",
                       font.size = 4, color = "blue"){
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



#' pointDRPlot
#'
#' Plot scatter for cells.
#'
#' @param cell.annotation A data.frame of cells' annotation containing the cells' coordinates and index to be colored.
#' @param sel.clusters An array of selected clusters to present. (The default is NULL and all clusters will be used.)
#' @param value The column name of cell.annotation, which is mapped to the colors of points.
#' @param colors An array of colors used to show the gredients or type of points. If NULL, the default colors will be used.
#' @param discrete A logical value indicating whether the value column is discrete or not.
#' @param limit.quantile A quantile threshold to limit the data and reduce the influence of outliers.
#' @param point.type A number indicating the shape type of points. "1" (default) means the point has a lightgrey border, and "2" means not.
#' @param legend.position The position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector).
#' @param legend.title The title of legends.
#' @inheritParams runScAnnotation
#'
#' @return A ggplot object for the scatter plot.
#' @export
#'
pointDRPlot <- function(cell.annotation, value,
                        sel.clusters = NULL,
                        coor.names = c("UMAP_1", "UMAP_2"),
                        colors = NULL,
                        discrete = T,
                        limit.quantile = 0,
                        point.type = 1,
                        legend.position = "right",
                        legend.title = NULL){
    if(!(all(coor.names %in% colnames(cell.annotation)))){
        stop("Error in 'pointDRPlot': 'coor.names' ", coor.names, " is not in colnames(cell.annotation).\n")
    }
    if(!(point.type %in% c(1, 2))){
        stop("Error in 'pointDRPlot': 'point.type' ", point.type, " is not allowed.\n")
    }

    if(is.null(legend.title)){
        legend.title <- value
    }
    if(is.null(colors)){
        if(discrete){
            colors <- getDefaultColors(length(unique(cell.annotation[[value]])))
        }else{
            colors <- c("white", "red")
        }
    }

    ratio <- diff(range(cell.annotation[, coor.names[1]])) / diff(range(cell.annotation[, coor.names[2]]))

    fill.value <- cell.annotation[[value]]
    if(!discrete){
        low.thres <- quantile(cell.annotation[[value]], limit.quantile)
        up.thres <- quantile(cell.annotation[[value]], 1 - limit.quantile)
        fill.value <- ifelse(fill.value < low.thres, low.thres,
                             ifelse(fill.value > up.thres, up.thres, fill.value))
    }

    coor.label <- coor.names
    if(all.equal(coor.label, c("tSNE_1", "tSNE_2")) == TRUE){
        coor.label[1] <- "t-SNE 1"
        coor.label[2] <- "t-SNE 2"
    }else if(all.equal(coor.label, c("UMAP_1", "UMAP_2")) == TRUE){
        coor.label[1] <- "UMAP 1"
        coor.label[2] <- "UMAP 2"
    }

    if(!is.null(sel.clusters)){
        sel.cell <- cell.annotation$Cluster %in% sel.clusters
        if(point.type == 1){
            p <- ggplot() +
                geom_point(cell.annotation[!sel.cell, ],
                           mapping = aes(x = .data[[coor.names[1]]],
                                         y = .data[[coor.names[2]]]),
                           fill = "white",
                           shape = 21, size = 0.8, stroke = 0.2, color = "white")
        }else if(point.type == 2){
            p <- ggplot() +
                geom_point(cell.annotation[!sel.cell, ],
                           mapping = aes(x = .data[[coor.names[1]]],
                                         y = .data[[coor.names[2]]],
                                         color = "white"),
                           shape = 16, size = 0.3)
        }
    }else{
        sel.cell <- rep(TRUE, dim(cell.annotation)[1])
        p <- ggplot()
    }

    if(point.type == 1){
        # color = "lightgrey"
        p <- p +
            geom_point(cell.annotation[sel.cell, ],
                       mapping = aes(x = .data[[coor.names[1]]],
                                     y = .data[[coor.names[2]]],
                                     fill = fill.value[sel.cell]),
                       shape = 21, size = 0.8, stroke = 0.2, color = "lightgrey") +
            coord_fixed(ratio = ratio) +
            ggplot_config(base.size = 6) +
            labs(x = coor.label[1], y = coor.label[2]) +
            labs(fill = legend.title) +
            theme(legend.position = legend.position)

        if(discrete){
            p <- p + scale_fill_manual(values = colors,
                                       guide = guide_legend(override.aes = list(size = 3),
                                                            keywidth = 0.1,
                                                            keyheight = 0.15,
                                                            default.unit = "inch"))
        }else{
            p <- p + scale_fill_gradientn(colors = colors)
        }
    }else if(point.type == 2){
        p <- p + geom_point(cell.annotation[sel.cell, ],
                       mapping = aes(x = .data[[coor.names[1]]],
                                     y = .data[[coor.names[2]]],
                                     color = fill.value[sel.cell]),
                       shape = 16, size = 0.3) +
            coord_fixed(ratio = ratio) +
            ggplot_config(base.size = 6) +
            labs(x = coor.label[1], y = coor.label[2]) +
            labs(color = legend.title) +
            theme(legend.position = legend.position)

        if(discrete){
            p <- p + scale_color_manual(values = colors,
                                        guide = guide_legend(override.aes = list(size = 3),
                                                             keywidth = 0.1,
                                                             keyheight = 0.15,
                                                             default.unit = "inch"))
        }else{
            p <- p + scale_color_gradientn(colors = colors)
        }
    }
    return(p)
}



#' clusterBarPlot
#'
#' @param cell.annotation A data.frame of cells' annotation containing the cells' Cluster and other information to be colored.
#' @param sel.col The column name of cell.annotation, which indicating the type of cells.
#' @param cell.colors An array of colors used to show the cells' type. If NULL, the default colors will be used.
#' @param legend.title The title of legends. If NULL, the value of "sel.col" will be used.
#' @param legend.position The position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector).
#' @param legend.ncol The number of column of legends.
#'
#' @return A bar plot.
#' @export
#'
clusterBarPlot <- function(cell.annotation, sel.col = "Cell.Type", cell.colors = NULL,
                           legend.title = NULL, legend.position = "bottom", legend.ncol = NULL){
    if(is.null(legend.title)){
        legend.title <- sel.col
    }
    if(is.null(cell.colors)){
        cell.colors <- getDefaultColors(length(unique(cell.annotation[[sel.col]])))
    }

    bar.df <- melt(table(cell.annotation[c(sel.col, "Cluster")]))
    # bar.df$Cluster = factor(bar.df$Cluster,
    #                         levels = 0:(length(unique(bar.df$Cluster))-1))
    bar.df$Cluster = factor(bar.df$Cluster)

    if(legend.position == "bottom" & is.null(legend.ncol)){
        legend.ncol <- 3
    }
    p <- ggplot(bar.df, aes(x = Cluster, y = value, fill = .data[[sel.col]])) +
        geom_bar(stat = "identity") +
        ggplot_config(base.size = 6) +
        scale_fill_manual(values = cell.colors) +
        labs(y = "Number of cells") +
        guides(fill = guide_legend(override.aes = list(size=1),
                                   title = legend.title, ncol = legend.ncol)) +
        theme(legend.position = legend.position)

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
#' @param expr A Seurat object.
#' @param cell.annotation A data.frame of cells' annotation.
#' @param bool.plotHVG A logical value indicating Whehter to plot highly variable genes.
#' @param diff.expr.genes A data.frame of differential expressed genes.
#' @inheritParams runScAnnotation
#'
#' @return A list of all plots generated by Seurat analyses.
#' @export
#'
plotSeurat <- function(expr,
                       cell.annotation = cell.annotation,
                       show.features = NULL, bool.add.features = T,
                       coor.names = c("UMAP_1", "UMAP_2"),
                       bool.plotHVG = T,
                       bool.runDiffExpr = T,
                       diff.expr.genes = NULL, n.markers = 5,
                       species = "human",
                       savePath){

    cat("[", paste0(Sys.time()), "] -----: Seurat plotting and saving\n")

    if(!dir.exists(file.path(savePath, "figures/singleMarkerPlot/"))){
        dir.create(file.path(savePath, "figures/singleMarkerPlot/"), recursive = T)
    }

    p.results <- list()

    # message(sprintf('------p.hvg------'))
    if(bool.plotHVG){
        top10 <- head(VariableFeatures(expr), 10)
        suppressWarnings(
            plot1 <- VariableFeaturePlot(expr, cols = c("grey", "#ec7d89"))
        )
        suppressWarnings(
            p.results[["p.hvg"]] <- LabelPoints(plot = plot1, points = top10, repel = TRUE,
                                                fontface = "italic") + NoLegend()
        )
    }


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
    clusters <- unique(cell.annotation$Cluster)
    clusters <- sort(clusters)

    def.colors <- getDefaultColors(n = length(clusters))
    cluster.colors <- c()
    for(i in 1:length(clusters)){
        cluster.colors[as.character(clusters[i])] = def.colors[i]
    }

    p.results[["p.cluster.tsne"]] <- pointDRPlot(cell.annotation, value = "Cluster",
                                                 coor.names = c("tSNE_1", "tSNE_2"),
                                                 colors = cluster.colors,
                                                 legend.title = "Cluster")

    if(all(c("UMAP_1", "UMAP_2") %in% colnames(cell.annotation))){
        p.results[["p.cluster.umap"]] <- pointDRPlot(cell.annotation, value = "Cluster",
                                                     coor.names = c("UMAP_1", "UMAP_2"),
                                                     colors = cluster.colors,
                                                     legend.title = "Cluster")
    }else{
        p.results[["p.cluster.umap"]] <- NULL
    }


    if(bool.runDiffExpr && !(is.null(diff.expr.genes))){
        # message(sprintf('------p.DE.heatmap------'))
        if("avg_logFC" %in% colnames(diff.expr.genes)){
            de.FC <- "avg_logFC"   # packageVersion("Seurat") < 4
        }else if("avg_log2FC" %in% colnames(diff.expr.genes)){
            de.FC <- "avg_log2FC"   # packageVersion("Seurat") >= 4
        }
        top.genes <- diff.expr.genes %>% group_by(cluster) %>% top_n(n = n.markers, wt = .data[[de.FC]])
        top.genes <- top.genes[order(top.genes$cluster, top.genes[[de.FC]], decreasing = c(F, T)), ]
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
                     show_colnames = F,
                     # legend = F,
                     silent = T)
    }

    # message(sprintf('------save images------'))
    ggsave(filename = file.path(savePath, "figures/hvg.png"), p.results[["p.hvg"]],
           width = 8, height = 4, dpi = 300)

    for(gene in names(p.results[["ps.markers"]])){
        ggsave(filename = paste0(savePath, "/figures/singleMarkerPlot/", gene, ".png"),
               p.results[["ps.markers"]][[gene]], width = 2, height = 2, dpi = 300)
    }

    if(length(p.results[["ps.markers"]]) > 0){
        markersPlot.height <- 2 * ceiling(length(p.results[["ps.markers"]]) / 4)
        ggsave(filename = file.path(savePath, "figures/markers-all.png"),
               p.results[["p.markers.all"]], width = 8, height = markersPlot.height, dpi = 300)
    }

    ggsave(filename = file.path(savePath, "figures/cluster-point-tsne.png"),
           p.results[["p.cluster.tsne"]], width = 6, height = 5, dpi = 300)

    if(!is.null(p.results[["p.cluster.umap"]])){
        ggsave(filename = file.path(savePath, "figures/cluster-point-umap.png"),
               p.results[["p.cluster.umap"]], width = 6, height = 5, dpi = 300)
    }

    if(bool.runDiffExpr && !(is.null(diff.expr.genes))){
        DEplot.height <- 0.5 + 0.1 * n.markers * length(unique(cell.annotation$Cluster))
        ggsave(filename = file.path(savePath, "figures/DE-heatmap.png"),
               p.results[["p.de.heatmap"]], width = 8, height = DEplot.height, dpi = 300)
    }

    # saveRDS(cell.annotation, file = file.path(savePath, "cell.annotation.RDS"))

    return(p.results)
}



#' runDoublet
#'
#' @param expr A Seurat object.
#' @param method The method to estimate doublet score. The default is "cxds".
#' @param pc.use An integer number indicating the number of PCs to use as input features. The default is 30.
#'
#' @return  An array of doublet scores.
#' @export
#'
#' @importFrom scds cxds bcds
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
runDoublet <- function(expr, method = "cxds", pc.use = 30){
    if(method == "cxds"){
        expr <- SingleCellExperiment(
            assays = list(counts = GetAssayData(expr, slot = "counts")))
        expr <- cxds(expr)
        doublet.score <- expr$cxds_score
    }else if(method == "bcds"){
        expr <- SingleCellExperiment(
            assays = list(counts = GetAssayData(expr, slot = "counts")))
        expr <- bcds(expr)
        doublet.score <- expr$bcds_score
    # }else if(method == "DoubletFinder"){
    #     cluster.results <- expr@meta.data$default
    #     homotypic.prop <- modelHomotypic(cluster.results)
    #     nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    #     expr <- doubletFinder_v3(expr, PCs = 1:pc.use,
    #                              pN = 0.25, pK = 0.09,
    #                              nExp = nExp_poi,
    #                              reuse.pANN = FALSE, sct = FALSE)
    #     doublet.score <- expr@meta.data[, paste0("pANN_0.25_0.09_", nExp)]
    }
    names(doublet.score) <- colnames(expr)
    return(doublet.score)
}



#' predCellType
#'
#' @param X.test A cells expression matrix (row for genes, column for cells).
#' @param ct.templates A list of gene weight vectors for each cell type.
#' @inheritParams runScAnnotation
#'
#' @return A list of predicted cell types and the relative correlations.
#' @export
#'
predCellType <- function(X.test, ct.templates = NULL, species = "human"){
    if(is.null(ct.templates)){
        ct.templates <- readRDS(system.file("rds", "cellTypeTemplates.RDS", package = "scCancer"))
        if(species == "mouse"){
            for(cur.ct in names(ct.templates)){
                cur.temp <- ct.templates[[cur.ct]]
                cur.genes <- getMouseGene(names(cur.temp), bool.name = T)
                cur.temp <- cur.temp[names(cur.genes)]
                names(cur.temp) <- cur.genes
                ct.templates[[cur.ct]] <- cur.temp
            }
        }
    }

    # cor.df <- list()
    # for(i in 1:length(ct.templates)){
    #     type <- names(ct.templates)[i]
    #     common.gene <- intersect(rownames(X.test), names(ct.templates[[type]]))
    #     s.cor <- cor(as.matrix(X.test[common.gene, ]),
    #                  ct.templates[[type]][common.gene],
    #                  method = "spearman")
    #     s.cor[is.na(s.cor)] <- 0
    #     cor.df[[type]] <- s.cor
    # }
    cor.df <- lapply(1:length(ct.templates), function(i){
        type <- names(ct.templates)[i]
        common.gene <- intersect(rownames(X.test), names(ct.templates[[type]]))
        s.cor <- cor(as.matrix(X.test[common.gene, ]),
                     ct.templates[[type]][common.gene],
                     method = "spearman")
        s.cor[is.na(s.cor)] <- 0
        return(s.cor)
    })
    names(cor.df) <- names(ct.templates)
    cor.df <- as.data.frame(cor.df)

    type.pred <- colnames(cor.df)[unlist(apply(cor.df, 1, which.max))]
    # type.pred[rowSums(cor.df > 0.25) == 0] <- "Unknown"

    thres <- apply(cor.df, 2, quantile, probs = 0.05)
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
#' @param expr A Seurat object.
#' @param cell.annotation A data.frame of cells' annotation.
#' @inheritParams runScAnnotation
#'
#' @return A list of updated Seurat object, cell.annotation, and the plots for cell type annotation.
#' @export
#'
runCellClassify <- function(expr, cell.annotation, coor.names = c("UMAP_1", "UMAP_2"),
                            savePath, ct.templates = NULL, species = "human"){
    if(!("Cell.Type" %in% names(cell.annotation))){
        # message("[", Sys.time(), "] -----: TME cell types annotation")
        cat("[", paste0(Sys.time()), "] -----: TME cell types annotation\n")
        t.results <- predCellType(X.test = expr@assays$RNA@data,
                                  ct.templates = ct.templates, species = species)

        expr[["Cell.Type"]] <- t.results$type.pred
        for(score in names(t.results$cor.df)){
            expr[[score]] <- t.results$cor.df[[score]]
        }

        cell.annotation$Cell.Type <- t.results$type.pred
        cell.annotation <- cbind(cell.annotation, t.results$cor.df)
    }else{
        cat("[", paste0(Sys.time()), "] -----: TME cell types combination\n")
    }

    # cell.colors <- c(
    #     "T.cells.CD4" = "#07a2a4",
    #     "T.cells.CD8" = "#9a7fd1",
    #     "B.cells" = "#588dd5",
    #     "NK.cells" = "#f5994e",
    #     "Myeloid.cells" = "#c05050",
    #     "Endothelial" = "#59678c",
    #     "Fibroblast" = "#c9ab00",
    #     "Epithelial" = "#7eb00a",
    #     "Unknown" = "grey")
    # cti = 1
    # new.types <- setdiff(cell.annotation$Cell.Type, names(cell.colors))
    # for(ct in new.types){
    #     cell.colors[ct] <- getDefaultColors(n = length(new.types), type = 3)[cti]
    #     cti = cti + 1
    # }

    cell.colors <- getCellTypeColor(cell.types = unique(cell.annotation$Cell.Type))

    # message(sprintf('------p.type------'))
    p.type <- pointDRPlot(cell.annotation, value = "Cell.Type",
                          coor.names = coor.names,
                          colors = cell.colors,
                          legend.title = "Cell type")

    # message(sprintf('------p.bar------'))
    p.bar <- clusterBarPlot(cell.annotation = cell.annotation,
                            cell.colors = cell.colors,
                            sel.col = "Cell.Type",
                            legend.title = "Cell type")

    ggsave(filename = file.path(savePath, "figures/cellType-point.png"),
           p.type, width = 5.2, height = 4, dpi = 300)
    ggsave(filename = file.path(savePath, "figures/cellType-bar.png"),
           p.bar, width = 6, height = 4, dpi = 300)

    return(list(expr = expr,
                cell.annotation = cell.annotation,
                p.results = list(p.type = p.type,
                                 p.bar = p.bar)))
}



#' getTumorCluster
#'
#' Identify tumor clusters according to the results of cell type prediction and cell malignancy estimatation.
#'
#' @param cell.annotation A data.frame of cells' annotation containing predicted cell typea and estimated cell malignant type.
#' @param epi.thres A threshold for epithelial cell percent to decide putative tumor clusters.
#' @param malign.thres A threshold for malignant cell percent to decide putative tumor clusters.
#'
#' @return A list of identified tumor clusters. If no clusters are found, return NULL.
#' @export
#'
getTumorCluster <- function(cell.annotation, epi.thres = 0.6, malign.thres = 0.8){
    bool.sel <- F
    epithe.clusters <- c()
    if("Cell.Type" %in% names(cell.annotation)){
        for(cluster in unique(cell.annotation$Cluster)){
            tmp <- subset(cell.annotation, Cluster == cluster)
            percent.epithe <- sum(tmp$Cell.Type == "Epithelial") / dim(tmp)[1]
            if(percent.epithe > epi.thres){
                epithe.clusters <- c(epithe.clusters, cluster)
            }
        }
        bool.sel <- T
    }else{
        epithe.clusters <- unique(cell.annotation$Cluster)
    }
    malign.clusters <- c()
    if("Malign.type" %in% names(cell.annotation)){
        for(cluster in unique(cell.annotation$Cluster)){
            tmp <- subset(cell.annotation, Cluster == cluster)
            percent.malign <- sum(tmp$Malign.type == "malignant") / dim(tmp)[1]
            if(percent.malign > malign.thres){
                malign.clusters <- c(malign.clusters, cluster)
            }
        }
        bool.sel <- T
    }else{
        malign.clusters <- unique(cell.annotation$Cluster)
    }

    tumor.clusters <- intersect(epithe.clusters, malign.clusters)
    if(length(tumor.clusters) == 0 || !bool.sel){
        cat("- Warning in 'getTumorCluster': Could not identify tumor clusters.\n")
        return(NULL)
    }
    tumor.clusters <- as.numeric(tumor.clusters)
    tumor.clusters <- tumor.clusters[order(tumor.clusters)]
    return(tumor.clusters)
}



#' runCellCycle
#'
#' Estimate cell cycle scores.
#'
#' @param expr A Seurat object.
#' @inheritParams runScAnnotation
#'
#' @return An array of cell cycle scores.
#' @export
#'
runCellCycle <- function(expr, species = "human"){
    cat("[", paste0(Sys.time()), "] -----: cell cycle score estimation\n")
    cellCycle.genes <- read.table(system.file("txt", "cellCycle-genes.txt", package = "scCancer"),
                                  header = F, stringsAsFactors = F)$V1
    if(species == "mouse"){
        cellCycle.genes <- getMouseGene(cellCycle.genes)
    }
    suppressWarnings(
        expr <- AddModuleScore(expr, features = list(cellCycle.genes),
                               assay = "RNA", name = "cellCycle")
    )
    return(expr[["cellCycle1"]]$cellCycle1)
}



#' runStemness
#'
#' Estimate cell stemness according to the Spearman correlation with stemness signature.
#'
#' @param X An expression matrix of gene by cell to estimate stemness.
#' @param stem.sig An array of stemness signature. The default is NULL, and a prepared signature will be used.
#' @inheritParams runScAnnotation
#'
#' @return An array of cell stemness scores.
#' @export
#'
runStemness <- function(X, stem.sig = NULL, species = "human"){
    cat("[", paste0(Sys.time()), "] -----: stemness score calculation\n")
    if(is.null(stem.sig)){
        stem.sig.file <- system.file("txt", "pcbc-stemsig.tsv", package = "scCancer")
        stem.sig <- read.delim(stem.sig.file, header = FALSE, row.names = 1)
        if(species == "mouse"){
            sig.genes <- getMouseGene(rownames(stem.sig), bool.name = T)
            stem.sig <- stem.sig[names(sig.genes), , drop=F]
            rownames(stem.sig) <- sig.genes
        }
    }

    common.genes <- intersect(rownames(stem.sig), rownames(X))
    X <- X[common.genes, ]
    stem.sig <- stem.sig[common.genes, ]

    s <- apply(X, 2, function(z) {cor(z, stem.sig, method = "sp", use = "complete.obs")})
    names(s) <- colnames(X)

    s <- s - min(s)
    s <- s / max(s)

    return(s)
}



#' getDefaultGeneSets
#'
#' @inheritParams runScAnnotation
#'
#' @return A list of gene sets (50 hallmark gene sets).
#' @export
#'
getDefaultGeneSets <- function(species = "human"){
    geneSets <- readLines(system.file("txt", "hallmark-pathways.txt", package = "scCancer"))
    geneSets <- strsplit(geneSets, "\t")
    geneSets <- as.data.frame(geneSets, stringsAsFactors = F)
    colnames(geneSets) <- geneSets[1, ]
    geneSets <- as.list(geneSets[2, ])
    geneSets <- sapply(geneSets, function(x) strsplit(x, ", "))
    if(species == "mouse"){
        for(set.name in names(geneSets)){
            geneSets[[set.name]] <- getMouseGene(geneSets[[set.name]])
        }
    }
    return(geneSets)
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
runGeneSets <- function(expr, geneSets, method = "average"){
    cat("[", paste0(Sys.time()), "] -----: gene set signatures analysis\n")
    if(class(geneSets) != "list"){
        cat("- Warning in 'runGeneSets': The 'geneSets' should be a list of several gene sets.\n")
        return(NULL)
    }
    if(method == "average"){
        suppressWarnings(
            expr <- AddModuleScore(expr, features = geneSets,
                                   assay = "RNA", name = "geneSets")
        )
        t.scores <- expr[[paste0("geneSets", 1:length(geneSets))]]
        t.scores <- scale(t.scores)
    }else if(method == "GSVA"){
        tmp.data <- as.matrix(GetAssayData(object = expr, slot = "scale.data"))
        tmp.data <- tmp.data[VariableFeatures(expr), ]
        t.scores <- gsva(tmp.data, geneSets)
        t.scores <- t(t.scores)
    }else{
        cat("- Warning in 'runGeneSets': The 'method' ", method, " is not allowed.\n", sep = "")
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
                  # legend = F,
                  silent = T)

    if(!is.null(savePath)){
        geneSetPlot.height <- 0.5 + 0.11 * length(gs.name)
        ggsave(filename = file.path(savePath, "figures/geneSet-heatmap.png"),
               p, width = 10, height = geneSetPlot.height, dpi = 300)
    }

    return(p)
}




#' runExprProgram
#'
#' Perform non-negative matrix factorization (NMF) to identify expression programs.
#'
#' @param expr A Seurat object.
#' @param rank An integer of decomposition rank used in NMF.
#' @param sel.clusters A vector of selected clusters to analyze. The default is NULL and all clusters will be used.
#' @inheritParams runScAnnotation
#'
#' @return A list of decomposed matrixes (W and H), and the relative genes of each programs.
#' @export
#'
#' @importFrom methods as
#'
runExprProgram <- function(expr, rank = 50, sel.clusters = NULL, clusterStashName = "default", savePath = NULL){
    cat("[", paste0(Sys.time()), "] -----: expression programs analysis\n")

    data <- as(object = expr[["RNA"]]@data, Class = "TsparseMatrix")
    if(!is.null(sel.clusters)){
        data <- data[, expr@meta.data[[clusterStashName]] %in% sel.clusters]
    }

    if(rank > min(dim(data))){
        rank <- min(rank, min(dim(data)))
        cat("- Warning in 'runExprProgram':
            The input rank is larger than the size of data for NMF, and use the minimum of them instead.\n")
    }

    ave.data <-  Matrix::rowSums(data) / Matrix::rowSums(data > 0)

    data@x <- data@x - ave.data[data@i + 1]
    data@x[data@x < 0] <- 0

    nmf.results <- nnmf(as.matrix(data), k = rank, verbose = 0)

    W <- nmf.results$W
    colnames(W) <- paste0("p", 1:dim(W)[2])
    H <- nmf.results$H
    rownames(H) <- paste0("p", 1:dim(H)[1])

    all.genes <- rownames(W)
    sel.W <- (W > quantile(W, 1 - 50/dim(W)[1]))
    for(pi in 1:dim(W)[2]){
        if(sum(sel.W[, pi]) > 10){
            tmp <- data.frame(program = colnames(W)[pi], gene = all.genes[sel.W[, pi]], value = W[sel.W[, pi], pi])
            tmp <- tmp[order(tmp$value, decreasing = T), ]
        }else{
            tmp <- data.frame(program = colnames(W)[pi], gene = all.genes, value = W[, pi])
            tmp <- tmp[order(tmp$value, decreasing = T), ]
            tmp <- tmp[1:10, ]
        }
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
#' @param sel.clusters A vector of selected clusters to analyze. The default is NULL and all clusters will be used.
#' @inheritParams runScAnnotation
#'
#' @return A heatmap for cells' expression programs.
#' @export
#' @importFrom NNLM nnmf
#'
plotExprProgram <- function(H, cell.annotation, bool.limit = T, sel.clusters = NULL, savePath = NULL){
    if(bool.limit){
        up.bound <- quantile(as.matrix(H), 0.995)
        H <- limitData(H, max = up.bound)
    }
    if(!is.null(sel.clusters)){
        cell.annotation <- subset(cell.annotation, Cluster %in% sel.clusters)
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
        exprProgPlot.height <- 0.5 + 0.11 * dim(H)[1]
        ggsave(filename = file.path(savePath, "figures/exprProgram-heatmap.png"),
               p, width = 10, height = exprProgPlot.height, dpi = 300)
    }

    # clusters <- unique(cell.annotation$Cluster)
    # clusters <- sort(clusters)
    #
    # def.colors <- getDefaultColors(n = length(clusters))
    # cluster.colors <- c()
    # for(i in 1:length(clusters)){
    #     # cluster.colors[as.character(clusters[i])] = def.colors[clusters[i]]
    #     cluster.colors[i] = def.colors[clusters[i]]
    # }
    # cluster.colors = list(Cluster = cluster.colors)
    # ha <- HeatmapAnnotation(df = data.frame(Cluster = cell.annotation$Cluster),
    #                         name ="Cluster", col = cluster.colors)
    #
    # p <- Heatmap(H, name = "H",
    #              col = c("#f9fcfb","#009b45"),
    #              top_annotation = ha,
    #              column_split = cell.annotation$Cluster,
    #              cluster_column_slices = F,
    #              show_column_names = F,
    #              show_heatmap_legend = F)
    # if(!is.null(savePath)){
    #     png(filename = file.path(savePath, "figures/exprProgram-heatmap.png"), width = 1300, height = 800)
    #     p
    #     dev.off()
    # }

    return(p)
}




#' runCellInteraction
#'
#' @param expr A Seurat object.
#' @param cellSetName The colunm name of `expr`'s `meta.data`, used to indicate the cell set annotation.
#' @inheritParams runScAnnotation
#'
#' @return A data frame which contains the cell sets ligand-receptor pairs and their scores.
#' @export
#'
runCellInteraction <- function(expr, cellSetName = "default", species = "human", savePath = NULL){
    cat("[", paste0(Sys.time()), "] -----: cell interaction analysis\n")

    pairsLigRec <- read.table(system.file("txt", "PairsLigRec.txt", package = "scCancer"),
                              sep = "\t", header = T,stringsAsFactors = F)
    if(species == "mouse"){
        hg.mm.Genes <- read.table(system.file("txt", "hg-mm-HomologyGenes.txt", package = "scCancer"),
                                  header = T, stringsAsFactors = F)

        hg.num <- table(hg.mm.Genes$hgGenes)
        hg.mm.Genes <- subset(hg.mm.Genes, !(hgGenes %in% names(hg.num)[hg.num > 1]))
        rownames(hg.mm.Genes) <- hg.mm.Genes$hgGenes

        ju.in <- pairsLigRec$Ligand %in% hg.mm.Genes$hgGenes & pairsLigRec$Receptor %in% hg.mm.Genes$hgGenes
        pairsLigRec <- pairsLigRec[ju.in, ]

        new.pairs <- list()
        new.pairs$Ligand <- hg.mm.Genes[pairsLigRec$Ligand, "mmGenes"]
        new.pairs$Receptor <- hg.mm.Genes[pairsLigRec$Receptor, "mmGenes"]
        new.pairs$Pair.Name <- paste0(new.pairs$Ligand, "_", new.pairs$Receptor)
        pairsLigRec <- data.frame(new.pairs)
        rm(new.pairs)
    }
    ju.pairs <- (pairsLigRec$Ligand %in% rownames(expr)) & (pairsLigRec$Receptor %in% rownames(expr))
    pairsLigRec <- pairsLigRec[ju.pairs, ]
    all.genes <- unique(c(pairsLigRec$Ligand, pairsLigRec$Receptor))

    cellSets <- unique(expr@meta.data[[cellSetName]])
    set.mean <- lapply(cellSets, FUN = function(x){
        return(Matrix::rowMeans(GetAssayData(expr, slot = "counts")[all.genes, colnames(expr)[expr[[cellSetName]] == x]]))
    })
    names(set.mean) <- paste0("S.", cellSets)
    set.mean <- data.frame(set.mean)

    set.rate <- lapply(cellSets, FUN = function(x){
        return(Matrix::rowMeans(GetAssayData(expr, slot = "counts")[all.genes, colnames(expr)[expr[[cellSetName]] == x]] > 0))
    })
    names(set.rate) <- paste0("S.", cellSets)
    set.rate <- data.frame(set.rate)

    final.pairs <- data.frame("Ligand" = c(), "Receptor" = c(),
                              "Ligand.CellSet" = c(), "Receptor.CellSet" = c(), "Score" = c())
    stat.df <- data.frame("Ligand.CellSet" = c(), "Receptor.CellSet" = c(),
                          "Num" = c(), "Sum" = c())
    for(i in cellSets){
        for(j in cellSets){
            if(i != j){
                cur.pairs <- data.frame(Ligand = pairsLigRec$Ligand,
                                        Receptor = pairsLigRec$Receptor,
                                        Ligand.CellSet = paste0("S.", i),
                                        Receptor.CellSet = paste0("S.", j),
                                        Score = set.mean[pairsLigRec$Ligand, paste0("S.", i)] *
                                            set.mean[pairsLigRec$Receptor, paste0("S.", j)],
                                        Detect.rate = set.rate[pairsLigRec$Ligand, paste0("S.", i)] *
                                            set.rate[pairsLigRec$Receptor, paste0("S.", j)])
                cur.pairs <- subset(cur.pairs, Score > 0)
                final.pairs <- rbind(final.pairs, cur.pairs)

                ju.stat <- cur.pairs$Score > 1 & cur.pairs$Detect.rate > 0.25
                cur.stat <- data.frame("Ligand.CellSet" = c(i),
                                       "Receptor.CellSet" = c(j),
                                       "Num" = c(sum(ju.stat)),
                                       "Sum" = c(sum(cur.pairs$Score[ju.stat])))
                stat.df <- rbind(stat.df, cur.stat)
            }
        }
    }

    final.pairs <- final.pairs[order(final.pairs$Score, decreasing = T), ]
    rownames(final.pairs) <- NULL
    colnames(stat.df) <- c("Ligand.CellSet", "Receptor.CellSet", "Num", "Sum")
    # stat.df <- subset(stat.df, Num > 0)

    if(dim(stat.df)[1] > 0){
        if(setequal(cellSets, 1:length(cellSets))){
            stat.df$Ligand.CellSet <- factor(paste0("S.", stat.df$Ligand.CellSet), levels = paste0("S.", 1:length(cellSets)))
            stat.df$Receptor.CellSet <- factor(paste0("S.", stat.df$Receptor.CellSet), levels = paste0("S.", 1:length(cellSets)))
        }
    }

    if(!is.null(savePath)){
        write.table(final.pairs, file = file.path(savePath, "InteractionScore.txt"),
                    quote = F, sep = "\t", row.names = F)
    }
    return(list(interaction.score = final.pairs,
                stat.df = stat.df))
}



#' plotCellInteraction
#'
#' @param stat.df A data.frame of cell sets interaction result.
#' @param cell.annotation  A data.frame of cells' annotation containing the cells' cluster and type.
#'
#' @return A plot showing the result of cell interaction.
#' @export
#'
#' @importFrom gridExtra grid.arrange
#'
plotCellInteraction <- function(stat.df, cell.annotation){
    p.inter <- ggplot(stat.df, aes(x = Ligand.CellSet, y = Receptor.CellSet,
                                             size = Num, color = Sum)) +
        geom_point() +
        coord_fixed() +
        scale_color_gradientn(colors = colorRampPalette(c("#FAFBFD", "#2e68b7"))(50)) +
        theme_light() +
        theme(axis.title.x = element_text(size = 13, vjust = -0.1),
              axis.title.y = element_text(size = 13, vjust = 1.5),
              axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = 11),
              axis.ticks = element_line(color = "black"),
              panel.border = element_rect(color = "black"))
    p.inter.leg <- cowplot::get_legend(p.inter)
    p.inter <- p.inter + theme(legend.position = "none")

    sel.col = "Cell.Type"
    cell.colors <- getCellTypeColor(cell.types = unique(cell.annotation[[sel.col]]))

    bar.num <- table(cell.annotation[c(sel.col, "Cluster")])
    bar.num <- t(t(bar.num) / colSums(bar.num))
    bar.df <- melt(bar.num)
    bar.df$Cluster = factor(bar.df$Cluster)

    p.type <- ggplot(bar.df, aes(x = Cluster, y = value, fill = .data[[sel.col]])) +
        geom_bar(stat = "identity") +
        coord_equal(length(unique(bar.df$Cluster)) / 5) +
        scale_fill_manual(values = cell.colors) +
        labs(y = "Type fraction", x = NULL) +
        guides(fill = guide_legend(override.aes = list(size = 1), title = sel.col)) +
        scale_x_discrete(position = "top") +
        theme_classic() +
        theme(axis.title.x = element_text(size = 13, vjust = -0.1),
              axis.title.y = element_text(size = 13, vjust = 1.5),
              axis.text.y = element_text(size = 11),
              axis.text.x = element_blank(),
              panel.background = element_rect(colour = "black"),
              axis.line = element_line(color = "black"))
    p.type.leg <- cowplot::get_legend(p.type)
    p.type <- p.type + theme(legend.position = "none")

    p.iMain <- rbind(ggplotGrob(p.inter), ggplotGrob(p.type), size = "last")
    # p.ic <- plot_grid(p.iMain, p.inter.leg, p.type.leg, ncol = 3, rel_widths = c(6, 1, 1))
    # ggsave(filename = "../iii.png",
    #        p.ic, width = 7, height = 6.5, dpi = 300)

    p.ic <- grid.arrange(
        grobs = list(p.iMain, p.inter.leg, p.type.leg),
        widths = c(4.5, 1),
        layout_matrix = rbind(c(1, 2), c(1, 3)))
    return(p.ic)
}




#' runScAnnotation
#'
#' According to the results of 'runScStatistics', perform cell and gene quality control.
#' Using the R package Seurat to perform basic operations (normalization, log-transformation,
#' highly variable genes identification, removing unwanted variance, scaling, centering,
#' dimension reduction, clustering, and differential expression analy-sis).
#' Perform some cancer-specific analyses: cancer micro-environmental cell type classification,
#' cell malignancy estimation, cell cycle analysis, cell stemness analysis,
#' gene set signature analysis, expression programs identification, and so on.
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
#' Must be some of c("mitochondrial", "ribosome", "dissociation")(default) or NULL.
#' @param nCell.min An integer number used to filter gene. The default is 3.
#' Genes with the number of expressed cells less than this threshold will be filtered.
#' @param bgPercent.max A float number used to filter gene. The default is 1 (no filtering).
#' Genes with the background percentage larger than this threshold will be filtered.
#' @param bool.rmContamination A logical value indicating whether to remove ambient RNA contamination based on 'SoupX'.
#' @param vars.add.meta A vector indicating the variables to be added to Seurat object's meta.data.
#' The default is c("mito.percent", "ribo.percent", "diss.percent").
#' @param vars.to.regress A vector indicating the variables to regress out in R package Seurat.
#' The default is c("nCount_RNA", "mito.percent", "ribo.percent").
#' @param pc.use An integer number indicating the number of PCs to use as input features. The default is 30.
#' @param resolution A float number used in function 'FindClusters' in Seurat. The default is 0.8.
#' @param clusterStashName A character string used as the name of cluster identies. The default is "default".
#' @param show.features A list or vector for genes to be plotted in 'markerPlot'.
#' @param bool.add.features A logical value indicating whether to add default features to 'show.features' or not.
#' @param bool.runDiffExpr A logical value indicating whether to perform differential expressed analysis.
#' @param n.markers An integer indicating the number of differential expressed genes showed in the plot. The defalut is 5.
#' @param species A character string indicating what species the sample belong to.
#' Only "human"(default) or "mouse" are allowed.
#' @param genome A character string indicating the version of the reference gene annotation information.
#' This information is mainly used to infer CNV profile and estimate malignancy.
#' Only 'hg19' (defalut) or 'hg38' are allowed for "human" species, and only "mm10" is allowed for "mouse" species.
#' @param hg.mm.mix  A logical value indicating whether the sample is a mix of
#' human cells and mouse cells(such as PDX sample).
#' If TRUE, the arguments 'hg.mm.thres' and 'mix.anno' should be set to corresponding values.
#' @param bool.runDoublet A logical value indicating whether to estimate doublet scores.
#' @param doublet.method The method to estimate doublet score. The default is "bcds".
#' "cxds"(co-expression based doublet scoring) and "bcds"(binary classification based doublet scoring) are allowed.
#' These methods are from R package "scds".
#' @param bool.runCellClassify A logical value indicating whether to predict the usual cell type. The default is TRUE.
#' @param roughlabel.path Path for generated rough label, useful when bool.runCellClassify = FALSE.
#' @param bool.runCellSubtypeClassify A logical value indicating whether to predict the usual cell subtype. The default is TRUE.
#' @param celltype.list A list of cell types for subtype annotation, which depends on either rough annotation or user's input.
#' @param ct.templates A list of vectors of several cell type templates.
#' The default is NULL and the templates prepared in this package will be used.
#' @param submodel.path A folder containing all models preset for cell subtype annotation(.csv file)
#' The default is NULL and all models prepared in this package will be used.
#' @param markers.path A folder containing all marker files preset for cell subtype annotation(.txt file)
#' The default is NULL and all marker files prepared in this package will be used.
#' @param subtype.umap A logical value indicating whether to generate umap plot group by cell subtypes. The default is FALSE.
#' @param coor.names A vector indicating the names of two-dimension coordinate used in visualization.
#' @param bool.runMalignancy A logical value indicating whether to estimate malignancy.
#' @param malignancy.method The method to be used in malignant cell identification.
#' inferCNV and xgboost are allowed. Recommend "xgboost" for large dataset, "both" for small dataset.
#' @param cnv.ref.data An expression matrix of gene by cell, which is used as the normal reference during estimating malignancy.
#' The default is NULL, and an immune cells or bone marrow cells expression matrix will be used for human or mouse species, respectively.
#' @param cnv.referAdjMat An adjacent matrix for the normal reference data.
#' The larger the value, the closer the cell pair is.
#' The default is NULL, and a SNN matrix of the default ref.data will be used.
#' @param cutoff A threshold used in the CNV inference.
#' @param bool.intraTumor A logical value indicating whether to use the identified tumor clusters to
#' perform following intra-tumor heterogeneity analyses.
#' @param p.value.cutoff A threshold to decide weather the bimodality distribution of malignancy score is significant.
#' @param bool.runCellCycle A logical value indicating whether to estimate cell cycle scores.
#' @param bool.runStemness A logical value indicating whether to estimate stemness scores.
#' @param bool.runGeneSets A logical value indicating whether to estimate gene sets signature scores.
#' @param geneSets A list of gene sets to be analyzed. The default is NULL and 50 hallmark gene sets from MSigDB will be used.
#' @param geneSet.method The method to be used in calculate gene set scores. Currently, only "average" and "GSVA" are allowed.
#' @param bool.runExprProgram A logical value indicating whether to run non-negative matrix factorization (NMF) to identify expression programs.
#' @param nmf.rank 	An integer of decomposition rank used in NMF.
#' @param bool.runInteraction A logical value indicating whether to run cell set ligand-receptor interaction analysis.
#' @param genReport A logical value indicating whether to generate a .html/.md report (suggest to set TRUE).
#'
#' @return A results list with all useful objects used in the function.
#' @export
#'
#' @import Matrix knitr ggplot2 Seurat
#' @importFrom markdown markdownToHTML
#' @importFrom pheatmap pheatmap
#' @importFrom stringr str_c
#' @importFrom dplyr "%>%"
runScAnnotation <- function(dataPath,
                            statPath,
                            savePath = NULL,
                            authorName = NULL,
                            sampleName = "sc",
                            bool.filter.cell = T,
                            bool.filter.gene = T,
                            anno.filter = c("mitochondrial", "ribosome", "dissociation"),
                            nCell.min = 3, bgPercent.max = 1,
                            bool.rmContamination = F,
                            # contamination.fraction = NULL,
                            vars.add.meta = c("mito.percent", "ribo.percent", "diss.percent"),
                            vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent"),
                            npcs = 50,
                            pc.use = 30,
                            resolution = 0.8,
                            clusterStashName = "default",
                            show.features = NULL, bool.add.features = T,
                            bool.runDiffExpr = T,
                            n.markers = 5,
                            species = "human",
                            genome = "hg19",
                            hg.mm.mix = F,
                            bool.runDoublet = T,
                            doublet.method = "bcds",
                            bool.runCellClassify = T,
                            roughlabel.path = NULL,
                            bool.runCellSubtypeClassify = T,
                            subtypeClassifyMethod = "Scoring",
                            celltype.list = NULL,
                            ct.templates = NULL,
                            submodel.path = NULL,
                            markers.path = NULL,
                            unknown.cutoff = 0.3,
                            subtype.umap = FALSE,
                            coor.names = c("UMAP_1", "UMAP_2"),
                            bool.runMalignancy = T,
                            malignancy.method = "both",
                            cnv.ref.data = NULL,
                            cnv.referAdjMat = NULL,
                            cutoff = 0.1,
                            p.value.cutoff = 0.5,
                            bool.intraTumor = T,
                            bool.runCellCycle = T,
                            bool.runStemness = T,
                            bool.runGeneSets = T,
                            geneSets = NULL,
                            geneSet.method = "average",
                            bool.runExprProgram = T,
                            nmf.rank = 50,
                            bool.runInteraction = T,
                            genReport = T){

    cat("[", paste0(Sys.time()), "] START: RUN scAnnotation\n")
    results <- as.list(environment())
    checkAnnoArguments(results)

    if(is.null(savePath)){
        savePath <- sataPath
    }

    if(species == "mouse" & genome == "hg19"){
        genome <- "mm10"
    }

    if(!dir.exists(file.path(savePath, "figures/"))){
        dir.create(file.path(savePath, "figures/"), recursive = T)
    }

    suppressWarnings( dataPath <- normalizePath(dataPath, "/") )
    suppressWarnings( statPath <- normalizePath(statPath, "/") )
    suppressWarnings( savePath <- normalizePath(savePath, "/") )
    results[["dataPath"]] <- dataPath
    results[["statPath"]] <- statPath
    results[["savePath"]] <- savePath


    ## --------- filter data ---------
    t.results <- prepareSeurat(
        dataPath = dataPath,
        statPath = statPath,
        savePath = savePath,
        sampleName = sampleName,
        bool.filter.cell = bool.filter.cell,
        bool.filter.gene = bool.filter.gene,
        anno.filter = anno.filter,
        nCell.min = nCell.min,
        bgPercent.max = bgPercent.max,
        hg.mm.mix = hg.mm.mix,
        bool.rmContamination = bool.rmContamination,
        # contamination.fraction = contamination.fraction,
        vars.add.meta = vars.add.meta,
        vars.to.regress = vars.to.regress
    )
    expr <- t.results$expr
    gene.manifest <- t.results$gene.manifest
    results[["bool.rmContamination"]] = t.results$bool.rmContamination
    results[["contamination.frac"]] = t.results$contamination.frac
    rm(t.results)
    gc()

    results[["filter.thres"]] = read.table(file.path(statPath, 'cell.QC.thres.txt'),
                                           header = T, stringsAsFactors = F)

    ## --------- seurat ---------
    t.results <- runSeurat(
        expr = expr,
        savePath = savePath,
        npcs = npcs,
        pc.use = pc.use,
        resolution = resolution,
        clusterStashName = clusterStashName,
        bool.runDiffExpr = bool.runDiffExpr
    )
    expr = t.results$expr
    cell.annotation = t.results$cell.annotation
    results[["diff.expr.genes"]] = t.results$diff.expr.genes
    rm(t.results)
    gc()


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


    ## --------- doublet ---------
    if(bool.runDoublet){
        cat("[", paste0(Sys.time()), "] -----: Doublet score estimation\n")
        doubletScore <- runDoublet(expr, method = doublet.method, pc.use = pc.use)
        expr[["doublet.score"]] <- doubletScore
        cell.annotation$doublet.score <- doubletScore
        results[["doublet.plot"]] <-
            pointDRPlot(cell.annotation,
                        value = "doublet.score",
                        coor.names = coor.names,
                        colors = c("white", "#214478"),
                        discrete = F,
                        legend.position = "right",
                        legend.title = "Doublet\n score")

        cell.annotation$nUMI <- Matrix::colSums(GetAssayData(expr, slot = "counts"))
        results[["nUMI.plot"]] <-
            pointDRPlot(cell.annotation,
                        value = "nUMI",
                        coor.names = coor.names,
                        colors = c("white", "#447821"),
                        discrete = F,
                        legend.position = "right",
                        legend.title = "nUMI")

        ggsave(filename = file.path(savePath, "figures/doublet-point.png"),
               results[["doublet.plot"]], width = 5, height = 4, dpi = 300)
        ggsave(filename = file.path(savePath, "figures/nUMI-point.png"),
               results[["nUMI.plot"]], width = 5, height = 4, dpi = 300)
    }


    ## --------- cell type ---------
    if(bool.runCellClassify){
        t.results <- runCellClassify(expr, cell.annotation,
                                     coor.names = coor.names,
                                     savePath = savePath,
                                     ct.templates = ct.templates,
                                     species = species)
        expr <- t.results$expr
        cell.annotation <- t.results$cell.annotation
        results[["cellType.plot"]] <- t.results$p.results
        rm(t.results)
        expr$Cell.Type %>%
            gsub("T.cells.CD4", "T.cells", .) %>%
            gsub("T.cells.CD8", "T.cells", .) -> expr$Cell.Type
        # saveRDS(expr, file = file.path(savePath, "expr-rough.RDS"))
        saveRDS(cell.annotation, file = file.path(savePath, "rough-cell-annotation.RDS"))
        saveRDS(expr$Cell.Type, file = file.path(savePath, "rough-labels.RDS"))
    }
    else{
        expr$Cell.Type <- readRDS(roughlabel.path)
    }
    ## --------- cell subtype ---------
    if(bool.runCellSubtypeClassify){
        if(is.null(celltype.list)){
            default.list <- c("T.cells", "Myeloid.cells", "B.cells", "Fibroblast", "Endothelial")
            celltype.list <- intersect(unique(expr$Cell.Type), default.list)
        }
        folder.name <- "cellSubtypeAnno"
        if(is.null(submodel.path)){
            # if(subtypeClassifyMethod == "Scoring"){
            #     submodel.path <- system.file("csv", package = "scCancer")
            # }
            # else{
            #     submodel.path <- file.path(system.file("rds", package = "scCancer"),
            #                                "cellSubtypeTemplates-XGBoost.rds")
            #     folder.name <- "cellSubtypeAnno-XGBoost"
            # }

            # submodel.path <- system.file("csv", package = "scCancer")
            submodel.path <- file.path(system.file("rds", package = "scCancer"), "cellSubtypeTemplates.rds")
        }
        if(is.null(markers.path)){
            markers.path <- system.file("txt", package = "scCancer")
        }
        if(!dir.exists(file.path(savePath, folder.name))){
            dir.create(file.path(savePath, folder.name), recursive = T)
        }
        t.results <- runCellSubtypeClassify(expr = expr,
                                            submodel.path = submodel.path,
                                            markers.path = markers.path,
                                            savePath = file.path(savePath, folder.name),
                                            celltype.list = celltype.list,
                                            dropout.modeling = FALSE,
                                            unknown.cutoff = unknown.cutoff,
                                            umap.plot = subtype.umap)
        results[["fine.labels"]] <- t.results[["fine.labels"]]
        results[["similarity.matrix"]] <- t.results[["similarity.matrix"]]
        rm(t.results)
        saveRDS(results[["fine.labels"]], file = file.path(savePath, folder.name, "fine-labels.RDS"))
        saveRDS(results[["similarity.matrix"]], file = file.path(savePath, folder.name, "similarity-matrix.RDS"))
    }

    ## --------- malignancy ---------
    results[["cnv.anno"]] <- FALSE
    results[["xgboost.anno"]] <- FALSE
    if(bool.runMalignancy){
        # if(species != "human"){
        #     cat("- Warning in 'runScAnnotation': To perform 'runMalignancy', the argument 'species' needs to be 'human'.\n")
        #     results[["bool.runMalignancy"]] = FALSE
        # }else{
        #
        # }
        if(malignancy.method == "inferCNV" | malignancy.method == "both"){
            cat("[", paste0(Sys.time()), "] -----: malignant cells identification with inferCNV\n")
            t.results <- runMalignancy(expr = expr,
                                       gene.manifest = gene.manifest,
                                       cell.annotation = cell.annotation,
                                       savePath = savePath,
                                       cutoff = cutoff, minCell = 3,
                                       p.value.cutoff = p.value.cutoff,
                                       coor.names = coor.names,
                                       ref.data = cnv.ref.data,
                                       referAdjMat = cnv.referAdjMat,
                                       species = species,
                                       genome = genome,
                                       hg.mm.mix = hg.mm.mix)
            expr <- t.results$expr
            cell.annotation <- t.results$cell.annotation
            results[["cnvList"]] <- t.results$cnvList
            results[["referScore"]] <- t.results$referScore
            results[["ju.exist.malign"]] <- t.results$ju.exist.malign
            results[["malign.thres"]] <- t.results$malign.thres
            # results[["bimodal.pvalue"]] <- t.results$bimodal.pvalue
            results[["malign.plot.cnv"]] <- t.results$p.results
            results[["cnv.anno"]] <- TRUE
            rm(t.results)
        }
        if(malignancy.method == "xgboost" | malignancy.method == "both"){
            cat("[", paste0(Sys.time()), "] -----: malignant cells identification with XGBoost model\n")
            t.results <- predMalignantCell(expr = expr,
                                           cell.annotation = cell.annotation,
                                           malignancy.method = "xgboost",
                                           savePath = savePath,
                                           MALIGNANT.THRES = 0.5)
            cell.annotation <- t.results$cell.annotation
            results[["malign.plot.xgboost"]] <- t.results$plot
            rm(t.results)
            results[["xgboost.anno"]] <- TRUE
        }
    }


    ## --------- select tumor clusters ---------
    if(bool.intraTumor){
        tumor.clusters <- getTumorCluster(cell.annotation = cell.annotation)
        results[["tumor.clusters"]] <- tumor.clusters

        if(is.null(tumor.clusters)){
            sel.clusters <- unique(cell.annotation$Cluster)
            sel.clusters <- sel.clusters[order(sel.clusters)]
        }else{
            sel.clusters <- tumor.clusters
        }
    }else{
        sel.clusters <- unique(cell.annotation$Cluster)
        sel.clusters <- sel.clusters[order(sel.clusters)]
    }


    ## --------- cell cycle ---------
    if(bool.runCellCycle){
        CellCycle.score <- runCellCycle(expr, species = species)
        cell.annotation$CellCycle.score <- CellCycle.score
        expr[["CellCycle.score"]] <- CellCycle.score
        results[["cellCycle.plot"]] <-
            pointDRPlot(cell.annotation,
                        sel.clusters = sel.clusters,
                        value = "CellCycle.score",
                        coor.names = coor.names,
                        colors = c("white", "#009b45"),
                        discrete = F,
                        legend.position = "right",
                        legend.title = "Cell cycle\n score")
        ggsave(filename = file.path(savePath, "figures/cellCycle-point.png"),
               results[["cellCycle.plot"]], width = 5, height = 4, dpi = 300)
    }


    ## --------- stemness ---------
    if(bool.runStemness){
        stem.scores <- runStemness(X = GetAssayData(object = expr, slot = "scale.data"), species = species)
        cell.annotation[["Stemness.score"]] <- stem.scores
        expr[["Stemness.score"]] <- stem.scores

        results[["stemness.plot"]] <-
            pointDRPlot(cell.annotation,
                        sel.clusters = sel.clusters,
                        value = "Stemness.score",
                        coor.names = coor.names,
                        colors = c("white", "#ff9000"),
                        discrete = F,
                        legend.position = "right",
                        legend.title = "Stemness\n score")
        ggsave(filename = file.path(savePath, "figures/stemness-point.png"),
               results[["stemness.plot"]], width = 5, height = 4, dpi = 300)
    }


    ## --------- gene sets ----------
    if(bool.runGeneSets){
        if(is.null(geneSets)){
            geneSets <- getDefaultGeneSets(species = species)
        }
        t.scores <- runGeneSets(expr = expr, geneSets = geneSets, method = geneSet.method)
        if(!is.null(t.scores)){
            cell.annotation <- cbind(cell.annotation, t.scores)
            for(gs in colnames(t.scores)){
                expr[[gs]] <- t.scores[, gs]
            }
            bool.limit <- T
            if(geneSet.method == "GSVA"){
                bool.limit <- F
            }
            results[["geneSet.plot"]] <-
                plotGeneSet(subset(cell.annotation, Cluster %in% sel.clusters),
                            prefix = "GS__",
                            bool.limit = bool.limit,
                            savePath = savePath)
            results[["geneSetPlot.height"]] <- 0.5 + 0.11 * dim(t.scores)[2]
            rm(t.scores)
        }else{
            bool.runGeneSets = FALSE
        }
    }


    ## ---------- expression programs ----------
    if(bool.runExprProgram){
        results[["exprProgram.results"]] <- runExprProgram(expr, rank = nmf.rank,
                                                           sel.clusters = sel.clusters,
                                                           savePath = savePath)
        results[["exprProgram.plot"]] <- plotExprProgram(H = results[["exprProgram.results"]]$H,
                                                         cell.annotation,
                                                         sel.clusters = sel.clusters,
                                                         savePath = savePath)
        results[["exprProgPlot.height"]] <- 0.5 + 0.11 * dim(results[["exprProgram.results"]]$H)[1]
    }


    ## ---------- cell interaction ----------
    if(bool.runInteraction){
        t.results <- runCellInteraction(expr, cellSetName = clusterStashName,
                                        species = species, savePath = savePath)
        results[["interaction.score"]] <- t.results$interaction.score

        # results[["inter.plot"]] <-
        #     ggplot(t.results$stat.df, aes(x = Ligand.CellSet, y = Receptor.CellSet, size = Num, color = Sum)) +
        #     geom_point() +
        #     coord_fixed() +
        #     scale_color_gradientn(colors = colorRampPalette(c("#FAFBFD", "#2e68b7"))(50)) +
        #     theme_light() +
        #     theme(axis.title.x = element_text(size = 13, vjust = -0.1),
        #           axis.title.y = element_text(size = 13, vjust = 1.5),
        #           axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        #           axis.text.y = element_text(size = 11),
        #           panel.border = element_rect(color = "black"))

        results[["inter.plot"]] <- plotCellInteraction(t.results$stat.df, cell.annotation)

        ggsave(filename = file.path(savePath, "figures/interaction-score.png"),
               results[["inter.plot"]], width = 7, height = 6.5, dpi = 300)
    }

    results[["expr"]] <- expr
    results[["cell.annotation"]] <- cell.annotation

    ## -------- save ---------
    saveRDS(expr, file = file.path(savePath, "expr.RDS"))
    write.table(cell.annotation, file = file.path(savePath, "cellAnnotation.txt"),
                quote = F, sep = "\t", row.names = F)

    if(genReport){
        cat("[", paste0(Sys.time()), "] -----: report generating\n")
        if(!dir.exists(file.path(savePath, 'report-figures/'))){
            dir.create(file.path(savePath, 'report-figures/'), recursive = T)
        }
        suppressWarnings(
            knit(system.file("rmd", "main-scAnno.Rmd", package = "scCancer"),
                 file.path(savePath,'report-scAnno.md'), quiet = T)
            # knit(system.file("rmd", "main-scAnno.Rmd", package = "scCancer2"), file.path(savePath,'report-scAnno.md'), quiet = T)
        )
        markdownToHTML(file.path(savePath,'report-scAnno.md'),
                       file.path(savePath, 'report-scAnno.html'))
    }

    cat("[", paste0(Sys.time()), "] END: Finish scAnnotation\n\n")

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
genAnnoReport <- function(results, savePath){
    cat("[", paste0(Sys.time()), "] -----: report generating\n")

    if(!dir.exists(savePath)){
        dir.create(savePath, recursive = T)
    }

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

