annoSpecies <- function(expr.data,
                        cell.manifest,
                        gene.manifest,
                        hg.mm.thres = 0.9,
                        mix.anno = c("human" = "hg19", "mouse" = "mm10")){
    gene.manifest$Species <- substr(gene.manifest$EnsemblID, 1, 4)

    nUMI <- Matrix::colSums(expr.data)

    hg.gene <- subset(gene.manifest, Species == mix.anno["human"])$Symbol
    umiHgFrac = Matrix::colSums(expr.data[hg.gene,]) / nUMI

    mm.gene <- subset(gene.manifest, Species == mix.anno["mouse"])$Symbol
    umiMmFrac = Matrix::colSums(expr.data[mm.gene,]) / nUMI

    cell.species <- rep("middle", length(umiHgFrac))
    cell.species[umiHgFrac >= hg.mm.thres] <- "human"
    cell.species[umiMmFrac >= hg.mm.thres] <- "mouse"

    cell.manifest$umiHgFrac = umiHgFrac
    cell.manifest$umiMmFrac = umiMmFrac
    cell.manifest$cell.species = cell.species

    return(list(cell.manifest = cell.manifest,
                gene.manifest = gene.manifest))
}



getSingleSpecies <- function(expr.data,
                             cell.manifest,
                             gene.manifest,
                             species = "human",
                             hg.mm.thres = 0.9,
                             mix.anno = c("human" = "hg19", "mouse" = "mm10")){
    anno.res <- annoSpecies(expr.data, cell.manifest, gene.manifest,
                            hg.mm.thres, mix.anno = mix.anno)
    cell.manifest <- subset(anno.res$cell.manifest, cell.species == species)
    gene.manifest <- subset(anno.res$gene.manifest, Species == mix.anno[species])

    if(dim(gene.manifest)[1] == 0){
        stop("Error in 'getSingleSpecies': no ", mix.anno[species], " genes data found.\n")
    }
    if(dim(cell.manifest)[1] == 0){
        stop("Error in 'getSingleSpecies': no ", species, " cells data found.\n")
    }

    gene.manifest$Symbol <- substr(gene.manifest$Symbol, 6, 60)
    gene.manifest$EnsemblID <- substr(gene.manifest$EnsemblID, 6, 60)

    rownames(expr.data) <- substr(rownames(expr.data), 6, 60)
    expr.data <- expr.data[gene.manifest$Symbol, cell.manifest$barcodes]

    return(list(expr.data = expr.data,
                cell.manifest = cell.manifest,
                gene.manifest = gene.manifest))
}



getGeneManifest <- function(data.path, only.expr = T) {
    gene.loc <- paste0(data.path, '/features.tsv.gz')
    if(!grepl("feature_bc_matrix", data.path)) {
        gene.loc <- paste0(data.path, '/genes.tsv')
    }
    gene.manifest <- read.delim(gene.loc, header = FALSE, stringsAsFactors = FALSE)
    colnames(gene.manifest) <-
        c("EnsemblID", "Symbol", "Type")[1:dim(gene.manifest)[2]]
    # gene.manifest$Symbol <- gsub("_", "-", make.unique(gene.manifest$Symbol))

    if(only.expr & "Type" %in% colnames(gene.manifest)){
        gene.manifest <- subset(gene.manifest, Type == "Gene Expression")
    }
    return(gene.manifest)
}



addGeneAnno <- function (gene.manifest, species = "human"){
    annotation <- rep("other", dim(gene.manifest)[1])

    if(species == "human"){
        mito.genes <- grep('^MT-', gene.manifest$Symbol, value = TRUE)
        ribo.genes <- grep('^RPL|^RPS|^MRPL|^MRPS', gene.manifest$Symbol, value = TRUE)
        diss.genes <- read.table(system.file("txt", "diss-genes.txt", package = "scCancer"),
                                 header = F, stringsAsFactors = F)$V1
    }else if(species == "mouse"){
        mito.genes <- grep('^mt-', gene.manifest$Symbol, value = TRUE)
        ribo.genes <- grep('^Rpl|^Rps|^Mrpl|^Mrps', gene.manifest$Symbol, value = TRUE)
        diss.genes <- read.table(system.file("txt", "diss-genes.txt", package = "scCancer"),
                                 header = F, stringsAsFactors = F)$V1
        diss.genes <- getMouseGene(diss.genes)
    }

    annotation[gene.manifest$Symbol %in% mito.genes] <- "mitochondrial"
    annotation[gene.manifest$Symbol %in% ribo.genes] <- "ribosome"
    annotation[gene.manifest$Symbol %in% diss.genes] <- "dissociation"

    gene.manifest$Annotation <- annotation
    return(gene.manifest)
}



#' prepareData
#'
#' @param samplePath A path containing the cell ranger processed data.
#' @inheritParams runScStatistics
#'
#' @return A list of expr.data, cell.manifest, gene.manifest, raw.data, min.nUMI, cr.version and run.emptydrop
#' @export
#'
prepareData <- function(samplePath,
                        species = "human",
                        hg.mm.mix = F,
                        hg.mm.thres = 0.9,
                        mix.anno = c("human" = "hg19", "mouse" = "mm10")) {
    raw.data = T
    data.path <- get10Xpath(samplePath, raw.data = raw.data)
    if(is.null(data.path)){
        raw.data = F
        data.path <- get10Xpath(samplePath, raw.data = raw.data)
        if(is.null(data.path)){
            stop("Cannot find the raw data or filtered data.\n")
        }else{
            cat("- Warning in 'prepareData': Cannot find the raw data, and use the filtered data instead.\n")
        }
    }

    expr.data <- Read10Xdata(data.dir = data.path, only.expr = T)
    # rownames(expr.data) <- gsub("_", "-", rownames(expr.data))
    cr.version <- getCRversion(data.path)

    run.emptydrop <- F
    if(raw.data){
        filter.path <- get10Xpath(samplePath, raw.data = F)
        if(!is.null(filter.path) & cr.version == "Cell Ranger (version >= 3)"){
            filtered.cell <- getBarcodes(filter.path)
        }else{
            if(requireNamespace("DropletUtils", quietly = TRUE)){
                run.emptydrop <- T
                e.out <- emptyDrops(expr.data)
                is.cell <- e.out$FDR <= 0.01
                filtered.cell <- colnames(expr.data)[which(is.cell)]
            }else{
                if(!is.null(filter.path)){
                    cat("- Warning in 'prepareData': Package 'DropletUtils' isn't installed and the pipeline will use the supplied CR2 filtered data instead.\n")
                    filtered.cell <- getBarcodes(filter.path)
                }else{
                    stop("Please install 'DropletUtils' package or provide filtered data.\n")
                }
            }
        }
    }

    ## gene.manifest
    gene.manifest <- getGeneManifest(data.path, only.expr = T)
    gene.manifest$Symbol <- make.unique(gene.manifest$Symbol)

    ## cell.manifest
    nGene = Matrix::colSums(expr.data > 0)
    nUMI = Matrix::colSums(expr.data)
    droplet.type <- rep("cell", dim(expr.data)[2])
    names(droplet.type) <- colnames(expr.data)

    if(raw.data){
        min.nUMI <- min(nUMI[colnames(expr.data) %in% filtered.cell])
        droplet.type[!(names(droplet.type) %in% filtered.cell)] <- "empty"
    }else{
        min.nUMI <- min(nUMI)
    }
    cell.manifest = data.frame(
        barcodes = colnames(expr.data),
        droplet.type = droplet.type,
        nUMI = nUMI,
        nGene = nGene
    )
    rm(nGene, nUMI, droplet.type)

    cell.manifest = cell.manifest[cell.manifest$nUMI > 0, ]
    expr.data <- expr.data[, cell.manifest$barcodes]

    ## hg.mm.mix : multi-species
    if(hg.mm.mix){
        anno.res <- annoSpecies(expr.data, cell.manifest, gene.manifest,
                                hg.mm.thres, mix.anno = mix.anno)
        cell.manifest <- subset(anno.res$cell.manifest, cell.species == species)
        gene.manifest <- subset(anno.res$gene.manifest, Species == mix.anno[species])

        if(dim(gene.manifest)[1] == 0){
            stop("Error in 'getSingleSpecies': no ", mix.anno[species], " genes data found.\n")
        }
        if(dim(cell.manifest)[1] == 0){
            stop("Error in 'getSingleSpecies': no ", species, " cells data found.\n")
        }

        gene.manifest$Symbol <- substr(gene.manifest$Symbol, 6, 60)
        gene.manifest$EnsemblID <- substr(gene.manifest$EnsemblID, 6, 60)
        rownames(expr.data) <- substr(rownames(expr.data), 6, 60)
        rm(anno.res)

        expr.data <- expr.data[gene.manifest$Symbol, cell.manifest$barcodes]
    }

    gene.manifest <- addGeneAnno(gene.manifest, species = species)
    geneType <- c("mitochondrial", "ribosome", "dissociation")
    for(tt in geneType) {
        cur.genes <- subset(gene.manifest, Annotation == tt, select = c("Symbol"))$Symbol
        cur.percent <- Matrix::colSums(expr.data[cur.genes,]) / Matrix::colSums(expr.data)
        cur.percent[is.na(cur.percent)] = 0
        cur.name = paste0(substr(tt, 1, 4), ".percent")
        cell.manifest[[cur.name]] = cur.percent
    }

    gene.manifest$Symbol <- gsub("_", "-", gene.manifest$Symbol)

    results <- list(expr.data = expr.data,
                    cell.manifest = cell.manifest,
                    gene.manifest = gene.manifest,
                    raw.data = raw.data,
                    min.nUMI = min.nUMI,
                    cr.version = cr.version,
                    run.emptydrop = run.emptydrop)

    return(results)
}




cellsPlot <- function(cell.manifest, plot.type = "histogram") {
    cur.color <- c("cell" = '#825f87', "empty" = '#aeacac')

    if(plot.type == "histogram"){
        p <- ggplot() +
            geom_histogram(
                data = subset(cell.manifest, droplet.type == "empty"),
                mapping = aes(x = nUMI, fill = "empty"),
                bins = 200, alpha = 0.7) +
            geom_histogram(
                data = subset(cell.manifest, droplet.type == "cell"),
                mapping = aes(x = nUMI, fill = "cell"),
                bins = 200, alpha = 0.7) +
            geom_freqpoly(data = cell.manifest, mapping = aes(x = nUMI), bins = 200) +
            labs(y = "Droplet number") +
            scale_fill_manual(
                name = 'Droplet type',
                guide = 'legend',
                # labels = c("cell", "empty"),
                values = cur.color) +
            scale_x_continuous(trans = 'log10') + scale_y_continuous(trans = 'log10') +
            guides(fill = guide_legend(override.aes = list(size = 1, alpha = 0.7))) +
            ggplot_config(base.size = 6) +
            theme(legend.position = "bottom")

    }else if(plot.type == "rankplot"){
        tmp.df <- cell.manifest[order(cell.manifest["nUMI"], decreasing = TRUE),]
        xi <- 1:dim(tmp.df)[1]
        tmp.df$xi = xi
        p <- ggplot(tmp.df, aes(x = xi, y = nUMI, color = droplet.type)) +
            geom_point(cex = 1, alpha = 0.5) +
            labs(x = "Droplets", y = "nUMI") +
            scale_color_manual(values = cur.color, name = "Droplet type") +
            scale_x_continuous(trans = 'log10') + scale_y_continuous(trans = 'log10') +
            guides(color = guide_legend(override.aes = list(size = 5))) +
            ggplot_config(base.size = 6) +
            theme(legend.position = "bottom")

    }else{
        stop("Error in 'cellPlot': invalid 'plot.type' argument.\n")
    }

    return(p)
}


calcThres <- function(cell.manifest, values = c("nUMI", "nGene")){
    high.thresholds <- list()
    for(value in values){
        outliers <- getOutliers(cell.manifest[[value]])
        high.thresholds[[value]] <- round(min(outliers), 3)
    }
    return(high.thresholds)
}


histPlot <- function(cell.manifest, value, xlines = c()){
    p <- ggplot(cell.manifest, aes(x = .data[[value]])) +
        geom_histogram(bins = 200, position = "stack", fill = "#a788ab") +
        labs(x = value, y = "Cell number") +
        ggplot_config(base.size = 6)

    for(x in xlines){
        p <- p + geom_vline(xintercept = x, colour = "red", linetype = "dashed")
    }
    return(p)
}


marginPlot <- function(cell.manifest, value, color = "#a788ab", xlines = c(), ylines = c()){
    p <- ggplot(cell.manifest, aes(y = .data[[value]], x = nUMI)) +
        geom_point(cex = 0.8, alpha = 0.1, color = color) +
        labs(y = value) +
        ggplot_config(base.size = 7)
    for(x in xlines){
        p <- p + geom_vline(xintercept = x, colour = "red", linetype = "dashed")
    }
    for(y in ylines){
        p <- p + geom_hline(yintercept = y, colour = "red", linetype = "dashed")
    }
    p <- ggMarginal(p, bins = 200, type = "histogram", fill = color, color = color)
    return(p)
}


getBgPercent <- function(cell.manifest.all, expr.data, bg.low = 1, bg.up = 10){
    cell.manifest.all$barcodes <- rownames(cell.manifest.all)
    bg.bars <- subset(cell.manifest.all, nUMI >= bg.low & nUMI <= bg.up)$barcodes

    if(length(bg.bars) == 0){
        result <- NULL
    }else{
        bg.sum <- Matrix::rowSums(expr.data[, bg.bars])
        bg.percent <- bg.sum / sum(bg.sum)
        result <- data.frame(row.names = rownames(expr.data),
                             est = bg.percent,
                             counts = bg.sum)
    }
    return(result)
}


getNcell <- function(cell.manifest, expr.data){
    tmp.mat <- expr.data[, cell.manifest$barcodes]
    nCell <- Matrix::rowSums(tmp.mat > 0)
    return(nCell)
}


getDetectRate <- function(nCell, n){
    detect.rate <- nCell / n
    return(detect.rate)
}


getExprProp <- function(cell.manifest, expr.data){
    tmp.mat <- expr.data[, cell.manifest$barcodes]
    expr.frac <- t(t(tmp.mat) / Matrix::colSums(tmp.mat))
    return(expr.frac)
}


getLowFrac <- function(expr.frac, bg.percent){
    expr.frac <- as(expr.frac, "dgTMatrix")
    expr.frac@x <- expr.frac@x - bg.percent[expr.frac@i + 1]
    low.frac <- Matrix::rowSums(expr.frac < 0) / dim(expr.frac)[2]
    return(low.frac)
}


getPropMedian <-function(expr.frac, nCell){
    prop.median <- rep(0, length(nCell))
    names(prop.median) <- names(nCell)
    tmp.sel.genes <- names(nCell)[nCell >= dim(expr.frac)[2] / 2]
    prop.median[tmp.sel.genes] <- apply(expr.frac[tmp.sel.genes, ], 1, median)
    return(prop.median)
}


updateGeneM <- function(gene.manifest,
                        bg.percent = NULL,
                        nCell = NULL,
                        detect.rate = NULL,
                        prop.median = NULL,
                        low.frac = NULL){

    addGeneM <- function(gene.manifest, values, name){
        if(!is.null(values)){
            gene.manifest[[name]] <- values
        }
        return(gene.manifest)
    }

    gene.manifest <- addGeneM(gene.manifest, nCell, "nCell")
    gene.manifest <- addGeneM(gene.manifest, prop.median, "prop.median")
    gene.manifest <- addGeneM(gene.manifest, bg.percent, "bg.percent")
    gene.manifest <- addGeneM(gene.manifest, detect.rate, "detect.rate")
    gene.manifest <- addGeneM(gene.manifest, low.frac, "low.frac")

    return(gene.manifest)
}




genePropPlot <- function(gene.manifest, expr.frac){
    show.num = 100
    if("bg.percent" %in% colnames(gene.manifest)){
        gene.manifest <- gene.manifest[order(gene.manifest$bg.percent, decreasing = T), ]
    }else{
        gene.manifest <- gene.manifest[order(gene.manifest$prop.median, decreasing = T), ]
    }
    gene.show <- head(gene.manifest$Symbol, show.num)

    rownames(gene.manifest) <- gene.manifest$Symbol
    gene.manifest <- gene.manifest[gene.show, ]

    rate.df.plot <- melt(as.data.frame(as.matrix(t(expr.frac[gene.show, ]))), id.vars = NULL)
    rate.df.plot$Annotation <- factor(gene.manifest[as.character(rate.df.plot$variable), ]$Annotation,
                                      levels = c("mitochondrial", "ribosome", "dissociation", "other"), ordered = T)
    rate.df.plot$variable <- factor(rate.df.plot$variable,
                                    levels = rev(gene.show), ordered = TRUE)

    sub.ix <- c(rep("1-50", 50), rep("51-100", 50))
    names(sub.ix) <- gene.show
    rate.df.plot$subPlot <- factor(sub.ix[as.character(rate.df.plot$variable)],
                                   levels = c("1-50", "51-100"), ordered = T)

    if("bg.percent" %in% colnames(gene.manifest)){
        bg.df <- data.frame(ix = c(1:50, 1:50) + 0.1,
                            frac = rev(gene.manifest[gene.show, ]$bg.percent),
                            type = "Background proportion")
        bg.df$subPlot <- factor(sub.ix[as.character(rev(gene.show))],
                                levels = c("1-50", "51-100"), ordered = T)
    }else{
        bg.df <- NULL
    }

    gene.colors <- c(other = "#d56f6d", ribosome = "#f29721", dissociation = "#65a55d", mitochondrial = "#3778bf")
    p <- ggplot() + coord_flip() +
        scale_color_manual(values = "red", name = "") +
        scale_fill_manual(values = gene.colors, name = "Gene:") +
        scale_y_continuous(breaks = c(0.0001, 0.001, 0.01, 0.1, 1),
                           labels = c("0.0001", "0.001", "0.01", "0.1", "1"),
                           trans = 'log10') +
        theme_bw() +
        theme(axis.title.y = element_blank(),
              legend.position = "top",
              legend.title = element_text(face="bold"))

    p1 <- p + geom_boxplot(data = subset(rate.df.plot, subPlot == "1-50"),
                           mapping = aes(x = variable, y = value, fill = Annotation),
                           outlier.size = 0.1, alpha = 0.8) +
        xlab("") + ylab(paste0("Gene proportion of total UMI (1-50)"))

    p2 <- p + geom_boxplot(data = subset(rate.df.plot, subPlot == "51-100"),
                           mapping = aes(x = variable, y = value, fill = Annotation),
                           outlier.size = 0.1, alpha = 0.8) +
        xlab("") + ylab(paste0("Gene proportion of total UMI (51-100)"))

    all.p <- p + geom_boxplot(data = rate.df.plot,
                              mapping = aes(x = variable, y = value, fill = Annotation),
                              outlier.size = 0.1, alpha = 0.8)

    if(!is.null(bg.df)){
        p1 <- p1 + geom_point(data = subset(bg.df, subPlot == "1-50"),
                              aes(x = ix, y = frac, color = type), shape = "*", size = 5)
        p2 <- p2 + geom_point(data = subset(bg.df, subPlot == "51-100"),
                              aes(x = ix, y = frac, color = type), shape = "*", size = 5)
        all.p <- all.p + geom_point(data = bg.df, aes(x = ix, y = frac, color = type), shape = "*", size = 5)
    }

    p <- grid_arrange_shared_legend(p1, p2, all.p = all.p, ncol = 2, nrow = 1)

    # return(list(rdf = rate.df.plot, bdf = bg.df, p1 = p1, p2 = p2, all.p = all.p, p = p))
    return(p)
}


bgCellScatter <- function(gene.manifest){
    tmp.gene.mani <- subset(gene.manifest, nCell > 10)
    tmp.gene.mani$Annotation <- factor(tmp.gene.mani$Annotation, ordered = T,
                                       levels = c("mitochondrial", "ribosome", "dissociation", "other"))

    gene.colors <- c(other = "#d56f6d", ribosome = "#f29721", dissociation = "#65a55d", mitochondrial = "#3778bf")
    p <- ggplot() +
        geom_point(tmp.gene.mani,
                   mapping = aes(x = prop.median, y = bg.percent, color = Annotation),
                   alpha = 0.3) +
        scale_color_manual(values = gene.colors, name = "Gene:") +
        coord_fixed() +
        scale_x_continuous(trans = 'log10',
                           limits = c(0.00001, 0.1),
                           breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1),
                           labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1")) +
        scale_y_continuous(trans = 'log10',
                           limits = c(0.00001, 0.1),
                           breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1),
                           labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1")) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
        ggplot_config(base.size = 6) +
        theme(legend.justification = c(0,1), legend.position = c(0.01,1)) +
        guides(color = guide_legend(override.aes = list(size = 2),
                                    keywidth = 0.1,
                                    keyheight = 0.15,
                                    default.unit = "inch")) +
        labs(x = "Median of gene proportion in cells", y = "Gene proportion in background")

    return(p)
}


bgDetScatter <- function(gene.manifest){
    tmp.gene.mani <- subset(gene.manifest, nCell > 10)
    tmp.gene.mani$Annotation <- factor(tmp.gene.mani$Annotation, ordered = T,
                                       levels = c("mitochondrial", "ribosome", "dissociation", "other"))

    gene.colors <- c(other = "#d56f6d", ribosome = "#f29721", dissociation = "#65a55d", mitochondrial = "#3778bf")
    p <- ggplot() +
        geom_point(tmp.gene.mani,
                   mapping = aes(x = detect.rate, y = bg.percent, color = Annotation),
                   alpha = 0.3) +
        scale_color_manual(values = gene.colors, name = "Gene") +
        coord_fixed(1/4) +
        scale_y_continuous(trans = 'log10',
                           limits = c(0.00001, 0.1),
                           breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1),
                           labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1")) +
        ggplot_config(base.size = 6) +
        theme(legend.justification = c(0,1), legend.position = c(0.01,1)) +
        guides(color = guide_legend(override.aes = list(size = 2),
                                    keywidth = 0.1,
                                    keyheight = 0.15,
                                    default.unit = "inch")) +
        labs(x = "Gene detection rate in cells", y = "Gene proportion in background")
    return(p)
}


#' runScStatistics
#'
#' perform cell QC, gene QC, visualization and give suggested thresholds.
#'
#' @param dataPath A path containing the cell ranger processed data.
#' Under this path, folders 'filtered_feature_bc_matrix' and 'raw_feature_bc_matrix' exist generally.
#' @param savePath A path to save the results files(suggest to create a foler named by sample name).
#' @param authorName A character string for authors name and will be shown in the report.
#' @param sampleName A character string giving a label for this sample.
#' @param species A character string indicating what species the sample belong to.
#' Must be one of "human"(default) and "mouse".
#' @param hg.mm.mix A logical value indicating whether the sample is a mix of
#' human cells and mouse cells(such as PDX sample).
#' If TRUE, the arguments 'hg.mm.thres' and 'mix.anno' should be set to corresponding values.
#' @param hg.mm.thres A float-point threshold within [0.5, 1] to identify human and mouse cells.
#' Cells with UMI percentage of single species larger than the threshold are labeled human or mouse cells.
#' The default is 0.6.
#' @param mix.anno A vector to indicate the prefix of genes from different species.
#' The default is c("human" = "hg19", "mouse" = "mm10").
#' @param bg.spec.genes A list of backgroud specific genes, which are used to remove ambient genes' influence.
#' @param bool.filterhighUMI A logical value indicating whether to filter cells with extremely high nUMI
#' @param bool.runSoupx A logical value indicating whether to estimate contamination fraction by SoupX.
#' @param genReport A logical value indicating whether to generate a .html/.md report (suggest to set TRUE).
#'
#' @return A results list with all useful objects used in the function.
#' @export
#'
#' @import Matrix knitr ggplot2 SoupX
#' @importFrom markdown markdownToHTML
#' @importFrom ggExtra ggMarginal
#' @importFrom reshape2 melt
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid unit.c grid.newpage grid.draw
#' @importFrom dplyr "%>%" top_n group_by
#' @importFrom cowplot get_legend plot_grid
#' @importFrom grDevices boxplot.stats colorRampPalette
#' @importFrom stats cor density filter median quantile sd
#' @importFrom utils read.delim read.table write.csv write.table
#'
runScStatistics <- function(dataPath, savePath,
                            authorName = NULL,
                            sampleName = "sc",
                            species = "human",
                            hg.mm.mix = F,
                            hg.mm.thres = 0.6,
                            mix.anno = c("human" = "hg19", "mouse" = "mm10"),
                            bg.spec.genes = NULL,
                            bool.filterhighUMI = T,
                            bool.runSoupx = F,
                            genReport = T){

    # message("[", Sys.time(), "] START: RUN scStatistics")
    cat("[", paste0(Sys.time()), "] START: RUN scStatistics\n")
    # results <- as.list(environment())
    checkStatArguments(as.list(environment()))

    # if(bool.runSoupx){
    #     bool.runSoupx <- F
    # }

    if(!dir.exists(file.path(savePath, "figures/"))){
        dir.create(file.path(savePath, "figures/"), recursive = T)
    }

    suppressWarnings( dataPath <- normalizePath(dataPath, "/") )
    suppressWarnings( savePath <- normalizePath(savePath, "/") )

    cat("[", paste0(Sys.time()), "] -----: data preparation\n")
    all <- prepareData(samplePath = dataPath,
                       species = species,
                       hg.mm.mix = hg.mm.mix,
                       hg.mm.thres = hg.mm.thres,
                       mix.anno = mix.anno)

    expr.data <- all$expr.data
    cell.manifest <- all$cell.manifest
    gene.manifest <- all$gene.manifest
    raw.data <- all$raw.data
    min.nUMI <- all$min.nUMI
    cr.version <- all$cr.version
    run.emptydrop <- all$run.emptydrop
    rm(all)

    cat("[", paste0(Sys.time()), "] -----: cell calling\n")
    if(raw.data){
        p.cells.1 <- cellsPlot(cell.manifest, plot.type = "histogram")
        p.cells.2 <- cellsPlot(cell.manifest, plot.type = "rankplot")
        suppressWarnings(
            ggsave(filename = file.path(savePath, "figures/cells-distr-hist.png"),
                   p.cells.1, dpi = 300, height = 3, width = 4)
        )
        suppressWarnings(
            ggsave(filename = file.path(savePath, "figures/cells-distr-rank.png"),
               p.cells.2, dpi = 300, height = 3, width = 4)
        )
    }else{
        p.cells.1 <- NULL
        p.cells.2 <- NULL
    }

    cat("[", paste0(Sys.time()), "] -----: nUMI & nGene distribution plot\n")
    cell.manifest.all <- cell.manifest
    cell.manifest <- subset(cell.manifest, droplet.type == "cell")
    cell.threshold <- calcThres(cell.manifest,
                                values = c("nUMI", "nGene", "mito.percent", "ribo.percent", "diss.percent"))
    if (bool.filterhighUMI == F){
        cell.threshold$nUMI <- max(cell.manifest[["nUMI"]])
    }
    p.nUMI <- histPlot(cell.manifest, value = "nUMI", xlines = c(cell.threshold$nUMI))
    p.nGene <- histPlot(cell.manifest, value = "nGene", xlines = c(200, cell.threshold$nGene))
    ggsave(filename = file.path(savePath, "figures/nUMI-distr.png"),
           p.nUMI, dpi = 300, height = 2.5, width = 4)
    ggsave(filename = file.path(savePath, "figures/nGene-distr.png"),
           p.nGene, dpi = 300, height = 2.5, width = 4)


    cat("[", paste0(Sys.time()), "] -----: mito & ribo & diss distribution plot\n")
    p.mito <- marginPlot(cell.manifest, value = "mito.percent", color = "#6d9fd5",
                         xlines = c(cell.threshold$nUMI), ylines = c(cell.threshold$mito.percent))
    p.ribo <- marginPlot(cell.manifest, value = "ribo.percent", color = "#f6b969",
                         xlines = c(cell.threshold$nUMI), ylines = c(cell.threshold$ribo.percent))
    p.diss <- marginPlot(cell.manifest, value = "diss.percent", color = "#94c08e",
                         xlines = c(cell.threshold$nUMI), ylines = c(cell.threshold$diss.percent))
    ggsave(filename = file.path(savePath, "figures/mito-distr.png"),
           p.mito, dpi = 300, height = 4, width = 4)
    ggsave(filename = file.path(savePath, "figures/ribo-distr.png"),
           p.ribo, dpi = 300, height = 4, width = 4)
    ggsave(filename = file.path(savePath, "figures/diss-distr.png"),
           p.diss, dpi = 300, height = 4, width = 4)


    cat("[", paste0(Sys.time()), "] -----: gene statistics\n")
    bg.result <- getBgPercent(cell.manifest.all, expr.data, bg.low = 1, bg.up = 10)
    bg.percent <- bg.result$est
    nCell <- getNcell(cell.manifest, expr.data)
    detect.rate <- getDetectRate(nCell, n = dim(cell.manifest)[1])
    expr.frac <- getExprProp(cell.manifest, expr.data)
    # low.frac <- getLowFrac(expr.frac, bg.percent)
    prop.median <- getPropMedian(expr.frac, nCell)

    gene.manifest <- updateGeneM(gene.manifest,
                                 bg.percent = bg.percent,
                                 nCell = nCell,
                                 detect.rate = detect.rate,
                                 prop.median = prop.median,
                                 low.frac = NULL)

    cat("[", paste0(Sys.time()), "] -----: gene proportion plot\n")
    suppressWarnings(
        p.geneProp <- genePropPlot(gene.manifest, expr.frac)
    )
    suppressWarnings(
        ggsave(filename = file.path(savePath, "figures/geneProp.png"),
               p.geneProp, dpi = 300, height = 8, width = 8)
    )

    if(!is.null(bg.percent) && raw.data){
        p.bg.cell <- bgCellScatter(gene.manifest)
        p.bg.detect <- bgDetScatter(gene.manifest)
        suppressWarnings(
            ggsave(filename = file.path(savePath, "figures/bg-cell-scatter.png"),
                   p.bg.cell, dpi = 300, height = 4, width = 4)
        )
        suppressWarnings(
            ggsave(filename = file.path(savePath, "figures/bg-detect-scatter.png"),
                   p.bg.detect, dpi = 300, height = 4, width = 4)
        )
    }else{
        p.bg.cell <- NULL
        p.bg.detect <- NULL
    }


    if(bool.runSoupx){
        cat("[", paste0(Sys.time()), "] -----: ambient genes (SoupX)\n")
        suppressWarnings( cls.path <- file.path(dataPath, "analysis/clustering/graphclust/clusters.csv") )
        bool.cls.info <- file.exists(cls.path)
        if(is.null(bg.percent) | !raw.data){
            bool.runSoupx <- F
            cat("- Warning in 'SoupX': Cannot identify background distribution. So skip this step.\n")
        }
        if(!bool.cls.info){
            bool.runSoupx <- F
            cat("- Warning in 'SoupX': Cannot find clustering information file produced by cell ranger (`", normalizePath(cls.path),
                "`). So skip this step.\n", sep = "")
        }
    }
    if(bool.runSoupx){
        soupx.obj = load10X(dataPath, verbose = F)
        soupx.obj = autoEstCont(soupx.obj)
        contamination.frac = soupx.obj$fit$rhoEst
        saveRDS(soupx.obj, file.path(savePath, "soupx-object.RDS"))
    }else{
        # bg.spec.genes = NULL
        soupx.obj = NULL
        contamination.frac = NULL
        # p.bg.genes = NULL
    }


    cat("[", paste0(Sys.time()), "] -----: results saving\n")
    filter.thres <- list(
        Index = c("nUMI", "nGene", "mito.percent", "ribo.percent", "diss.percent"),
        Low.threshold = c(0, 200, -Inf, -Inf, -Inf),
        High.threshold = c(cell.threshold$nUMI,
                           cell.threshold$nGene,
                           cell.threshold$mito.percent,
                           cell.threshold$ribo.percent,
                           cell.threshold$diss.percent)
    )
    filter.thres <- data.frame(filter.thres)

    write.table(gene.manifest, file = file.path(savePath, "geneManifest.txt"),
                quote = F, sep = "\t", row.names = F)
    write.table(cell.manifest.all, file = file.path(savePath, "cellManifest-all.txt"),
                quote = F, sep = "\t", row.names = F)
    write.table(filter.thres, file = file.path(savePath, "cell.QC.thres.txt"),
                quote = F, sep = "\t", row.names = F)


    nList <- c(dim(cell.manifest.all)[1], dim(cell.manifest)[1])
    tmp <- cell.manifest[cell.manifest$nUMI < cell.threshold$nUMI, ]
    nList <- c(nList, dim(tmp)[1])
    tmp <- tmp[tmp$nGene >= 200, ]
    nList <- c(nList, dim(tmp)[1])
    tmp <- tmp[tmp$nGene < cell.threshold$nGene, ]
    nList <- c(nList, dim(tmp)[1])
    tmp <- tmp[tmp$mito.percent < cell.threshold$mito.percent, ]
    nList <- c(nList, dim(tmp)[1])
    tmp <- tmp[tmp$ribo.percent < cell.threshold$ribo.percent, ]
    nList <- c(nList, dim(tmp)[1])
    tmp <- tmp[tmp$diss.percent < cell.threshold$diss.percent, ]
    nList <- c(nList, dim(tmp)[1])
    rm(tmp)

    results <- list(dataPath = dataPath,
                    savePath = savePath,
                    authorName = authorName,
                    sampleName = sampleName,
                    species = species,
                    hg.mm.mix = hg.mm.mix,
                    mix.anno = mix.anno,
                    hg.mm.thres = hg.mm.thres,
                    expr.data = expr.data,
                    cell.manifest = cell.manifest,
                    gene.manifest = gene.manifest,
                    cell.threshold = cell.threshold,
                    filter.thres = filter.thres,
                    p.cells.1 = p.cells.1,
                    p.cells.2 = p.cells.2,
                    p.nUMI = p.nUMI,
                    p.nGene = p.nGene,
                    p.mito = p.mito,
                    p.ribo = p.ribo,
                    p.diss = p.diss,
                    p.geneProp = p.geneProp,
                    p.bg.cell = p.bg.cell,
                    p.bg.detect = p.bg.detect,
                    min.nUMI = min.nUMI,
                    cr.version = cr.version,
                    raw.data = raw.data,
                    run.emptydrop = run.emptydrop,
                    bool.runSoupx = bool.runSoupx,
                    # bg.spec.genes = bg.spec.genes,
                    soupx.obj = soupx.obj,
                    contamination.frac = contamination.frac,
                    nList = nList)

    ## generate report
    if(genReport){
        cat("[", paste0(Sys.time()), "] -----: report generating\n")
        # results$cell.manifest <- subset(cell.manifest.all, droplet.type == "cell")

        if(!dir.exists(file.path(savePath, 'report-figures/'))){
            dir.create(file.path(savePath, 'report-figures/'), recursive = T)
        }
        suppressWarnings(
            knit(system.file("rmd", "main-scStat.Rmd", package = "scCancer"),
                 file.path(savePath,'report-scStat.md'), quiet = T)
        )
        markdownToHTML(file.path(savePath,'report-scStat.md'),
                       file.path(savePath, 'report-scStat.html'))

        # results$cell.manifest <- cell.manifest.all
    }

    saveRDS(results, file = file.path(savePath, "scStatistics-results.RDS"))
    cat("[", paste0(Sys.time()), "] END: Finish scStatistics\n\n")
    return(results)
}



#' genStatReport
#'
#' @param results A list generated by 'runScStatistics'
#' @inheritParams runScStatistics
#'
#' @return NULL
#' @export
#'
genStatReport <- function(results, savePath){
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
        knit(system.file("rmd", "main-scStat.Rmd", package = "scCancer"),
             file.path(savePath,'report-scStat.md'), quiet = T)
    )
    markdownToHTML(file.path(savePath,'report-scStat.md'),
                   file.path(savePath, 'report-scStat.html'))
}

