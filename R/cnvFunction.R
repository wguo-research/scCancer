prepareCNV <- function(expr.data,
                       gene.manifest,
                       cell.annotation,
                       ref.data = NULL,
                       species = "human",
                       genome = "hg19",
                       hg.mm.mix = F){
    ## gene.chr
    if(species == "human"){
        if(genome == "hg38"){
            gene.chr <- read.table(system.file("txt", "gene-chr-hg38.txt", package = "scCancer"),
                                   col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
                                   stringsAsFactors = F)
        }else if(genome == "hg19"){
            gene.chr <- read.table(system.file("txt", "gene-chr-hg19.txt", package = "scCancer"),
                                   col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
                                   stringsAsFactors = F)
        }else{
            stop("Error in 'runInferCNV': ", genome, " is not allowed for 'genome'.\n")
        }
    }else if(species == "mouse"){
        if(genome == "mm10"){
            gene.chr <- read.table(system.file("txt", "gene-chr-mm10.txt", package = "scCancer"),
                                   col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
                                   stringsAsFactors = F)
        }else{
            stop("Error in 'runInferCNV': ", genome, " is not allowed for 'genome'.\n")
        }
    }else{
        stop("Error in 'runInferCNV': ", species, " is not allowed for 'species'.\n")
    }

    ## reference.data
    if(is.null(ref.data)){
        if(species == "human"){
            ref.data <- readRDS(system.file("rds", "cnvRef_Data-HM.RDS", package = "scCancer"))
        }else if(species == "mouse"){
            ref.data <- readRDS(system.file("rds", "cnvRef_Data-boneMarrow-MS.RDS", package = "scCancer"))
        }
    }
    ref.anno <- data.frame(cellName = colnames(ref.data),
                           cellAnno = "Reference",
                           stringsAsFactors = F)

    ## combine data
    com.genes <- intersect(rownames(expr.data), rownames(ref.data))
    ref.data <- ref.data[com.genes, ]
    expr.data <- expr.data[com.genes, ]
    expr.data <- cbind(as.matrix(expr.data), ref.data)
    rownames(expr.data) <- gene.manifest[rownames(expr.data), ]$EnsemblID

    ## combine cell.anno
    cell.anno <- data.frame(cellName = cell.annotation$barcodes,
                            cellAnno = cell.annotation$Cell.Type,
                            stringsAsFactors = F)
    cell.anno <- rbind(cell.anno, ref.anno)
    rownames(cell.anno) <- cell.anno$cellName

    ## common genes between expr.data and gene.chr
    com.genes <- intersect(rownames(expr.data), gene.chr$EnsemblID)
    gene.chr <- subset(gene.chr, EnsemblID %in% com.genes)

    # gene.chr <- gene.chr[with(gene.chr, order(CHR, C_START, C_STOP)), ]

    expr.data <- expr.data[gene.chr$EnsemblID, ]
    rownames(gene.chr) <- gene.chr$EnsemblID

    return(list(expr.data = expr.data,
                gene.chr = gene.chr,
                cell.anno = cell.anno))
}



rmGeneForCNV <- function(cnvList, cutoff = 0.1, minCell = 3){
    gene.mean <- Matrix::rowMeans(cnvList$expr.data)
    gene.sum <- Matrix::rowSums(cnvList$expr.data > 0)
    genes.sel <- rownames(cnvList$expr.data)[gene.mean >= cutoff & gene.sum >= minCell]

    cnvList$expr.data <- cnvList$expr.data[genes.sel, ]
    cnvList$gene.chr <- cnvList$gene.chr[genes.sel, ]

    return(cnvList)
}


normalizeDataForCNV <- function(cnvList){
    expr.data <- cnvList$expr.data

    cs = Matrix::colSums(expr.data)
    expr.data <- t(t(expr.data) / cs)
    normalize_factor <- 10^round(log10(mean(cs)))
    expr.data <- expr.data * normalize_factor

    cnvList$expr.data <- expr.data
    return(cnvList)
}



anscombeTransform <- function(cnvList){
    cnvList$expr.data <- 2 * sqrt(cnvList$expr.data + 3/8)
    return(cnvList)
}



logForCNV <- function(cnvList){
    cnvList$expr.data <- log2(cnvList$expr.data + 1)
    return(cnvList)
}



getAverageBounds <- function(cnvList){
    lower.bound <- mean(apply(cnvList$expr.data, 2, min))
    upper.bound <- mean(apply(cnvList$expr.data, 2, max))
    threshold = mean(abs(c(lower.bound, upper.bound)))
    return(threshold)
}



boundForCNV <- function(cnvList, threshold){
    cnvList$expr.data[cnvList$expr.data > threshold] <- threshold
    cnvList$expr.data[cnvList$expr.data < (-1 * threshold)] <- -1 * threshold
    return(cnvList)
}



smoothOne <- function(ori.data, window.len = window.len){
    half.window <- (window.len - 1) / 2

    pad.data <- c(rep(0, half.window), ori.data, rep(0, half.window))
    bool.data <- c(rep(0, half.window), rep(1, length(ori.data)), rep(0, half.window))

    kernel.vec <- c(1:half.window, half.window + 1, half.window:1)

    sum.data <- filter(pad.data, kernel.vec, sides = 2)
    num.data <- filter(bool.data, kernel.vec, sides = 2)
    sum.data <- sum.data[!is.na(sum.data)]
    num.data <- num.data[!is.na(num.data)]

    smo.data <- sum.data / num.data
    return(smo.data)
}


smoothByChr <- function(cnvList, window.len = 101){
    chrList <- cnvList$gene.chr$CHR
    chrs <- as.character(unique(cnvList$gene.chr$CHR))

    if(window.len < 2){
        cat("- Warning in 'smoothBychr': Window length < 2, returning original data.\n")
        return(cnvList)
    }

    expr.data <- cnvList$expr.data
    for(chr in chrs) {
        # print(chr)
        cur.genes.ix <- which(chrList == chr)
        cur.data <- expr.data[cur.genes.ix, , drop=F]

        if(length(cur.genes.ix) > 1) {
            if(window.len %% 2 == 0){
                window.len = window.len + 1
                cat("- Warning in 'smoothBychr': Window length is even, adding one to 'window.len'.\n")
            }

            smooth.data <- apply(cur.data, 2, smoothOne, window.len = window.len)
            rownames(smooth.data) <- rownames(cur.data)
            colnames(smooth.data) <- colnames(cur.data)

            expr.data[cur.genes.ix, ] <- smooth.data
        }
    }
    cnvList$expr.data <- expr.data
    return(cnvList)
}



centerAcrossChr <- function(cnvList, method = "median"){
    expr.data <- cnvList$expr.data
    if (method == "median") {
        row_median <- apply(expr.data, 2, function(x) { median(x, na.rm=T) } )
        expr.data <- t(apply(expr.data, 1, "-", row_median))
    }
    else {
        row_means <- apply(expr.data, 2, function(x) { mean(x, na.rm=T) } )
        expr.data <- t(apply(expr.data, 1, "-", row_means))
    }
    cnvList$expr.data <- expr.data
    return(cnvList)
}



subtractRefExpr <- function(cnvList, inv_log = TRUE){
    ref.cellNames <- subset(cnvList$cell.anno, cellAnno == "Reference")$cellName

    if (inv_log) {
        ref.means = log2(Matrix::rowMeans(2^cnvList$expr.data[, ref.cellNames] - 1) + 1)
    } else {
        ref.means = Matrix::rowMeans(cnvList$expr.data[, ref.cellNames])
    }

    cnvList$expr.data <- cnvList$expr.data - ref.means
    return(cnvList)
}


invertLog2 <- function(cnvList){
    cnvList$expr.data <- 2^cnvList$expr.data
    return(cnvList)
}



denoiseByRefMeanSd <- function(cnvList, sd_amplifier=1.5){
    expr.data <- cnvList$expr.data
    ref.cellNames <- subset(cnvList$cell.anno, cellAnno == "Reference")$cellName

    mean.ref.vals <- mean(expr.data[, ref.cellNames])
    mean.ref.sd <- mean(apply(expr.data[, ref.cellNames], 2, function(x) sd(x, na.rm=T))) * sd_amplifier

    up.bound <- mean.ref.vals + mean.ref.sd
    low.bound <- mean.ref.vals - mean.ref.sd

    expr.data[expr.data > low.bound & expr.data < up.bound] <- mean.ref.vals

    cnvList$expr.data <- expr.data

    return(cnvList)
}



removeOutliers <- function(cnvList){
    expr.data <- cnvList$expr.data
    up.bound <- mean(apply(expr.data, 2, max))
    low.bound <- mean(apply(expr.data, 2, min))

    expr.data[expr.data < low.bound] <- low.bound
    expr.data[expr.data > up.bound] <- up.bound

    cnvList$expr.data <- expr.data

    return(cnvList)
}



runCNV <- function(expr.data,
                   gene.manifest,
                   cell.annotation,
                   cutoff = 0.1, minCell = 3,
                   ref.data = NULL,
                   species = "human",
                   genome = "hg19",
                   hg.mm.mix = F){

    cnvList <- prepareCNV(expr.data = expr.data,
                          gene.manifest = gene.manifest,
                          cell.annotation,
                          ref.data = ref.data,
                          species = species,
                          genome = genome,
                          hg.mm.mix = hg.mm.mix)
    cnvList <- rmGeneForCNV(cnvList, cutoff = cutoff, minCell = minCell)
    cnvList <- normalizeDataForCNV(cnvList)
    cnvList <- anscombeTransform(cnvList)
    cnvList <- logForCNV(cnvList)
    threshold <- getAverageBounds(cnvList)
    cnvList <- boundForCNV(cnvList, threshold)
    cnvList <- smoothByChr(cnvList, window.len = 101)
    cnvList <- centerAcrossChr(cnvList, method = "median")
    cnvList <- subtractRefExpr(cnvList)
    cnvList <- invertLog2(cnvList)
    cnvList <- denoiseByRefMeanSd(cnvList, sd_amplifier = 1.0)
    cnvList <- removeOutliers(cnvList)

    return(cnvList)
}


getMalignScore <- function(cnvList, cell.type = "Observation", method = "smooth", adjMat = NULL){
    if(cell.type == "Observation"){
        cell.names <- subset(cnvList$cell.anno, cellAnno != "Reference")$cellName
    }else if(cell.type == "Reference"){
        cell.names <- subset(cnvList$cell.anno, cellAnno == "Reference")$cellName
    }

    cur.data <- cnvList$expr.data[, cell.names]

    if(is.null(adjMat) & method == "smooth"){
        cat("- Warning in 'getMalignScore': Adjacent matrix is not provided, and use 'direct' method instead.\n")
        method <- "direct"
    }
    if(method == "smooth"){
        thres <- quantile(adjMat@x, 1- (dim(adjMat)[1] * 10 / length(adjMat@x)))

        indexes <- as.matrix((adjMat > thres) + 0)
        tt <- 0.5 / (rowSums(indexes) - 1)
        tt[is.infinite(tt)] <- 0

        indexes <- indexes * tt
        indexes <- indexes * (1 - diag(rep(1, dim(indexes)[1])))
        diagValue <- rep(0.5, dim(indexes)[1])
        diagValue[tt == 0] <- 1

        indexes <- t(indexes + diag(diagValue))

        new.cur.data <- as.matrix(cur.data) %*% indexes
        malignScore <- colSums((new.cur.data - 1)^2)
        malignScore <- malignScore / dim(new.cur.data)[1]

    }else if(method == "direct"){
        malignScore <- colSums((cur.data - 1)^2)
        malignScore <- malignScore / dim(cur.data)[1]
    }

    names(malignScore) <- colnames(cur.data)

    return(malignScore)
}



malignPlot <- function(obserScore, referScore, malign.thres = NULL){
    scoreDF <- data.frame(malignScore = c(obserScore, referScore),
                          sets = c(rep("Observation", length(obserScore)),
                                   rep("Reference", length(referScore))))
    p <- ggplot() +
        geom_histogram(data = subset(scoreDF, sets == "Observation"),
                       mapping = aes(x = malignScore, fill = "Observation"),
                       bins = 150, alpha = 0.6) +
        geom_histogram(data = subset(scoreDF, sets == "Reference"),
                       mapping = aes(x = malignScore, fill = "Reference"),
                       bins = 150, alpha = 0.6) +
        labs(x = "Malignancy score", y = "Droplets count") +
        scale_fill_manual(name = "Cells sets", guide = "legend",
                          values = c("Observation"="#2e68b7", "Reference"="grey")) +
        theme_classic() +
        ggplot_config(base.size = 7) +
        theme(legend.justification = c(1.12,1.12), legend.position = c(1,1))
    if(!is.null(malign.thres)){
        p <- p + geom_vline(xintercept = malign.thres, colour = "red", linetype = "dashed")
    }
    return(p)
}



getBimodalThres <- function(scores){
    x.density <- density(scores)
    d.x.density <- diff(x.density$y)
    d.sign <- (d.x.density > 0) + 0

    ext.pos <- which(d.sign[2:length(d.sign)] - d.sign[1:(length(d.sign)-1)] != 0)
    ext.density <- x.density$y[ext.pos]
    y.max <- max(ext.density)
    if(length(ext.pos) >= 3){
        del.ix <- c()
        for(ei in 2:length(ext.density)){
            if(abs(ext.density[ei] - ext.density[ei - 1]) < y.max * 0.001){
                del.ix <- c(del.ix, ei - 1, ei)
            }
        }
        sel.ix <- !(1:length(ext.density) %in% unique(del.ix))
        ext.density <- ext.density[sel.ix]
        ext.pos <- ext.pos[sel.ix]
    }

    if(length(ext.pos) >= 3){
        t.ext.density <- c(0, ext.density, 0)
        ext.height <- sapply(2:(length(ext.pos) + 1), FUN = function(x){
            return(min(abs(t.ext.density[x] - t.ext.density[x-1]), abs(t.ext.density[x] - t.ext.density[(x+1)])))
        })
        ext <- data.frame(x = ext.pos, y = ext.density, height = ext.height)
        max.ix <- order(ext.density, decreasing = T)
        if(ext.height[max.ix[2]] / ext.height[max.ix[1]] > 0.01){
            cut.df <- ext[c(max.ix[2]:max.ix[1]), ]
            threshold <- x.density$x[cut.df[which.min(cut.df$y), ]$x]
        }else{
            threshold <- NULL
        }
    }else{
        threshold <- NULL
    }

    return(threshold)
}


#
# getBimodalThres <- function(scores){
#     x.density <- density(scores)
#     d.x.density <- diff(x.density$y)
#     d.sign <- (d.x.density > 0) + 0
#
#     ext.pos <- which(d.sign[2:length(d.sign)] - d.sign[1:(length(d.sign)-1)] != 0)
#     if(length(ext.pos) >= 3){
#         ext.density <- x.density$y[ext.pos]
#         t.ext.density <- c(0, ext.density, 0)
#         ext.height <- sapply(2:(length(ext.pos) + 1), FUN = function(x){
#             return(min(abs(t.ext.density[x] - t.ext.density[x-1]), abs(t.ext.density[x] - t.ext.density[(x+1)])))
#         })
#         ext <- data.frame(x = ext.pos, y = ext.density, height = ext.height)
#
#         max.ix <- order(ext.density, decreasing = T)
#         if(ext.height[max.ix[2]] / ext.height[max.ix[1]] > 0.1){
#             cut.df <- ext[c(max.ix[2]:max.ix[1]), ]
#             threshold <- x.density$x[cut.df[which.min(cut.df$y), ]$x]
#         }else{
#             threshold <- NULL
#         }
#     }else{
#         threshold <- NULL
#     }
#     return(threshold)
# }



#' plotMalignancy
#'
#' @param cell.annotation A data.frame of cells' annotation containing the cells'
#' malignancy score (`Malign.score`) and type (`Malign.type`).
#' @inheritParams runScAnnotation
#'
#' @return A plot list.
#' @export
#'
plotMalignancy <- function(cell.annotation,
                           malignancy.method,
                           coor.names = c("tSNE_1", "tSNE_2"),
                           savePath = NULL){
    ## scatter plot of malignancy
    p.malignType.Point <- pointDRPlot(cell.annotation, value = "Malign.type",
                                      coor.names = coor.names,
                                      # colors = c("malignant" = "#f57e87", "nonMalignant" = "#66d5a5"),
                                      colors = c("malignant" = "#9467BD", "nonMalignant" = '#2CA02C'),
                                      legend.position = "right",
                                      legend.title = "Malignancy\n type")

    p.malignScore.Point <- pointDRPlot(cell.annotation, value = "Malign.score",
                                       coor.names = coor.names,
                                       # colors = c("white", "#f57e87"),
                                       colors = c("white", "#9467BD"),
                                       discrete = F,
                                       limit.quantile = 0.1,
                                       legend.position = "right",
                                       legend.title = "Malignancy\n score")

    p.malignType.bar <- clusterBarPlot(cell.annotation = cell.annotation,
                                       # cell.colors = c("malignant" = "#f57e87", "nonMalignant" = "#66d5a5"),
                                       cell.colors = c("malignant" = "#9467BD", "nonMalignant" = '#2CA02C'),
                                       sel.col = "Malign.type",
                                       legend.title = "Malignancy type")

    ## save
    if(!is.null(savePath)){
        ggsave(filename = file.path(savePath, "figures",
                                    paste0(malignancy.method, "-malignType-point.png")),
               p.malignType.Point, width = 5, height = 3.8, dpi = 300)
        ggsave(filename = file.path(savePath, "figures",
                                    paste0(malignancy.method, "-malignScore-point.png")),
               p.malignScore.Point, width = 5, height = 3.8, dpi = 300)
        ggsave(filename = file.path(savePath, "figures",
                                    paste0(malignancy.method, "-malignType-bar.png")),
               p.malignType.bar, width = 6, height = 3, dpi = 300)
    }

    return(list(p.malignType.Point = p.malignType.Point,
                p.malignScore.Point = p.malignScore.Point,
                p.malignType.bar = p.malignType.bar))
}




#' runMalignancy
#'
#' @param expr A Seurat object.
#' @param gene.manifest A data.frame of genes' manifest.
#' @param cell.annotation A data.frame of cells' annotation.
#' @param cutoff The cut-off for min average read counts per gene among
#' reference cells. The default is 0.1.
#' @param minCell An integer number used to filter gene. The default is 3.
#' @param p.value.cutoff The p-value to decide whether the distribution of
#' malignancy score is bimodality.
#' @param ref.data An expression matrix of gene by cell, which is used as the normal reference.
#' The default is NULL, and an immune cells or bone marrow cells expression matrix will be used for human or mouse species, respectively.
#' @param referAdjMat An adjacent matrix for the normal reference data.
#' The larger the value, the closer the cell pair is.
#' The default is NULL, and a SNN matrix of the default ref.data will be used.
#' @inheritParams runScAnnotation
#'
#' @return A list of cnvList, reference malignancy score, seurat object,
#' cell.annotatino, bimodal.pvalue, malign.thres, and all generated plots.
#' @export
#'
runMalignancy <- function(expr,
                          gene.manifest,
                          cell.annotation,
                          savePath,
                          cutoff = 0.1, minCell = 3,
                          p.value.cutoff = 0.5,
                          coor.names = c("tSNE_1", "tSNE_2"),
                          ref.data = NULL,
                          referAdjMat = NULL,
                          species = "human",
                          genome = "hg19",
                          hg.mm.mix = F){
    if(!dir.exists(file.path(savePath, 'malignancy-inferCNV/'))){
        dir.create(file.path(savePath, 'malignancy-inferCNV/'), recursive = T)
    }

    expr.data <- expr@assays$RNA@counts
    cnvList <- runCNV(expr.data = expr.data,
                      gene.manifest = gene.manifest,
                      cell.annotation = cell.annotation,
                      cutoff = cutoff, minCell = minCell,
                      ref.data = ref.data,
                      species = species,
                      genome = genome,
                      hg.mm.mix = hg.mm.mix)

    if(is.null(ref.data)){
        if(species == "human"){
            referAdjMat <- readRDS(system.file("rds", "cnvRef_SNN-HM.RDS", package = "scCancer"))
        }else if(species == "mouse"){
            referAdjMat <- readRDS(system.file("rds", "cnvRef_SNN-boneMarrow-MS.RDS", package = "scCancer"))
        }
    }
    referScore.smooth <- getMalignScore(cnvList, "Reference", method = "smooth", adjMat = referAdjMat)
    obserScore.smooth <- getMalignScore(cnvList, "Observation", method = "smooth",
                                        adjMat = expr@graphs$RNA_snn)
    up.refer <- quantile(referScore.smooth, 0.995)
    low.refer <- quantile(referScore.smooth, 0.005)
    referScore.smooth <- (referScore.smooth - low.refer) / (up.refer - low.refer)
    obserScore.smooth <- (obserScore.smooth - low.refer) / (up.refer - low.refer)

    all.thres <- getBimodalThres(scores = c(referScore.smooth, obserScore.smooth))
    malign.thres <- getBimodalThres(scores = obserScore.smooth)

    ju.exist.malign <- !is.null(all.thres) | !is.null(malign.thres)

    ## malignancy type
    if(!is.null(all.thres)){
        malign.type <- rep("malignant", length(obserScore.smooth))
        names(malign.type) <- names(obserScore.smooth)
        if(!is.null(malign.thres)){
            malign.type[names(obserScore.smooth)[obserScore.smooth < malign.thres]] <- "nonMalignant"
        }
    }else{
        malign.type <- rep("nonMalignant", length(obserScore.smooth))
        names(malign.type) <- names(obserScore.smooth)
        if(!is.null(malign.thres)){
            malign.type[names(obserScore.smooth)[obserScore.smooth >= malign.thres]] <- "malignant"
        }
    }
    p.malignScore <- malignPlot(obserScore.smooth, referScore.smooth,
                                malign.thres = malign.thres)

    ## add score and type to cell.annotation
    cell.annotation$Malign.score <- obserScore.smooth[rownames(cell.annotation)]
    cell.annotation$Malign.type <- malign.type[rownames(cell.annotation)]
    expr[["Malign.score"]] <- cell.annotation$Malign.score
    expr[["Malign.type"]] <- cell.annotation$Malign.type

    ## plot
    p.results <- plotMalignancy(cell.annotation = cell.annotation,
                                malignancy.method = "inferCNV",
                                coor.names = coor.names,
                                savePath = savePath)
    p.results[["p.malignScore"]] <- p.malignScore
    ggsave(filename = file.path(savePath, "figures/inferCNV-malignScore.png"),
           p.malignScore, width = 5, height = 4, dpi = 300)

    ## save results
    write.table(cnvList$expr.data[, names(obserScore.smooth)],
                file = file.path(savePath, "malignancy-inferCNV/inferCNV-observation.txt"),
                quote = F, sep = "\t", row.names = T)
    write.table(cnvList$expr.data[, names(referScore.smooth)],
                file = file.path(savePath, "malignancy-inferCNV/inferCNV-reference.txt"),
                quote = F, sep = "\t", row.names = T)
    write.table(data.frame(referScore.smooth),
                file = file.path(savePath, "malignancy-inferCNV/refer-malignScore.txt"),
                quote = F, sep = "\t", row.names = T)

    results <- list(
        cnvList = cnvList,
        referScore = referScore.smooth,
        expr = expr,
        cell.annotation = cell.annotation,
        ju.exist.malign = ju.exist.malign,
        # bimodal.pvalue = bimodal.pvalue,
        malign.thres = malign.thres,
        p.results = p.results
    )
    return(results)
}


