prepareCNV <- function(dataPath, statPath,
                       cell.annotation,
                       hg.mm.mix = F,
                       species = "human"){
    ## gene.chr
    if(species == "human"){
        gene.chr <- read.table(system.file("txt", "gene.chr.txt", package = "scCancer"),
                               col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
                               stringsAsFactors = F)
    }else{
        stop("Error in 'runInferCNV': ", species, " is not 'human'.\n")
    }

    ## expr.data
    exprList <- getFilterData(
        dataPath = dataPath,
        statPath = statPath,
        bool.filter.cell = T,
        bool.filter.gene = T,
        anno.filter = c(),
        nCell.min = 3,
        bgPercent.max = 1,
        hg.mm.mix = hg.mm.mix
    )
    expr.data <- exprList$expr.data
    gene.manifest <- exprList$gene.manifest
    rownames(gene.manifest) <- gene.manifest$Symbol


    ## reference.data
    ref.data <- readRDS(system.file("rds", "cnvRefData-immune.RDS", package = "scCancer"))
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
        warning("Window length < 2, returning original data.\n")
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
                warning("Window length is even, adding one to 'window.len'.\n")
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



runCNV <- function(dataPath, statPath, savePath,
                   cell.annotation,
                   cutoff = 0.1, minCell = 3,
                   hg.mm.mix = F){

    # cell.annotation <- read.table(paste0(savePath, "/cellAnnotation.txt"),
    #                               sep = "\t", header = T, stringsAsFactors = F)

    cnvList <- prepareCNV(dataPath = dataPath,
                          statPath = statPath,
                          cell.annotation,
                          hg.mm.mix = hg.mm.mix,
                          species = "human")
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

    if(method == "smooth"){
        if(is.null(adjMat)){
            stop("Please provide adjacent matrix.\n")
        }else{
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
        }

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
        geom_vline(aes(xintercept = malign.thres), colour = "red", linetype = "dashed") +
        labs(x = "Malignancy score", y = "Droplets count") +
        scale_fill_manual(name = "Cells sets", guide = "legend",
                          values = c("Observation"="#2e68b7", "Reference"="grey")) +
        ggplot_config(base.size = 8) +
        theme_classic() +
        theme(legend.justification = c(1.12,1.12), legend.position = c(1,1))

    return(p)
}


getMalignThres <- function(malignScore){
    d.t <- dip.test(malignScore)
    p.value <- d.t$p.value

    score.density <- density(malignScore)
    d.score.density <- diff(score.density$y)
    d.sign <- (d.score.density > 0) + 0
    # threshold <- score.density$x[which(d.sign[2:length(d.sign)] - d.sign[1:(length(d.sign)-1)] == 1)[1]+1]

    ext.pos <- which(d.sign[2:length(d.sign)] - d.sign[1:(length(d.sign)-1)] != 0)
    ext.density <- score.density$y[ext.pos]

    if(length(ext.pos) >= 3){
        left.gap <- ext.density[1:(length(ext.density)-2)] - ext.density[2:(length(ext.density)-1)]
        right.gap <- ext.density[3:length(ext.density)] - ext.density[2:(length(ext.density)-1)]
        thres.pos <- ext.pos[which.max(apply(data.frame(left.gap, right.gap), 1, min)) + 1] + 1
    }else{
        if(which(d.sign[2:length(d.sign)] - d.sign[1:(length(d.sign)-1)] == 1) > 0){
            thres.pos <- ext.pos[1]
        }else{
            thres.pos <- 1
        }
    }
    return(list(threshold = score.density$x[thres.pos],
                p.value = p.value))
}



#' runMalignancy
#'
#' @param cell.annotation A data.frame of cells' annotation
#' @param expr A Seurat object.
#' @param cutoff The cut-off for min average read counts per gene among
#' reference cells. The default is 0.1)
#' @param minCell An integer number used to filter gene. The default is 3.
#' @param p.value.cutoff The p-value to decide whether the distribution of
#' malignancy score is bimodality.
#' @inheritParams runScAnnotation
#'
#' @return A list of cnvList, reference malignancy score, seurat object,
#' cell.annotatino, bimodal.pvalue, malign.thres, and all generated plots.
#' @export
#'
#' @importFrom diptest dip.test
#'
#' @examples
runMalignancy <- function(dataPath, statPath, savePath,
                          cell.annotation,
                          expr,
                          cutoff = 0.1, minCell = 3,
                          p.value.cutoff = 0.5,
                          hg.mm.mix = F){
    if(!dir.exists(file.path(savePath, 'malignancy/'))){
        dir.create(file.path(savePath, 'malignancy/'), recursive = T)
    }

    cnvList <- runCNV(dataPath = dataPath,
                      statPath = statPath,
                      savePath = savePath,
                      cell.annotation = cell.annotation,
                      cutoff = cutoff, minCell = minCell,
                      hg.mm.mix = hg.mm.mix)

    referAdjMat <- readRDS(system.file("rds", "cnvRef_SNN-immune.RDS", package = "scCancer"))
    referScore.smooth <- getMalignScore(cnvList, "Reference", method = "smooth", adjMat = referAdjMat)
    obserScore.smooth <- getMalignScore(cnvList, "Observation", method = "smooth",
                                        adjMat = expr@graphs$RNA_snn)

    ## malignant type
    ju.exist.malign <- dip.test(c(referScore.smooth, obserScore.smooth))$p.value < 0.95
    tmp <- getMalignThres(obserScore.smooth)
    malign.thres <- tmp$threshold
    bimodal.pvalue <- tmp$p.value

    ## score hist plot
    if(ju.exist.malign){
        malign.type <- rep("malignant", length(obserScore.smooth))
        names(malign.type) <- names(obserScore.smooth)

        if(bimodal.pvalue < p.value.cutoff){
            p.malignScore <- malignPlot(obserScore.smooth, referScore.smooth,
                                        malign.thres = malign.thres)
            malign.type[names(obserScore.smooth)[obserScore.smooth < malign.thres]] <- "nonMalignant"
        }else{
            p.malignScore <- malignPlot(obserScore.smooth, referScore.smooth,
                                        malign.thres = min(obserScore.smooth))
        }
    }else{
        p.malignScore <- malignPlot(obserScore.smooth, referScore.smooth,
                                    malign.thres = max(obserScore.smooth))
        malign.type <- rep("nonMalignant", length(obserScore.smooth))
        names(malign.type) <- names(obserScore.smooth)
    }

    ## add score and type to cell.annotation
    cell.annotation$Malign.score <- obserScore.smooth[rownames(cell.annotation)]
    cell.annotation$Malign.type <- malign.type[rownames(cell.annotation)]
    expr[["Malign.score"]] <- cell.annotation$Malign.score
    expr[["Malign.type"]] <- cell.annotation$Malign.type

    ## scatter plot of malignancy
    p.malignType.Point <- pointDRPlot(cell.annotation, value = "Malign.type",
                                      colors = c("malignant" = "#f57e87", "nonMalignant" = "#66d5a5"),
                                      legend.position = "bottom",
                                      legend.title = "Malignancy type")

    p.malignScore.Point <- pointDRPlot(cell.annotation, value = "Malign.score",
                                       colors = c("white", "#f57e87"),
                                       discrete = F,
                                       limit.quantile = 0.1,
                                       legend.position = "bottom",
                                       legend.title = "Malignancy score")

    ## save
    ggsave(filename = file.path(savePath, "figures/malignScore.png"),
           p.malignScore, width = 8, height = 4, dpi = 800)
    ggsave(filename = file.path(savePath, "figures/malignType-point.png"),
           p.malignType.Point, width = 6, height = 6.2, dpi = 800)
    ggsave(filename = file.path(savePath, "figures/malignScore-point.png"),
           p.malignScore.Point, width = 6, height = 6.2, dpi = 800)

    write.table(cnvList$expr.data[, names(obserScore.smooth)],
                file = file.path(savePath, "malignancy/inferCNV-observation.txt"),
                quote = F, sep = "\t", row.names = T)
    write.table(cnvList$expr.data[, names(referScore.smooth)],
                file = file.path(savePath, "malignancy/inferCNV-reference.txt"),
                quote = F, sep = "\t", row.names = T)
    write.table(data.frame(referScore.smooth),
                file = file.path(savePath, "malignancy/refer-malignScore.txt"),
                quote = F, sep = "\t", row.names = T)

    results <- list(
        cnvList = cnvList,
        referScore = referScore.smooth,
        expr = expr,
        cell.annotation = cell.annotation,
        ju.exist.malign = ju.exist.malign,
        bimodal.pvalue = bimodal.pvalue,
        malign.thres = malign.thres,
        p.results = list(p.malignScore = p.malignScore,
                         p.malignType.Point = p.malignType.Point,
                         p.malignScore.Point = p.malignScore.Point)
    )
    return(results)
}


