
get10Xpath <- function (samplePath, raw.data = F){
    prefix <- ifelse(raw.data, 'raw', 'filtered')
    cur.path <- paste0(samplePath, '/')
    res.path <- paste0(cur.path, prefix, '_feature_bc_matrix')
    if (!dir.exists(res.path)){
        res.path <- paste0(cur.path, prefix, '_gene_bc_matrices/hg19/')
    }
    if (!dir.exists(res.path)){
        res.path <- paste0(cur.path, prefix, '_gene_bc_matrices/hg38/')
    }
    if (!dir.exists(res.path)){
        res.path <- paste0(cur.path, prefix, '_gene_bc_matrices/mm10/')
    }
    if(!dir.exists(res.path)){
        res.path <- NULL
    }
    return(res.path)
}


ExtractField <- function (string, field = 1, delim = "_"){
    fields <- as.numeric(x = unlist(x = strsplit(
        x = as.character(x = field), split = ",")))
    if (length(x = fields) == 1) {
        return(strsplit(x = string, split = delim)[[1]][field])
    }
    return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}


getCRversion <- function(data.path){
    version <- "Cell Ranger (version 2)"
    if(grepl("feature_bc_matrix", data.path)){
        version <- "Cell Ranger (version >= 3)"
    }
    return(version)
}


getBarcodes <- function(data.path){
    barcode.loc <- paste0(data.path, "/barcodes.tsv")
    if(grepl("feature_bc_matrix", data.path)){
        barcode.loc <- paste0(barcode.loc, ".gz")
    }
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
        cell.names <- as.vector(x = as.character(
            x = sapply(
                X = cell.names,
                FUN = ExtractField,
                field = 1,
                delim = "-"
            )
        ))
    }
    return(cell.names)
}


#' Read10Xdata
#'
#' Read expression matrix data from 10X. This function is modified from Seurat package.
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files provided by 10X.
#' A vector or named vector can be given in order to load several data directories.
#' If a named vector is given, the cell barcode names will be prefixed with the name.
#' @param gene.column An integer indicating which column of genes.tsv or features.tsv to use for gene names; default is 2.
#' @param unique.features Make feature names unique (default TRUE).
#' @param only.expr Whether to read expression data only if have multiple features (default TRUE).
#'
#' @return If the 10X data only has expression data or the argument 'only.expr' is TRUE,
#' a sparse matrix containing the expression data will be returned.
#' Otherwise, if the 10X data has multiple data types,
#' a list containing a sparse matrix of the data from each type will be returned.
#'
#' @export
#'
Read10Xdata <- function (data.dir = NULL, gene.column = 2,
                         unique.features = TRUE, only.expr = TRUE)  {
    full.data <- list()
    for (i in seq_along(data.dir)) {
        run <- data.dir[i]
        if (!dir.exists(paths = run)) {
            stop("Directory provided does not exist")
        }
        if (!grepl("\\/$", run)) {
            run <- paste(run, "/", sep = "")
        }
        barcode.loc <- file.path(run, "barcodes.tsv")
        gene.loc <- file.path(run, "genes.tsv")
        features.loc <- file.path(run, "features.tsv.gz")
        matrix.loc <- file.path(run, "matrix.mtx")
        pre_ver_3 <- file.exists(gene.loc)
        if (!pre_ver_3) {
            addgz <- function(s) {
                return(paste0(s, ".gz"))
            }
            barcode.loc <- addgz(s = barcode.loc)
            matrix.loc <- addgz(s = matrix.loc)
        }
        if (!file.exists(barcode.loc)) {
            stop("Barcode file missing")
        }
        if (!pre_ver_3 && !file.exists(features.loc)) {
            stop("Gene name or features file missing")
        }
        if (!file.exists(matrix.loc)) {
            stop("Expression matrix file missing")
        }
        data <- readMM(file = matrix.loc)
        cell.names <- readLines(barcode.loc)
        if (all(grepl(pattern = "\\-1$", x = cell.names))) {
            cell.names <- as.vector(x = as.character(
                x = sapply(
                    X = cell.names,
                    FUN = ExtractField,
                    field = 1,
                    delim = "-"
                )
            ))
        }
        if (is.null(x = names(x = data.dir))) {
            if (i < 2) {
                colnames(x = data) <- cell.names
            } else {
                colnames(x = data) <- paste0(i, "_", cell.names)
            }
        } else {
            colnames(x = data) <-
                paste0(names(x = data.dir)[i], "_", cell.names)
        }
        feature.names <- read.delim(
            file = ifelse(
                test = pre_ver_3,
                yes = gene.loc,
                no = features.loc
            ),
            header = FALSE,
            stringsAsFactors = FALSE
        )
        if (any(is.na(x = feature.names[, gene.column]))) {
            warning("Some features names are NA. Replacing NA names with ID from the opposite column requested",
                    call. = FALSE, immediate. = TRUE)
            na.features <- which(x = is.na(x = feature.names[,
                                                             gene.column]))
            replacement.column <- ifelse(test = gene.column ==
                                             2, yes = 1, no = 2)
            feature.names[na.features, gene.column] <- feature.names[na.features,
                                                                     replacement.column]
        }
        if (unique.features) {
            fcols = ncol(x = feature.names)
            if (fcols < gene.column) {
                stop(paste0("gene.column was set to ", gene.column,
                            " but feature.tsv.gz (or genes.tsv) only has ",
                            fcols, " columns.", " Try setting the gene.column argument to a value <= to ",
                            fcols, "."))
            }
            rownames(x = data) <- make.unique(names = feature.names[,
                                                                    gene.column])
        }
        # In cell ranger 3.0, a third column specifying the type of data was added
        # and we will return each type of data as a separate matrix
        if (ncol(x = feature.names) > 2) {
            data_types <- factor(x = feature.names$V3)
            lvls <- levels(x = data_types)
            if (length(x = lvls) > 1 && length(x = full.data) ==
                0) {
                message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
            }
            expr_name <- "Gene Expression"
            if (expr_name %in% lvls) {
                # Return Gene Expression first
                lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
            }
            data <- lapply(
                X = lvls,
                FUN = function(l) {
                    return(data[data_types == l, ])
                }
            )
            names(x = data) <- lvls
        } else{
            data <- list(data)
        }
        full.data[[length(x = full.data) + 1]] <- data
    }
    # Combine all the data from different directories into one big matrix, note this
    # assumes that all data directories essentially have the same features files
    list_of_data <- list()
    for (j in 1:length(x = full.data[[1]])) {
        list_of_data[[j]] <-
            do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
        list_of_data[[j]] <-
            as(object = list_of_data[[j]], Class = "CsparseMatrix")
    }
    names(x = list_of_data) <- names(x = full.data[[1]])

    if (only.expr){
        return(list_of_data[[1]])
    }else{
        # If multiple features, only return a list, otherwise a matrix.
        if (length(x = list_of_data) == 1) {
            return(list_of_data[[1]])
        } else {
            return(list_of_data)
        }
    }
}



#' ggplot_config
#'
#' @param base.size The size of text.
#'
#' @return A theme.
#' @export
#'
ggplot_config <- function(base.size = 8){
    p <- theme_classic() +
        theme(plot.title = element_text(size = 2 * base.size),
              axis.title.x = element_text(size = 2 * base.size, vjust = -0.2),
              axis.title.y = element_text(size = 2 * base.size, vjust = 0.2),
              axis.text.x = element_text(size = 2 * base.size),
              axis.text.y = element_text(size = 2 * base.size),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.title = element_text(size = 2 * base.size - 2),
              legend.text = element_text(size = 1.5 * base.size)
        )
    return(p)
}



getOutliers <- function(x){
    x.med <- median(x)
    outs <- boxplot.stats(x)$out
    outliers <- subset(outs, outs > x.med)
    return(outliers)
}



getCellix <- function(cell.manifest, filter.thres, arg){
    ixs <- lapply(arg, FUN = function(x) {
        ix <- which(cell.manifest[[x]] >= filter.thres[x, 'Low.threshold'] &
                    cell.manifest[[x]] < filter.thres[x, 'High.threshold'])
        return(ix)
    })
    res.ix <- ixs[[1]]
    for(i in 1:length(arg)){
        res.ix <- intersect(res.ix, ixs[[i]])
    }
    return(res.ix)
}



grid_arrange_shared_legend <- function(..., all.p, ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(all.p + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
    grid.newpage()
    grid.draw(combined)

    # return gtable invisibly
    invisible(combined)
}



#' getDefaultMarkers
#'
#' Return default markers of several common cell types.
#'
#' @inheritParams runScAnnotation
#'
#' @return A list of default markers of several common cell types.
#' @export
#'
getDefaultMarkers <- function(species = "human"){
    # feature.def <- list(
    #     "T cell" = c("CD3D"),
    #     "B cell" = c("CD79A"),
    #     "NK cell" = c("NKG7"),
    #     "Monocyte" = c("LYZ"),
    #     "Endothelial" = c("PLVAP"),
    #     "Myofibroblast" = c("ACTA2"),
    #     "Epithelial" = c("EPCAM", "KRT8"))

    if(species == "human"){
        feature.def <- list(
            "T cell" = c("PTPRC", "CD3D", "CD4", "CD8A", "CD8B"),
            "B cell" = c("CD79A"),
            "NK cell" = c("NKG7"),
            "Myeloid cell" = c("LYZ"),
            "Endothelial" = c("PLVAP"),
            "Fibroblast" = c("ACTA2"),
            "Epithelial" = c("EPCAM", "KRT8"))
    }else if(species == "mouse"){
        feature.def <- list(
            "T cell" = c("Ptprc", "Cd3d", "Cd4", "Cd8a", "Cd8b"),
            "B cell" = c("Cd79a"),
            "NK cell" = c("Nkg7"),
            "Myeloid cell" = c("Lyz1", "Lyz2"),
            "Endothelial" = c("Plvap"),
            "Fibroblast" = c("Acta2"),
            "Epithelial" = c("Epcam", "Krt8"))
    }

    return(feature.def)
}



#' getDefaultColors
#'
#' @param n The number of colors.
#' @param type The type of color style. Only 1, 2, or 3 is allowed.
#'
#' @return A vector of colors.
#' @export
#'
getDefaultColors <- function(n = NULL, type = 1){
    if(type == 1){
        colors <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb",
                    "#d7652d", "#7cd5c8", "#c49a3f", "#507d41", "#5d8d9c",
                    "#90353b", "#674c2a", "#1B9E77", "#c5383c", "#0081d1",
                    "#ffd900", "#502e71", "#c8b693", "#aed688", "#f6a97a",
                    "#c6a5cc", "#798234", "#6b42c8", "#cf4c8b", "#666666",
                    "#feb308", "#ff1a1a", "#1aff1a", "#1a1aff", "#ffff1a")
    }else if(type == 2){
        if(n <= 8){
            colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                        "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
        }else if(n <= 14){
            colors <- c("#437BFE", "#FEC643", "#43FE69", "#FE6943", "#C643FE",
                        "#43D9FE", "#B87A3D", "#679966", "#993333", "#7F6699",
                        "#E78AC3", "#333399", "#A6D854", "#E5C494")
        }
        else if(n <= 20){
            colors <- c("#87b3d4", "#d5492f", "#6bd155", "#683ec2", "#c9d754",
                        "#d04dc7", "#81d8ae", "#d34a76", "#607d3a", "#6d76cb",
                        "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
                        "#4e857e", "#6e282c", "#d293c8", "#393a2a", "#997579")
        }else if(n <= 30){
            colors <- c("#628bac", "#ceda3f", "#7e39c9", "#72d852", "#d849cc",
                        "#5e8f37", "#5956c8", "#cfa53f", "#392766", "#c7da8b",
                        "#8d378c", "#68d9a3", "#dd3e34", "#8ed4d5", "#d84787",
                        "#498770", "#c581d3", "#d27333", "#6680cb", "#83662e",
                        "#cab7da", "#364627", "#d16263", "#2d384d", "#e0b495",
                        "#4b272a", "#919071", "#7b3860", "#843028", "#bb7d91")
        }else{
            colors <- c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
                        "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
                        "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
                        "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
                        "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
                        "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
                        "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
                        "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
                        "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
                        "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c")
        }
    }else if(type == 3){
        # colors <- c("#07a2a4", "#9a7fd1", "#588dd5", "#f5994e",
        #             "#c05050", "#59678c", "#c9ab00", "#7eb00a")
        colors <- c("#c14089", "#6f5553", "#E5C494", "#738f4c", "#bb6240",
                    "#66C2A5", "#2dfd29", "#0c0fdc")
    }
    if(!is.null(n)){
        if(n <= length(colors)){
            colors <- colors[1:n]
        }else{
            step <- 16777200 %/% (n - length(colors)) - 2
            add.colors <- paste0("#", as.hexmode(seq(from = sample(1:step, 1),
                                                     by = step, length.out = (n-length(colors)))))
            colors <- c(colors, add.colors)
        }
    }
    return(colors)
}


#' getCellTypeColor
#'
#' @param cell.types A vector of cell types.
#'
#' @return A vector of colors.
#' @export
#'
getCellTypeColor <- function(cell.types){
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
    new.types <- setdiff(cell.types, names(cell.colors))
    for(ct in new.types){
        cell.colors[ct] <- getDefaultColors(n = length(new.types), type = 3)[cti]
        cti = cti + 1
    }
    return(cell.colors)
}


limitData <- function(data, min = NULL, max = NULL){
    data2 <- data
    if(!is.null(min)){
        data2[data2 < min] <- min
    }
    if(!is.null(max)){
        data2[data2 > max] <- max
    }
    return(data2)
}



getClusterInfo <- function(cell.annotation){
    cluster.info <- cell.annotation[order(cell.annotation$Cluster), 'Cluster', drop = F]
    cluster.info$Cluster <- as.factor(cluster.info$Cluster)

    num.cluster <- table(cluster.info$Cluster)
    # num.cluster <- num.cluster[as.character(1 : length(num.cluster))]
    num.cluster <- num.cluster[as.character(unique(cluster.info$Cluster))]
    cluster.pos <- cumsum(num.cluster)

    def.colors <- getDefaultColors()
    clusters <- unique(cell.annotation$Cluster)
    clusters <- sort(clusters)
    cluster.colors <- c()
    for(i in 1:length(clusters)){
        cluster.colors[as.character(clusters[i])] = def.colors[clusters[i]]
    }
    cluster.colors = list(Cluster = cluster.colors)

    return(list(cluster.info = cluster.info,
                cluster.colors = cluster.colors,
                cluster.pos = cluster.pos))
}



getMouseGene <- function(hg.genes, bool.name = F, deduplicate = T){
    hg.mm.HomologyGenes <- read.table(system.file("txt", "hg-mm-HomologyGenes.txt", package = "scCancer"),
                                  header = T, stringsAsFactors = F)
    hg.mm.HomologyGenes <- subset(hg.mm.HomologyGenes, hgGenes %in% hg.genes)

    if(deduplicate){
        hg.num <- table(hg.mm.HomologyGenes$hgGenes)
        hg.mm.HomologyGenes <- subset(hg.mm.HomologyGenes, !(hgGenes %in% names(hg.num)[hg.num > 1]))
        mm.num <- table(hg.mm.HomologyGenes$mmGenes)
        hg.mm.HomologyGenes <- subset(hg.mm.HomologyGenes, !(mmGenes %in% names(mm.num)[mm.num > 1]))
    }

    mm.genes <- hg.mm.HomologyGenes$mmGenes

    if(bool.name){
        names(mm.genes) <- hg.mm.HomologyGenes$hgGenes
    }
    return(mm.genes)
}




#' runSurvival
#'
#' According to the marker genes or signatures expression high/low levels,
#' patient are divided into two groups and then survival analysis is performed.
#' The survival curves can be plotted.
#'
#' @param features The names of marker genes or signatures to be analyzed.
#' @param data The data used to perform survival analysis.
#' It should be an expression or signature matrix with gene or signature by patient.
#' The row names are the features' anmes. The columns are patients' labels.
#' @param surv.time The survival time of patients. It should be in accord with the columns of data.
#' @param surv.event The status indicator of patients. 0=alive, 1=dead. It should be in accord with the columns of data.
#' @param cut.off The percentage threshold to divide patients into two groups.
#' The default is 0.5, which means the patients are divided by median.
#' Other values, such as 0.4, means the first 40 percent patients are set "Low" group
#' and the last 40 percent are set "High" group (the median 20 percent are discarded).
#' @param savePath The path to save the survival plots of genes or signatures (the default is NULL and the plots will be return without saving).
#'
#'
#' @return A list of survival curves plots.
#' @export
#'
#' @import survival survminer
#'
runSurvival <- function(features, data, surv.time, surv.event, cut.off = 0.5, savePath = NULL){
    data <- as.matrix(data)
    cut.off <- min(cut.off, 1 - cut.off)

    ps <- list()
    for(feat in features){
        if(feat %in% rownames(data)){
            dw.thres <- quantile(data[feat, ], cut.off)
            up.thres <- quantile(data[feat, ], 1-cut.off)
            p.df <- data.frame(sample = colnames(data),
                               surv.time = surv.time,
                               surv.event = surv.event)
            p.df$expr <- sapply(data[feat, ], function(x){
                if(x >= up.thres){
                    return("High")
                }else if(x < dw.thres){
                    return("Low")
                }else{
                    return("Med")
                }
            })
            surv.df <<- subset(p.df, expr != "Med")
            surv_object <<- Surv(time = surv.df$surv.time, event = surv.df$surv.event)
            fit <- survfit(surv_object ~ expr, data = surv.df)
            p.surv <- ggsurvplot(fit, pval = TRUE,
                                 palette = c("#f57e87", "#66d5a5"),
                                 legend.title = paste0(feat, ":"))
            if(!is.null(savePath)){
                if(!dir.exists(savePath)){
                    dir.create(savePath, recursive = T)
                }
                ggsave(filename = paste0(savePath, "surv-", feat, ".png"), p.surv$plot,
                       width = 3.5, height = 3.5, dpi = 300)
            }
            ps[[feat]] <- p.surv$plot
        }else{
            cat("- Warning in 'runSurvival':", feat, "not found.\n")
        }
    }
    return(ps)
}


#' generate10Xdata
#'
#' Generate a 10X-like data folder based on the data matrix and gene information,
#' which can be used directly to perform scCancer analysis.
#'
#' @param matrix A gene-cell matrix or data.frame.
#' @param gene.info A data.frame of gene information. It should contain two columns,
#' the first is gene Ensemble ID, and the second is gene symbol.
#' The order of the genes should be consistant with the row order of 'matrix'.
#' @param outPath A path to save the output files.
#' @param overwrite If TRUE and the output file already exists, the file is
#' silently overwritten, otherwise an exception is thrown. The default is "FALSE".
#'
#' @return NULL
#' @export
#'
#' @import Matrix R.utils
#'
generate10Xdata <- function(matrix, gene.info, outPath, overwrite = F){
    if(!dir.exists(paste0(outPath, "/filtered_feature_bc_matrix/"))){
        dir.create(paste0(outPath, "/filtered_feature_bc_matrix/"), recursive = T)
    }

    barcode.gz <- gzfile(paste0(outPath, "/filtered_feature_bc_matrix/barcodes.tsv.gz"), "w")
    write.table(colnames(matrix), barcode.gz, quote = F, col.names = F, row.names = F, sep = "\t")
    close(barcode.gz)

    gene.info[, 3] <- "Gene Expression"
    feature.gz <- gzfile(paste0(outPath, "/filtered_feature_bc_matrix/features.tsv.gz"), "w")
    write.table(gene.info, feature.gz, quote = F, col.names = F, row.names = F, sep = "\t")
    close(feature.gz)

    writeMM(as(as.matrix(matrix),"CsparseMatrix"), file = paste0(outPath, "/filtered_feature_bc_matrix/matrix.mtx"))
    gzip(paste0(outPath, "/filtered_feature_bc_matrix/matrix.mtx"), overwrite = overwrite)
}



#' extractFiles
#'
#' Extract files from each sample's folder and rename them with sample's name.
#'
#' @param savePath A path of samples' result folder.
#' @param sampleNames A vector of samples' names (the subfolder names in 'savePath').
#' @param outputPath A path to saving the extracted reports.
#' @param files The name of files you want to extract. The default is c("report-scStat.html", "report-scAnno.html").
#' @param subfolders The name of subfolders for the files you want to extract. The default is NULL.
#' It can be a character string, which means all files are under the subfolder.
#' It can also be a character string vector with same length as "files", which are corresponding to "files".
#'
#' @return NULL
#' @export
#'
extractFiles <- function(savePath, sampleNames, outputPath,
                         files = c("report-scStat.html", "report-scAnno.html"),
                         subfolders = NULL){
    message("[", Sys.time(), "] -----: extract files")
    if((!is.null(subfolders)) & (length(subfolders) != 1) & (length(subfolders) != length(files))){
        stop("The lengths of files and subfolders are not equal.")
    }

    if(!dir.exists(file.path(outputPath))){
        dir.create(file.path(outputPath), recursive = T)
    }

    for(sampleName in sampleNames){
        cur.path <- paste0(savePath, "/", sampleName, "/")
        ori.files <- paste0(cur.path, subfolders, "/", files)
        new.files <- paste0(outputPath, "/", sampleName, "-", files)
        file.copy(ori.files, new.files, overwrite = T)
    }
}



#' checkStatArguments
#'
#'
#' @param argList A list of arguments passed into 'runScStatistics".
#'
#' @return NULL
#' @export
#'
checkStatArguments <- function(argList){
    if(!dir.exists(argList$dataPath)){
        stop("No such directory for the 'dataPath':",argList$dataPath ,".\n")
    }

    if(!(argList$species %in% c("human", "mouse"))){
        stop("The parameter 'species' should be one of the c(\"human\", \"mouse\").\n")
    }

    if(!is.numeric(argList$hg.mm.thres)){
        stop("The parameter 'hg.mm.thres' should be a float-point number within [0.5, 1].\n")
    }else if(argList$hg.mm.thres < 0.5 | argList$hg.mm.thres > 1){
        stop("The parameter 'hg.mm.thres' should be within [0.5, 1].\n")
    }
}


#' checkAnnoArguments
#'
#' @param argList A list of arguments passed into 'runScAnnotation".
#'
#' @return NULL
#' @export
#'
checkAnnoArguments <- function(argList){
    if(!dir.exists(argList$dataPath)){
        stop("No such directory for the 'dataPath':",argList$dataPath ,".\n")
    }

    if(!dir.exists(argList$statPath)){
        stop("No such directory for the 'statPath':",argList$statPath ,".\n")
    }

    if(!(argList$species %in% c("human", "mouse"))){
        stop("The parameter 'species' should be one of the c(\"human\", \"mouse\").\n")
    }

    if(!(argList$genome %in% c("hg19", "hg38", "mm10"))){
        stop("The parameter 'genome' should be one of the c(\"hg19\", \"hg38\", \"mm10\").\n")
    }

    if(!(all(argList$anno.filter %in% c("mitochondrial", "ribosome", "dissociation"))) &
       !(is.null(argList$anno.filter))){
        stop("The parameter 'anno.filter' should be some of c(\"mitochondrial\", \"ribosome\", \"dissociation\") or NULL.\n")
    }

    if(!(argList$doublet.method %in% c("cxds", "bcds"))){
        stop("The parameter 'doublet.method' should be one of the c(\"cxds\", \"bcds\").\n")
    }

    if(!(all(argList$coor.names == c("tSNE_1", "tSNE_2")) |
         all(argList$coor.names == c("UMAP_1", "UMAP_2")))){
        stop("The parameter 'coor.names' should be c(\"tSNE_1\", \"tSNE_2\") or c(\"UMAP_1\", \"UMAP_2\").\n")
    }

    if(!(argList$geneSet.method %in% c("average", "GSVA"))){
        stop("The parameter 'geneSet.method' should be one of the c(\"average\", \"GSVA\").\n")
    }
}



#' checkCombArguments
#'
#' @param argList A list of arguments passed into 'runScCombination".
#'
#' @return NULL
#' @export
#'
checkCombArguments <- function(argList){
    if(length(argList$single.savePaths) != length(argList$sampleNames)){
        stop("The length of parameter 'single.savePaths' and 'sampleNames' should be equal.\n")
    }
    if(!(argList$comb.method %in% c("Harmony", "NormalMNN", "SeuratMNN", "Raw", "Regression", "LIGER"))){
        stop("The parameter 'comb.method' should be one of the c(\"Harmony\", \"NormalMNN\", \"SeuratMNN\", \"Raw\", \"Regression\", \"LIGER\").\n")
    }
}
