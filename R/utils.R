
get10Xpath <- function (samplePath, raw.data = F){
    prefix <- ifelse(raw.data, 'raw', 'filtered')
    cur.path <- paste0(samplePath, '/')
    res.path <- paste0(cur.path, prefix, '_feature_bc_matrix')
    if (!dir.exists(res.path)){
        res.path <- paste0(cur.path, prefix, '_gene_bc_matrices/hg19/')
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
    version <- "Cell Ranger V2"
    if(grepl("feature_bc_matrix", data.path)){
        version <- "Cell Ranger V3"
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



Read10Xdata <- function (data.dir = NULL, gene.column = 2) {
    full.data <- list()
    for (i in seq_along(data.dir)) {
        run <- data.dir[i]
        if (!dir.exists(run)) {
            stop("Directory provided does not exist.\n")
        }
        if (!grepl("\\/$", run)) {
            run <- paste(run, "/", sep = "")
        }
        barcode.loc <- paste0(run, "barcodes.tsv")
        gene.loc <- paste0(run, "genes.tsv")
        features.loc <- paste0(run, "features.tsv.gz")
        matrix.loc <- paste0(run, "matrix.mtx")
        # Flag to indicate if this data is from CellRanger >= 3.0
        pre_ver_3 <- file.exists(gene.loc)
        if (!pre_ver_3) {
            addgz <- function(s) {
                return(paste0(s, ".gz"))
            }
            barcode.loc <- addgz(s = barcode.loc)
            matrix.loc <- addgz(s = matrix.loc)
        }
        if (!file.exists(barcode.loc)) {
            stop("Barcode file missing.\n")
        }
        if (!pre_ver_3 && !file.exists(features.loc)) {
            stop("Gene name or features file missing.\n")
        }
        if (!file.exists(matrix.loc)) {
            stop("Expression matrix file missing.\n")
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
        rownames(x = data) <-
            make.unique(names = feature.names[, gene.column])
        # In cell ranger 3.0, a third column specifying the type of data was added
        # and we will return each type of data as a separate matrix
        if (ncol(x = feature.names) > 2) {
            # if (length(x = full.data) == 0) {
            #     message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
            # }
            data_types <- factor(x = feature.names$V3)
            lvls <- levels(x = data_types)
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
        # Fix for Issue #913
        list_of_data[[j]] <-
            as(object = list_of_data[[j]], Class = "dgCMatrix")
    }
    names(x = list_of_data) <- names(x = full.data[[1]])
    # If multiple features, will return a list, otherwise
    # a matrix.
    if (length(x = list_of_data) == 1) {
        return(list_of_data[[1]])
    } else {
        return(list_of_data)
    }
}



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
            "Myeloblasts" = c("LYZ"),
            "Endothelial" = c("PLVAP"),
            "Myofibroblast" = c("ACTA2", "S100A4"),
            "Epithelial" = c("EPCAM", "KRT8"))
    }else if(species == "mouse"){
        feature.def <- list(
            "T cell" = c("Ptprc", "Cd3d", "Cd4", "Cd8a", "Cd8b"),
            "B cell" = c("Cd79a"),
            "NK cell" = c("Nkg7"),
            "Myeloblasts" = c("Lyz1", "Lyz2"),
            "Endothelial" = c("Plvap"),
            "Myofibroblast" = c("Acta2", "S100a4"),
            "Epithelial" = c("Epcam", "Krt8"))
    }

    return(feature.def)
}



getDefaultColors <- function(n = NULL){
    colors <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb",
                "#d7652d", "#7cd5c8", "#c49a3f", "#507d41", "#5d8d9c",
                "#90353b", "#674c2a", "#1B9E77", "#c5383c", "#0081d1",
                "#ffd900", "#502e71", "#c8b693", "#aed688", "#f6a97a",
                "#c6a5cc", "#798234", "#6b42c8", "#cf4c8b", "#666666")
    if(!is.null(n)){
        if(n <= 25){
            colors <- colors[1:n]
        }else{
            step <- 16777200 %/% n - 2
            colors <- paste0("#", as.hexmode(seq(from = sample(1:step, 1), by = step, length.out = n)))
        }
    }
    return(colors)
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

