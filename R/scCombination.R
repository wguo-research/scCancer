
#' runScCombination
#'
#' Perform multi-samples analyses.
#'
#' @param single.savePaths A vecotr of paths containing the results files of step 'runScAnnotation' for each sample.
#' @param sampleNames A vector of labels for all samples.
#' @param combName A label for the combined samples.
#' @param comb.method The method to combine samples. The default is "NormalMNN". "Harmony", "NormalMNN", "SeuratMNN", "Raw", "Regression" and "LIGER" are optional.
#' @param harmony.theta The parameter 'theta' of function "RunHarmony" in the harmony package.
#' @param harmony.lambda The parameter 'lambda' of function "RunHarmony" in the harmony package.
#' @param harmony.sigma The parameter 'sigma' of function "RunHarmony" in the harmony package.
#' @param sample.colors The colors used for samples. The default is NULL, and the pre-set colors will be used.
#' @inheritParams runScAnnotation
#'
#' @return A results list with all useful objects used in the function.
#' @export
#'
#' @import harmony rliger
#'
runScCombination <- function(single.savePaths, sampleNames, savePath, combName,
                             authorName = NULL,
                             comb.method = "NormalMNN",
                             harmony.theta = NULL,
                             harmony.lambda = NULL,
                             harmony.sigma = 0.1,
                             vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent"),
                             pc.use = 30,
                             resolution = 0.8,
                             clusterStashName = "comb.cluster",
                             show.features = NULL, bool.add.features = T,
                             bool.runDiffExpr = T,
                             n.markers = 5,
                             sample.colors = NULL,
                             species = "human",
                             genome = "hg19",
                             hg.mm.mix = F,
                             bool.runCellClassify = T,
                             ct.templates = NULL,
                             coor.names = c("tSNE_1", "tSNE_2"),
                             bool.runMalignancy = T,
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
                             genReport = T){

    message("[", Sys.time(), "] START: RUN ScCombination")
    results <- as.list(environment())
    checkCombArguments(results)

    if(species == "mouse" & genome == "hg19"){
        genome <- "mm10"
    }

    if(!dir.exists(file.path(savePath, "figures/"))){
        dir.create(file.path(savePath, "figures/"), recursive = T)
    }
    suppressWarnings( savePath <- normalizePath(savePath, "/") )
    results[["savePath"]] <- savePath


    message("[", Sys.time(), "] -----: sample data combination")
    expr.list <- list()
    sample.ident <- c()
    for(i in 1:length(sampleNames)){
        sampleName <- sampleNames[i]
        cur.path <- single.savePaths[i]
        print(sampleName)
        expr.list[[sampleName]] <- readRDS(paste0(cur.path, "/expr.RDS"))
        sample.ident <- c(sample.ident, rep(sampleName, dim(expr.list[[sampleName]])[2]))
    }
    sample.ident <- as.factor(sample.ident)

    bool.plotHVG = T
    if(comb.method == "SeuratMNN"){
        message("[", Sys.time(), "] -----: combine data by Seurat MNN")
        suppressWarnings( expr.anchors <- FindIntegrationAnchors(object.list = expr.list,
                                                                 dims = 1:pc.use) )
        expr <- IntegrateData(anchorset = expr.anchors,
                              dims = 1:pc.use, verbose = F)
        expr <- ScaleData(expr, verbose = FALSE)
        DefaultAssay(expr) <- "integrated"
        expr[["sample.ident"]] <- sample.ident
        bool.plotHVG = F

        saveRDS(expr.anchors@anchors, file = file.path(savePath, "anchors.RDS"))

    }else if(comb.method == "Raw"){
        message("[", Sys.time(), "] -----: combine raw matrix data")
        suppressWarnings( expr <- merge(expr.list[[1]], expr.list[2:length(expr.list)]) )
        expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000, verbose = F)
        expr <- ScaleData(object = expr, vars.to.regress = vars.to.regress, verbose = F)
        expr[["sample.ident"]] <- sample.ident

    }else if(comb.method == "Regression"){
        message("[", Sys.time(), "] -----: combine data and regress out sample source")
        suppressWarnings( expr <- merge(expr.list[[1]], expr.list[2:length(expr.list)]) )
        expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000, verbose = F)
        expr[["sample.ident"]] <- sample.ident
        expr <- ScaleData(object = expr,
                          vars.to.regress = c("sample.ident", vars.to.regress),
                          verbose = F)

    }else if(comb.method == "Harmony"){
        message("[", Sys.time(), "] -----: combine data by Harmony")

        items <- unique(unlist(lapply(names(expr.list), function(x){
            grep("^GS__", names(expr.list[[x]]@meta.data), value = T)
        })))
        items <- c("doublet.score", "Cell.Type", "Malign.score",
                   "Malign.type", "CellCycle.score", "Stemness.score", items)

        ju.mat <- sapply(names(expr.list), function(x){
            !(items %in% names(expr.list[[x]]@meta.data))
        })
        comb.metadata <- lapply(items[rowSums(ju.mat) == 0], function(x){
            tmp <- do.call(c, lapply(names(expr.list), function(y){
                expr.list[[y]]@meta.data[[x]]
            }))
        })
        names(comb.metadata) <- items[rowSums(ju.mat) == 0]
        comb.metadata <- data.frame(comb.metadata)

        share.genes <- Reduce(intersect,  lapply(expr.list, rownames))
        for(s.name in names(expr.list)){
            expr.list[[s.name]] <- GetAssayData(expr.list[[s.name]], slot = "counts")[share.genes, ]
        }
        comb.data <- do.call(cbind, expr.list)
        rm(expr.list)

        expr <- CreateSeuratObject(counts = comb.data,  min.cells = 5) %>%
            Seurat::NormalizeData(verbose = FALSE) %>%
            FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>%
            ScaleData(verbose = FALSE) %>%
            RunPCA(pc.genes = expr@var.genes, verbose = FALSE)
        expr[["sample.ident"]] <- sample.ident
        expr <- expr %>% RunHarmony("sample.ident", plot_convergence = TRUE,
                                    theta = harmony.theta,
                                    lambad = harmony.lambda,
                                    sigma = harmony.sigma,
                                    verbose = F)

        expr@meta.data <- cbind(expr@meta.data, comb.metadata)

        bool.plotHVG <- F

    }else if(comb.method == "LIGER"){
        message("[", Sys.time(), "] -----: combine data by LIGER")

        items <- unique(unlist(lapply(names(expr.list), function(x){
            grep("^GS__", names(expr.list[[x]]@meta.data), value = T)
        })))
        items <- c("doublet.score", "Cell.Type", "Malign.score",
                   "Malign.type", "CellCycle.score", "Stemness.score", items)

        ju.mat <- sapply(names(expr.list), function(x){
            !(items %in% names(expr.list[[x]]@meta.data))
        })
        comb.metadata <- lapply(items[rowSums(ju.mat) == 0], function(x){
            tmp <- do.call(c, lapply(names(expr.list), function(y){
                expr.list[[y]]@meta.data[[x]]
            }))
        })
        names(comb.metadata) <- items[rowSums(ju.mat) == 0]
        comb.metadata <- data.frame(comb.metadata)

        for(e.i in 1:length(expr.list)){
            s.name <- names(expr.list)[e.i]
            expr.list[[s.name]] <- RenameCells(expr.list[[s.name]],
                                               new.names = paste0(colnames(expr.list[[s.name]]), "-", e.i))
            expr.list[[s.name]] <- GetAssayData(expr.list[[s.name]], slot = "counts")
        }
        expr = createLiger(expr.list)
        expr = normalize(expr)
        expr = selectGenes(expr, var.thresh = 0.1)
        expr = scaleNotCenter(expr)

        expr = optimizeALS(expr, k = 20)
        expr = quantileAlignSNF(expr)
        expr = runTSNE(expr)
        expr = ligerToSeurat(expr, use.liger.genes = T)

        expr = ScaleData(expr, verbose = FALSE)
        expr[["sample.ident"]] <- sample.ident
        expr@reductions$inmf@assay.used <- "RNA"

        expr@meta.data <- cbind(expr@meta.data, comb.metadata)

        bool.plotHVG = F

    }else if(comb.method == "NormalMNN"){
        message("[", Sys.time(), "] -----: combine data by normal cell MNN")
        suppressWarnings( expr.anchors <- FindIntegrationAnchors(object.list = expr.list,
                                                                 dims = 1:pc.use) )
        anchors <- expr.anchors@anchors

        anchors$cellType1 <- "NULL"
        anchors$cellType2 <- "NULL"
        anchors$malignType1 <- "NULL"
        anchors$malignType2 <- "NULL"
        anchors$malignScore1 <- -1
        anchors$malignScore2 <- -1
        for(oi in expr.anchors@reference.objects){
            cur.ix <- which(anchors$dataset1 == oi)
            anchors$cellType1[cur.ix] <- expr.list[[oi]]@meta.data$Cell.Type[anchors$cell1[cur.ix]]
            anchors$malignType1[cur.ix] <- expr.list[[oi]]@meta.data$Malign.type[anchors$cell1[cur.ix]]
            anchors$malignScore1[cur.ix] <- expr.list[[oi]]@meta.data$Malign.score[anchors$cell1[cur.ix]]

            cur.ix <- which(anchors$dataset2 == oi)
            anchors$cellType2[cur.ix] <- expr.list[[oi]]@meta.data$Cell.Type[anchors$cell2[cur.ix]]
            anchors$malignType2[cur.ix] <- expr.list[[oi]]@meta.data$Malign.type[anchors$cell2[cur.ix]]
            anchors$malignScore2[cur.ix] <- expr.list[[oi]]@meta.data$Malign.score[anchors$cell2[cur.ix]]
        }

        anchors.new <- subset(anchors, cellType1 != "Epithelial" & cellType1 != "Unknown" & cellType2 != "Epithelial" & cellType2 != "Unknown")
        if(dim(anchors)[1] == 0){
            anchors.new <- anchors
            cat("- Warning in 'runScCombination': Cannot find the nomral cell anchors, and use initial anchors instead.\n")
        }
        expr.anchors@anchors <- anchors.new

        expr <- IntegrateData(anchorset = expr.anchors,
                              dims = 1:pc.use, verbose = F)
        expr <- ScaleData(expr, verbose = FALSE)
        DefaultAssay(expr) <- "integrated"
        expr[["sample.ident"]] <- sample.ident
        bool.plotHVG = F

        saveRDS(anchors.new, file = file.path(savePath, "anchors.RDS"))
    }
    results[["bool.plotHVG"]] <- bool.plotHVG

    ## --------- seurat ---------
    t.results <- runSeurat(
        expr = expr,
        savePath = savePath,
        pc.use = pc.use,
        resolution = resolution,
        clusterStashName = clusterStashName,
        bool.runDiffExpr = bool.runDiffExpr,
        comb.method = comb.method
    )
    expr = t.results$expr
    cell.annotation = t.results$cell.annotation
    results[["diff.expr.genes"]] = t.results$diff.expr.genes
    rm(t.results)
    gc()

    for(item in c("doublet.score", "Cell.Type", "Malign.score",
                  "Malign.type", "CellCycle.score", "Stemness.score")){
        if(item %in% names(expr@meta.data)){
            cell.annotation[[item]] <- expr@meta.data[[item]]
        }
    }
    for(item in grep("^GS__", names(expr@meta.data), value = T)){
        cell.annotation[[item]] <- expr@meta.data[[item]]
    }

    results[["seurat.plots"]] <- plotSeurat(
        expr = expr,
        cell.annotation = cell.annotation,
        show.features = show.features,
        bool.add.features = bool.add.features,
        coor.names = coor.names,
        bool.plotHVG = bool.plotHVG,

        bool.runDiffExpr = bool.runDiffExpr,
        diff.expr.genes = results[["diff.expr.genes"]],
        n.markers = n.markers,

        species = species,
        savePath = savePath
    )

    results[["DEplot.height"]] <- 0.5 + 0.1 * n.markers * length(unique(cell.annotation$Cluster))
    results[["markersPlot.height"]] <- 2 * ceiling(length(results[["seurat.plots"]]$ps.markers) / 4)


    ## --------- sample source ---------
    message("[", Sys.time(), "] -----: plot sample source")
    cell.annotation$sample <- expr@meta.data$sample.ident
    if(is.null(sample.colors)){
        sample.colors <- getDefaultColors(n = length(unique(cell.annotation$sample)),
                                          type = 2)
    }

    if(setequal(sampleNames, unique(cell.annotation$sample))){
        cell.annotation$sample <- factor(cell.annotation$sample, levels = sampleNames)
    }else{
        cell.annotation$sample <- factor(cell.annotation$sample)
    }
    p.sample <- pointDRPlot(cell.annotation, value = "sample",
                            coor.names = coor.names,
                            colors = sample.colors,
                            point.type = 2,
                            legend.position = "right",
                            legend.title = "Sample")
    p.bar.sample <- clusterBarPlot(cell.annotation = cell.annotation,
                                   cell.colors = sample.colors,
                                   sel.col = "sample",
                                   legend.position = "bottom",
                                   legend.title = "Sample")

    ggsave(filename = file.path(savePath, "figures/sampleSource-point.png"),
           p.sample, width = 7, height = 5, dpi = 300)
    ggsave(filename = file.path(savePath, "figures/sampleSource-bar.png"),
           p.bar.sample, width = 6, height = 3, dpi = 300)
    results[["p.sample"]] <- p.sample
    results[["p.bar.sample"]] <- p.bar.sample


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
    }


    ## --------- malignancy ---------
    if(bool.runMalignancy){
        if(!(all(c("Malign.score", "Malign.type") %in% names(cell.annotation)))){
            message("[", Sys.time(), "] -----: cells malignancy annotation through inferCNV")
            for(i in 1:length(sampleNames)){
                cur.manifest <- read.table(paste0(single.savePaths[i], "/geneManifest.txt"),
                                           header = T, sep = "\t", stringsAsFactors = F)
                if(i == 1){
                    gene.manifest <- cur.manifest
                }else{
                    new.genes <- subset(cur.manifest, !(EnsemblID %in% gene.manifest$EnsemblID))
                    gene.manifest <- rbind(gene.manifest, new.genes)
                }
            }
            # rownames(gene.manifest) <- gene.manifest$EnsemblID
            rownames(gene.manifest) <- gene.manifest$Symbol
            # Run inferCNV
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
            results[["bimodal.pvalue"]] <- t.results$bimodal.pvalue
            results[["malign.plot"]] <- t.results$p.results
            rm(t.results)
        }else{
            message("[", Sys.time(), "] -----: cells malignancy combination")
            results[["malign.plot"]] <- plotMalignancy(cell.annotation = cell.annotation,
                                                       malignancy.method = "",
                                                       coor.names = coor.names,
                                                       savePath = savePath)
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
        if(!("CellCycle.score" %in% names(cell.annotation))){
            CellCycle.score <- runCellCycle(expr, species = species)
            cell.annotation$CellCycle.score <- CellCycle.score
            expr[["CellCycle.score"]] <- CellCycle.score
        }else{
            message("[", Sys.time(), "] -----: cell cycle score combination")
        }

        # CellCycle.score <- runCellCycle(expr, species = species)
        # cell.annotation$CellCycle.score <- CellCycle.score
        # expr[["CellCycle.score"]] <- CellCycle.score

        results[["cellCycle.plot"]] <-
            pointDRPlot(cell.annotation,
                        sel.clusters = sel.clusters,
                        value = "CellCycle.score",
                        coor.names = coor.names,
                        colors = c("white", "#009b45"),
                        discrete = F,
                        legend.position = "right",
                        legend.title = "Cell cycle score")
        ggsave(filename = file.path(savePath, "figures/cellCycle-point.png"),
               results[["cellCycle.plot"]], width = 5, height = 4, dpi = 300)
    }


    ## --------- stemness ---------
    if(bool.runStemness){
        if(!("Stemness.score" %in% names(cell.annotation))){
            stem.scores <- runStemness(X = GetAssayData(object = expr, slot = "scale.data"), species = species)
            cell.annotation[["Stemness.score"]] <- stem.scores
            expr[["Stemness.score"]] <- stem.scores
        }else{
            message("[", Sys.time(), "] -----: stemness score combination")
        }

        results[["stemness.plot"]] <-
            pointDRPlot(cell.annotation,
                        sel.clusters = sel.clusters,
                        value = "Stemness.score",
                        coor.names = coor.names,
                        colors = c("white", "#ff9000"),
                        discrete = F,
                        legend.position = "right",
                        legend.title = "Stemness")
        ggsave(filename = file.path(savePath, "figures/stemness-point.png"),
               results[["stemness.plot"]], width = 5, height = 4, dpi = 300)
    }


    ## --------- gene sets ----------
    if(bool.runGeneSets){
        if(is.null(geneSets)){
            geneSets <- getDefaultGeneSets(species = species)
        }
        if(geneSet.method == "GSVA" | !all(paste0("GS__", names(geneSets)) %in% names(cell.annotation))){
            t.scores <- runGeneSets(expr = expr, geneSets = geneSets, method = geneSet.method)
            if(!is.null(t.scores)){
                cell.annotation <- cbind(cell.annotation, t.scores)
            }
        }else{
            message("[", Sys.time(), "] -----: gene set signatures combination")
            t.scores <- cell.annotation[, paste0("GS__", names(geneSets))]
        }

        if(!is.null(t.scores)){
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
                                                           clusterStashName = clusterStashName,
                                                           savePath = savePath)
        results[["exprProgram.plot"]] <- plotExprProgram(H = results[["exprProgram.results"]]$H,
                                                         cell.annotation,
                                                         sel.clusters = sel.clusters,
                                                         savePath = savePath)
        results[["exprProgPlot.height"]] <- 0.5 + 0.11 * dim(results[["exprProgram.results"]]$H)[1]
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
            knit(system.file("rmd", "main-scAnnoComb.Rmd", package = "scCancer"),
                 file.path(savePath,'report-scAnnoComb.md'), quiet = T)
        )
        markdownToHTML(file.path(savePath,'report-scAnnoComb.md'),
                       file.path(savePath, 'report-scAnnoComb.html'))
    }

    message("[", Sys.time(), "] END: Finish ScCombination\n\n")

    return(results)
}

