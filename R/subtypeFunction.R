# Scibet method
# Li, C., Liu, B., Kang, B. et al.
# SciBet as a portable and fast single cell type identifier.
# Nat Commun 11, 1818 (2020). https://doi.org/10.1038/s41467-020-15523-2

#' Entropy_test
#' @param expr Rows should be cells and the last column should be "label".
#' @param k Number of genes to select.
#'
#' @return  'out'  selected gene list(class "character)
#' @export
Entropy_test <- function(expr, k = 1000) {
  labels <- factor(expr$label)
  labels_set <- levels(labels)
  label_total <- matrix(0,ncol(expr)-1,1)
  for(i in labels_set){
    label_TPM <- expr[expr$label==i,][,-ncol(expr)]
    label_mean <- colMeans(label_TPM)
    temp <- t(matrix(label_mean,1))
    colnames(temp) <- i
    label_total <- cbind(label_total, temp)
  }
  label_total <- label_total[,-1]
  #E-test
  log_E <- log2(rowMeans(label_total + 1))
  E_log <- rowMeans(log2(label_total + 1))
  t_scores <- log_E - E_log
  out <- cbind(label_total, t_scores)
  rownames(out) <- colnames(expr)[-ncol(expr)]
  out <- out[order(-out[, "t_scores"]),]
  select <- out[1:k,]
  out <- rownames(select)
  return(out)
}

#' HRG
#' @param expr Rows should be cells and the last column should be "label".
#' @param k Number of genes to select.
#'
#' @return  'out'  selected gene list(class "character)
#' @export
HRG <- function(counts, k = 1000) {
  object.seurat <- CreateSeuratObject(counts)
  all.genes = rownames(object.seurat)
  object.seurat <- ScaleData(object.seurat, features=all.genes,
                             verbose=FALSE)
  object.seurat <- RunPCA(object.seurat, features=all.genes,
                          verbose=FALSE)
  object.seurat <- FindRegionalGenes(object.seurat,dims=1:10,
                                     nfeatures=k, overlap_stop=0.99)
  out <- head(RegionalGenes(object.seurat), k)
  # Replace dashes in Seurat's feature names with underscores
  out <- gsub('[-]', '_', out)
  return(out)
}


#' SelectGenes
#' @name SelectGenes
#' @param expr Rows should be cells and the last column should be "label".
#'
#' @return geneset A user-defined set of genes for prediction.
#' @export
#'
SelectGenes <- function(expr, method="Entropy", k = 2000){
  if(method == "Entropy"){
    geneset <- Entropy_test(expr, k)
  }
  if(method == "HRG"){
    geneset <- HRG(data.frame(t(expr[, 1:dim(expr)[2]-1])), k)
  }
  return(geneset)
}


#' SelectCells
#' @name SelectCells
#' @param object Seurat object of training set.
#' @return representative.index A user-defined set of cells for training.
#' @export
#' @import org.Hs.eg.db
#' @importFrom garnett check_markers
SelectCells <- function(object,
                        marker_file_path,
                        cutoff = 0.75,
                        database = org.Hs.eg.db){
  marker_check <- check_markers(object, marker_file_path,
                                db=database,
                                cds_gene_id_type = "SYMBOL",
                                marker_file_gene_id_type = "SYMBOL")
  # print(plot_markers(marker_check))

  # select representative samples for training
  set.seed(unclass(Sys.time()))
  training_sample <- select_fine_samples(cds = object,
                                         marker_file = marker_file_path,
                                         db=database,
                                         cds_gene_id_type = "SYMBOL",
                                         cutoff = cutoff,
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL",
                                         return_initial_assign = TRUE)
  representative.index <- which(training_sample$Assignment != "None")
  return(list(index = representative.index,
              markers = marker_check))
}

#' Train
#' @name Train
#' @usage Train(expr, geneset, k=1000)
#' @param expr The reference expression dataframeThe expression dataframe, with rows being cells, and columns being genes. The last column should be "label".
#' @param geneset A user-defined set of genes for prediction.
#'
#' @return reference model(cellType-genes prob-matrix)
#' @export
#'
Train <- function(expr, geneset=NULL){
    if(is.null(geneset)){
        print("All scRNA-seq features...")
        path <- paste0(system.file("txt", package = "scCancer"), "/single_cell_features.tsv")
        gene.list <- read.table(path)$V2
        index_remain <- which(colnames(expr) %in% intersect(colnames(expr), gene.list))
        geneset <- colnames(expr)[index_remain]
    }
    labels <- factor(expr$label)
    label.level <- levels(labels)
    expr_select <- expr[,geneset]
    label_total <- matrix(0, length(geneset), 1)
    lambda_all <- matrix(0, length(geneset), 1)

    # Calculate zero-ratio lambda for every possible cell type.
    dropout.ratio <- c()
    for(j in 1:dim(expr_select)[2]){
        dropout.ratio <- c(dropout.ratio, length(which(expr_select[,j] == 0)) / dim(expr_select)[1])
    }

    # all possible labels
    for(i in label.level){
        label_TPM <- expr_select[labels == i,]
        # Calculate average expression.
        label_mean <- colSums(log2(label_TPM + 1)) # log normalization
        temp <- t(matrix(label_mean, 1))
        colnames(temp) <- i
        label_total <- cbind(label_total, temp)
        temp <- t(matrix(dropout.ratio, 1))
        colnames(temp) <- i
        lambda_all <- cbind(lambda_all, temp)
    }
    label_total <- label_total[,-1]
    rownames(label_total) <- geneset
    lambda_all <- lambda_all[,-1]
    rownames(lambda_all) <- geneset
    # calculate probability and log transform
    label_t <- t(label_total)
    # log2[(X + 1) / sum(X + 1)]
    prob <- log2(label_t + 1) - log2(rowSums(label_t) + length(geneset))
    prob <- t(prob)
    return(cbind(prob, lambda_all))
}

#' pro.core
#' @name pro.core
#' @usage pro.core(scibet.core)
#' @param scibet.core A SciBet core
#' @return A processed SciBet core
#' @export
pro.core <- function(scibet.core){
    cell.type <- unname(unlist(scibet.core[,1]))
    scibet.core <- as.data.frame(t(scibet.core[,-1]))
    colnames(scibet.core) <- cell.type
    return(as.matrix(scibet.core))
}

Dropout_Sampling <- function(lambda){
    nondropout <- as.matrix(1 - lambda)
    z <- c()
    for(i in 1:nrow(nondropout)){
        pos.prob <- nondropout[i,1]
        newz <- sample(0:1, ncol(nondropout), replace = TRUE,
                       prob = c(1 - pos.prob, pos.prob))
        z <- rbind(z, newz)
    }
    rownames(z) <- rownames(nondropout)
    colnames(z) <- colnames(nondropout)
    return(as.matrix(z))
}

#' MLEstimate
#' @name MLEstimate
#' @param test Rows should be cells and columns should be genes.
#' @param prob trained scibet model  Rows should be genes
#' and columns should be cell types.The genes must be matched with test_r
#' @return  'cellType' and 'cell*label' likelihood matrix
MLEstimate <- function(test,
                       prob,
                       lambda = NULL,
                       weighted.markers = NULL,
                       dropout.modeling = FALSE){
    #\sigma yi*logpi (log(p1p2) = logp1 + logp2, log(p^y) = ylogp)
    if (length(weighted.markers) > 0){
        weighted.markers <- weighted.markers[!is.na(weighted.markers)]
        indexes <- which(rownames(prob) %in% names(weighted.markers))
        for (index in indexes){
            prob[index,] <- weighted.markers[rownames(prob)[index]] * prob[index,]
        }
        # prob[which(rownames(prob) %in% markers),] <- 2 * prob[which(rownames(prob) %in% markers),]
        # prob <- apply(prob, 2, function(p){return(p / sum(p))})
    }
    if (dropout.modeling){
        # z~Bernoulli(1-lambda)
        # z <- Dropout_Sampling(lambda)
        # likelihoods <- as.matrix(test) %*% (as.matrix(prob) * z)
        # Expectation of log-likelihoods.
        likelihoods <- as.matrix(test) %*% (as.matrix(prob) * as.matrix(1 - lambda))
    }
    else{
        likelihoods <- as.matrix(test) %*% as.matrix(prob)
    }
    cellType <- c()
    for(i in 1:nrow(likelihoods)){
        index <- which.max(likelihoods[i,])
        cellType <- c(cellType,colnames(likelihoods)[index])
    }
    return(list(cellType = cellType,
                likelihoods = likelihoods))
}

#' Test
#' @name Test
#' @usage Test(prob, test, result="list")
#' @param prob A matrix generated by Train(cellType-genes prob-matrix)
#' @param test testset with rows being cells, and columns being genes.
#' @return cellType prediction result
#' @export
Test <- function(prob, lambda, test_set,
                 weighted.markers = NULL,
                 dropout.modeling = FALSE,
                 average.expr = NULL){
    # metacell
    if(!is.null(average.expr)){
        test <- average.expr
    }
    else{
        test <- test_set[,-ncol(test_set)]
    }
    genes <- rownames(prob)
    common.genes <- intersect(genes, colnames(test))
    test.normalized <- log1p(as.matrix(test[,common.genes])) / log(2)
    result <- MLEstimate(test = test.normalized,
                          prob = prob[common.genes, ],
                          lambda = lambda[common.genes, ],
                          weighted.markers = weighted.markers,
                          dropout.modeling = dropout.modeling)
    predict <- result[["cellType"]]
    likelihoods <- result[["likelihoods"]]
    # name by clusters
    if(!is.null(average.expr)){
        names(predict) <- seq(from = 0, to = length(predict) - 1)
    }
    # calculate accuracy if possible
    else{
        if(!is.null(test_set$label)){
            correct <- sum(test_set$label == predict)
            message("Accuracy: ", correct / length(predict))
        }
    }
    return(list(predict = predict, likelihoods = likelihoods))
}

#' CrossTest
#' @name CrossTest
#' @usage CrossTest(prob, test_set)
#' @param prob A matrix generated by Train(cellType-genes prob-matrix)
#' @param test_set testset with rows being cells, and columns being genes.
#' @return cellType prediction result
#' @export
CrossTest <- function(prob, test_set){
    # test <- test_set[,-ncol(test_set)]
    genes <- rownames(prob)
    common.genes <- intersect(genes, colnames(test_set))
    testa <- log1p(as.matrix(test_set[,common.genes])) / log(2)
    predict <- MLEstimate(testa, prob[common.genes, ])
    return(predict)
}


#' MarkerScore
#' @name MarkerScore
#' @export
#' @import org.Hs.eg.db
#' @importFrom DESeq2 estimateSizeFactors
#' @importFrom garnett check_markers
MarkerScore <- function(test_set,
                        marker_file_path,
                        cutoff = 0.3,
                        database = org.Hs.eg.db,
                        metacell = FALSE){
    # check markers on test set, set unknown labels
    set.seed(unclass(Sys.time()))
    object <- visualization_pipeline(test_set, metacell = metacell)[["object"]]
    average.expr <- NULL
    clustering <- NULL
    if(metacell){
        average.expr <- data.frame(t(AverageExpression(object)[[1]]))
        clustering <- as.character(object@meta.data[["RNA_snn_res.100"]])
    }
    object <- as.CellDataSet(object)
    object <- estimateSizeFactors(object)
    # adjust cutoff parameter
    suppressWarnings(marker_check <- check_markers(object,
                                  marker_file_path,
                                  db = database,
                                  cds_gene_id_type = "SYMBOL",
                                  marker_file_gene_id_type = "SYMBOL"))
    weighted.markers <- 1 + log(1 + marker_check$marker_score)
    names(weighted.markers) <- marker_check$marker_gene
    suppressWarnings(result <- select_fine_samples(cds = object,
                                  marker_file = marker_file_path,
                                  db=database,
                                  cds_gene_id_type = "SYMBOL",
                                  cutoff = cutoff,
                                  num_unknown = 50,
                                  marker_file_gene_id_type = "SYMBOL",
                                  return_initial_assign = TRUE))
    unknown.index <- which(result$Assignment == "None")
    return(list(average = average.expr,
                clustering = clustering,
                unknown = unknown.index,
                markers = weighted.markers))
}

AssignUnknown <- function(test_set, predict, unknown.index){
    # calculate accuracy after "unknown"s are excluded
    correct <- 0
    predict.unknown <- predict
    for (i in 1:length(predict)){
        if (i %in% unknown.index){
            predict.unknown[i] <- "unknown"
            next
        }
        if(is.null(test_set)){
            next
        }
        if (test_set$label[i] == predict[i]){
            correct <- correct + 1
        }
    }
    accuracy <- correct / (length(predict.unknown) - length(unknown.index))
    return(list(predict.unknown = predict.unknown,
                accuracy = accuracy))
}
