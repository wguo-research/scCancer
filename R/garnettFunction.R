# Garnett method
# Hannah A. Pliner, Jay Shendure & Cole Trapnell (2019).
# Supervised classification enables rapid annotation of cell atlases. Nature Methods

# Function select_fine_samples is modified for:
# 1. selecting core dataset for model training(at least 1 marker expressing)
# 2. assigning unknown label(no marker expressing)

# Part1. garnett parser for markers files
# --------------------------------------------------------------------
Lexer <- R6::R6Class(
  "Lexer",
  public = list(
    tokens = c("KEYWORD", "NAME", "NEWLINE", "NUM",
               "SIMP_KEY", "COM"),
    literals = c(">", ":", ","),
    t_COM = function(re='#.*',t) {}, # Comment character
    t_KEYWORD =
      c('expressed below|expressed above|expressed between|subtype of'),
    t_SIMP_KEY = c('celltype|expressed|not expressed|references'),
    t_NAME =
      '[a-zA-Z0-9_+/\\-\\.|=`~\\*&<^%?@!$();:]*[a-zA-Z][a-zA-Z0-9_+/\\-\\.|=`~\\*&<^%?@!$();]*',
    t_NUM = '([0-9]*\\.[0-9]+)|([0-9]+)',
    t_ignore = " \t",
    t_NEWLINE = "\n|\r\n",
    t_error = function(t) {
      cat(sprintf(" Marker file error. Illegal character '%s'\n", t$value[1]))
      t$lexer$skip(1)
      return(t)
    }
  )
)

#' Parsing the Garnett marker file
#'
#' Garnett uses a marker file to allow users to specify cell type definitions.
#' While the marker file is designed to be easy to construct and human-readable,
#' it is parsed by Garnett automatically, and so it needs to follow certain
#' formatting constraints.
#'
#' The following describes the constraints necessary in the input to the
#' \code{marker_file} argument of \code{\link{train_cell_classifier}} and
#' \code{\link{check_markers}}.
#'
#' @section Elements of a cell type description:
#' The basic structure of the Garnett marker file is a series of entries, each
#' describing elements of a cell type. After the cell name, each additional
#' line will be a descriptor, which begins with a keyword, followed by a colon
#' (':'). After the colon, a series of specifications can be added, separated
#' by commas (','). Descriptors may spill onto following lines so long as you
#' do not split a specification across multiple lines (i.e. if breaking up a
#' long descriptor across multiple lines, all but the last line should end with
#' a comma). Each new descriptor should begin on a new line. A generic cell
#' type entry looks like this:
#'
#' ```
#' > cell type name
#' descriptor: spec1, spec2,
#' spec3, spec4
#' descriptor2: spec1
#' ```
#'
#' The following are the potential descriptors:
#' \describe{
#'   \item{cell name}{\strong{Required} Each cell type must have a unique name,
#'   and the name should head the cell type description. To indicate a new cell
#'   type, use the \code{>} symbol, followed by the cell name, followed by a
#'   new line. For example, \code{> T cell}.}
#'   \item{expressed:}{\strong{Required} After the cell name, the minimal
#'   requirement for each cell type is the name of a single marker gene. The
#'   line in the marker file will begin with \code{expressed:}, followed by one
#'   or more gene names separated by commas. The last gene name of the
#'   descriptor is not followed by a comma. Gene IDs can be of any type
#'   (ENSEMBL, SYMBOL, etc.) that is present in the Bioconductor
#'   \code{\link[AnnotationDbi]{AnnotationDb-class}} package for your species.
#'   (See available packages on the
#'   \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor website}).
#'   For example, for human, use \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}. To
#'   see available gene ID types, you can run \code{columns(db)}. You will
#'   specify which gene ID type you used when calling
#'   \code{\link{train_cell_classifier}}.} If your species does not have an
#'   annotation dataset of type \code{\link[AnnotationDbi]{AnnotationDb-class}},
#'   you can set \code{db = 'none'}, however Garnett will then not convert gene
#'   ID types, so CDS and marker file gene ID types need to be the same.
#'   \item{not expressed:}{In addition to
#'   specifying genes that the cell type should express, you can also specify
#'   genes that your cell type should not express. Details on specifying genes
#'   are the same as for \code{expressed:}.}
#'   \item{subtype of:}{When present, this descriptor specifies that a cell
#'   type is a subtype of another cell type that is also described in the
#'   marker file. A biological example would be a CD4 T cell being a subtype of
#'   a T cell. This descriptor causes the cell type to be classified on a
#'   separate sub-level of the classification hierarchy, after the
#'   classification of its parent type is done (i.e. first T cells are
#'   discriminated from other cell types, then the T cells are subclassified
#'   using any cell types with the descriptor \code{subtype of: T cell}).
#'   \code{subtype of:} can only include a single specification, and the
#'   specification must be the exact name of another cell type specified in
#'   this marker file.}
#'   \item{references:}{This descriptor is not required, but is highly
#'   recommended. The specifications for this descriptor should be links/DOIs
#'   documenting how you chose your marker genes. While these specifications
#'   will not influence cell type classification, they will be packaged with
#'   the built classifier so that future users of the classifier can trace the
#'   origins of the markers/ }
#'   \item{*meta data:}{This wildcard descriptor allows you to specify any
#'   other property of a cell type that you wish to specify. The keyword will
#'   be the name of the column in your \code{pData} (meta data) table that you
#'   wish to specify, and the specifications will be a list of acceptable
#'   values for that meta data. An example use of this would be
#'   \code{tissue: liver, kidney}, which would specify that training cells for
#'   this cell type must have "liver" or "kidney" as their entry in the
#'   "tissue" column of the \code{pData} table.}
#'   \item{expressed below:}{While we recommend that you use \code{expressed:}
#'   and \code{not expressed:} to specify the cell type's marker genes, because
#'   these terms utilize the entirety of Garnett's built-in normalization and
#'   standardization, you can also specify expression using the following
#'   logical descriptors
#'   \code{expressed below:, expressed above:, expressed between:}.
#'   Note that no normalization occurs with these descriptors; they are used as
#'   logical gates only. To specify \code{expressed below:}, use the gene name,
#'   followed by a space, followed by a number. This will only allow training
#'   cells that have this gene expressed below the given value \strong{in the
#'   units of the expression matrix provided}. For example,
#'   \code{expressed below: MYOD1 7, MYH3 2}.}
#'   \item{expressed above:}{Similar to \code{expressed below:}, but will only
#'   allow training cells expressing the given gene above the value provided.}
#'   \item{expressed between:}{Similar to \code{expressed below:}, but provide
#'   two values separated by spaces. For example
#'   \code{expressed between: ACT5 2 5.5, ACT2 1 2.7}. This descriptor will
#'   only allow training cells expressing the given gene between the two values
#'   provided.}
#' }
#' @section Checking your marker file:
#'
#' Because only specific expressed markers are useful for Garnett
#' classification, we recommend that you always check your marker file for
#' ambiguity before proceeding with classification. To do this, we have
#' provided the functions \code{\link{check_markers}} and
#' \code{\link{plot_markers}} to facilitate marker checking. See that manual
#' pages for those functions for details.
#'
#' @seealso \code{\link{train_cell_classifier}}
#'
Parser <- R6::R6Class(
  "Parser",
  public = list(
    tokens = c("KEYWORD", "NAME", "NEWLINE",
               "SIMP_KEY", "NUM"),
    literals = c(">", ":", ","),

    # Parsing rules
    linenum = 1,
    p_file_1 = function(doc="file : NEWLINE file", p) {
      p$set(1, p$get(3))
    },
    p_file_2 = function(doc="file : celltype
                                | file celltype", p) {
      if(p$length() == 2) {
        cts = new.env(hash=TRUE)
        name_order <- list()
        cts[["name_order"]] <- name_order
        cell_ont <- p$get(2)
        cts[[cell_ont@name]] <- cell_ont
        cts[["name_order"]] <- c(cts[["name_order"]], cell_ont@name)
      } else {
        cts <- p$get(2)
        cell_ont <- p$get(3)
        if(cell_ont@name %in% names(cts)) {
          stop(paste0("Cell type '", cell_ont@name,
                      "' is defined a second time at or near line ",
                      self[["linenum"]], "."))
        }
        cts[[cell_ont@name]] <- cell_ont
        cts[["name_order"]] <- c(cts[["name_order"]], cell_ont@name)
      }
      p$set(1, cts)
    },
    p_header_start = function(doc="header_start : '>' NAME
                                                | '>' NUM
                                                | header_start NUM
                                                | header_start NAME", p) {
      if (is.character(p$get(2))) {
        ct <- new("cell_rules")
        ct@name <- p$get(3)
      } else {
        ct <- p$get(2)
        ct@name <- paste(ct@name, p$get(3))
      }
      p$set(1, ct)
    },

    p_header_1 = function(doc="header : header_start NEWLINE", p) {
      self[["linenum"]] <- self[["linenum"]] + 1
      p$set(1, p$get(2))
    },
    p_celltype_1 = function(doc="celltype : header", p) {
      p$set(1, p$get(2))
    },
    p_celltype_2 = function(doc="celltype : celltype rule", p) {
      ct <- p$get(2)
      rule <- p$get(3)

      if (rule[[1]] == "exp") {
        ct@expressed <- c(ct@expressed, rule[2:length(rule)])
      } else if (rule[[1]] == "nexp") {
        ct@not_expressed <- c(ct@not_expressed, rule[2:length(rule)])
      } else if (rule[[1]] == "ref") {
        ct@references <- c(ct@references, rule[2:length(rule)])
      } else if (rule[[1]]  == "sub") {
        ct@parenttype <- rule[[2]]
      } else if (rule[[1]] == "meta") {
        if (nrow(ct@meta) == 0) {
          ct@meta <- rule[[2]]
        } else {
          ct@meta <- rbind(ct@meta, rule[[2]])
        }
      } else if (rule[[1]] %in% c("exp_ab", "exp_bel", "exp_bet")) {
        ct@gene_rules <- c(ct@gene_rules, rule[2:length(rule)])
      }
      p$set(1, ct)
    },
    p_rule_head_1 = function(doc="rule_head : SIMP_KEY ':'", p) {
      p$set(1, p$get(2))
    },
    p_rule_head_comp_1 = function(doc="rule_head_comp : KEYWORD ':'", p) {
      p$set(1, p$get(2))
    },
    p_rule_head_2 = function(doc="rule_head : NAME ':'", p) {
      p$set(1, p$get(2))
    },
    p_rule_1 = function(doc="rule : rule_head expression NEWLINE", p) {
      self[["linenum"]] <- self[["linenum"]] + 1
      if(p$get(2) == "expressed") {
        p$set(1, c("exp", p$get(3)))
      } else if(p$get(2) == "not expressed") {
        p$set(1, c("nexp", p$get(3)))
      } else if(p$get(2) == "references") {
        p$set(1, c("ref", p$get(3)))
      } else {
        p$set(1, list("meta", data.frame(name = p$get(2),
                                         spec = p$get(3),
                                         stringsAsFactors = FALSE)))
      }
    },
    p_rule_2 = function(doc="rule : rule_head_comp comp_expression NEWLINE",
                        p) {
      self[["linenum"]] <- self[["linenum"]] + 1
      if(!is.list(p$get(3))) p$set(3, list(p$get(3)))
      if(p$get(2) == "expressed above") {
        grs <- lapply(p$get(3), function(x) {
          if (length(x) != 2) stop(paste0("Syntax error in marker file at or",
                                          " near line ", self[["linenum"]],
                                          ": expressed above needs one ",
                                          "value"))
          x <- new("gene_rule", gene_name = x[1],
                   lower = as.numeric(x[2]), upper = Inf)
        })
        p$set(1, c("exp_ab", grs))
      } else if(p$get(2) == "expressed below") {
        grs <- lapply(p$get(3), function(x) {
          if (length(x) != 2) stop(paste0("Syntax error in marker file at or",
                                          " near line ", self[["linenum"]],
                                          ": expressed below needs one ",
                                          "value"))
          x <- new("gene_rule", gene_name = x[1], upper = as.numeric(x[2]),
                   lower = -1)
        })
        p$set(1, c("exp_bel", grs))
      } else if(p$get(2) == "expressed between") {
        grs <- lapply(p$get(3), function(x) {
          if (length(x) != 3) stop(paste0("Syntax error in marker file at or",
                                          " near line ", self[["linenum"]],
                                          ": expressed between needs two",
                                          " values"))
          n1 <- as.numeric(x[2])
          n2 <- as.numeric(x[3])
          if (n1 > n2) stop(paste0("expressed between: ", x[1],
                                   " has expression values out of order ",
                                   "- swapping."))
          x <- new("gene_rule", gene_name = x[1], lower = min(n1, n2),
                   upper = max(n1, n2))
        })
        p$set(1, c("exp_bet", grs))
      } else if(p$get(2) == "subtype of") {
        p$set(1, list("sub", paste(unlist(p$get(3)), collapse = " ")))
      }
    },
    p_rule_3 = function(doc="rule : rule NEWLINE", p) {
      p$set(1, p$get(2))
      self[["linenum"]] <- self[["linenum"]] + 1
    },
    p_expression_1 = function(doc="expression : NAME
                                              | NUM", p) {
      p$set(1, p$get(2))
    },
    p_expression_2 = function(doc="expression : expression NAME
                                              | expression NUM", p) {
      if(length(p$get(2)) == 1) {
        p$set(1, paste(p$get(2), p$get(3)))
      } else {
        p$set(1, c(p$get(2)[1:(length(p$get(2)) - 1)],
                   paste(p$get(2)[length(p$get(2))], p$get(3))))
      }

    },
    p_expression_3 = function(doc="expression : expression ',' NAME
                                              | expression ',' NUM", p) {
      p$set(1, c(p$get(2), p$get(4)))
    },
    p_expression_4 = function(doc="expression : expression ',' NEWLINE NAME
                                              | expression ',' NEWLINE NUM",
                              p) {
      p$set(1, c(p$get(2), p$get(5)))
    },
    p_comp_expression_1 = function(doc="comp_expression : NAME
                                                        | NUM" , p) {
      p$set(1, p$get(2))
    },
    p_comp_expression_2 = function(doc="comp_expression : comp_expression NUM
                                                        | comp_expression NAME",
                                   p) {
      p$set(1, c(p$get(2), p$get(3)))
    },
    p_comp_expression_3 =
      function(doc="comp_expression : comp_expression ',' NAME NUM", p) {
        l <- p$get(2)
        if(!is.list(l)) l <- list(l)
        l[[length(l) + 1]] <- c(p$get(4), p$get(5))
        p$set(1, l)
      },
    p_comp_expression_4 =
      function(doc="comp_expression : comp_expression ',' NEWLINE NAME NUM",
               p) {
        l <- p$get(2)
        if(!is.list(l)) l <- list(l)
        l[[length(l) + 1]] <- c(p$get(5), p$get(6))
        p$set(1, l)
      },
    p_comp_expression_5 =
      function(doc="comp_expression : comp_expression ',' NAME NUM NUM", p) {
        l <- p$get(2)
        if(!is.list(l)) l <- list(l)
        l[[length(l) + 1]] <- c(p$get(4), p$get(5), p$get(6))
        p$set(1, l)
      },
    p_comp_expression_6 =
      function(doc="comp_expression : comp_expression ',' NEWLINE NAME NUM NUM",
               p) {
        l <- p$get(2)
        if(!is.list(l)) l <- list(l)
        l[[length(l) + 1]] <- c(p$get(5), p$get(6), p$get(7))
        p$set(1, l)
      },
    p_error = function(p) {
      if(is.null(p)) stop("Marker file error. Syntax error at EOF")
      else           stop(sprintf("Marker file error. Syntax error '%s' at or near line number %s\n",
                                  p$value, self[["linenum"]]))
    }
  )
)
# --------------------------------------------------------------------

# Part2. Methods for gene selection
# --------------------------------------------------------------------
cds_to_other_id <- function(cds,
                            db,
                            input_file_gene_id_type,
                            new_gene_id_type,
                            verbose = FALSE) {
  matrix <- exprs(cds)
  fdata <- fData(cds)

  new_g <- convert_gene_ids(row.names(fdata),
                            db,
                            input_file_gene_id_type,
                            new_gene_id_type)
  lstart <- length(new_g)
  new_g <- new_g[!is.na(new_g)]
  new_g <- new_g[!duplicated(new_g)]
  lend <- length(new_g)
  if((lstart-lend)/lstart > .7) warning(paste("More than 70% of IDs were lost",
                                              "when converting to",
                                              new_gene_id_type, "IDs. Did you",
                                              "specify the correct gene ID",
                                              "types and the correct db?"))
  if(verbose) message(paste("After converting CDS to", new_gene_id_type,"IDs,",
                            lstart - lend, "IDs were lost"))

  matrix <- matrix[names(new_g),]
  fdata <- fdata[names(new_g),, drop=FALSE]
  row.names(matrix) <- new_g
  row.names(fdata) <- new_g

  pd = new("AnnotatedDataFrame", data = pData(cds))
  fd = new("AnnotatedDataFrame", data = fdata)
  cds = suppressWarnings(newCellDataSet(matrix,
                                        phenoData=pd,
                                        featureData=fd,
                                        expressionFamily=cds@expressionFamily,
                                        lowerDetectionLimit=cds@lowerDetectionLimit))

  return(cds)
}

convert_gene_ids <- function(gene_list,
                             db,
                             start_type,
                             end_type) {

  tryCatch({suppressMessages(AnnotationDbi::mapIds(db, keys = gene_list,
                                                   column = end_type, start_type))},
           error = function(e) {
             msg <- paste0("Garnett cannot convert the gene IDs using the ",
                           "db and types provided. Please check that your db, ",
                           "cds_gene_id_type and marker_file_gene_id_type ",
                           "parameters are correct. Please note that the ", "
                         cds_gene_id_type refers to the type of the ",
                           "row.names of the feature (gene) table in your cds. ",
                           "Conversion error: ", e)
             stop(msg)
           })
}

#' cell_rules class
#'
#' Representation of cell type derived from marker file.
#'
#' @slot name character. Name of the cell type.
#' @slot gene_names character. A list of all of the genes included in the
#'  definition.
#' @slot expressed character. A list of genes defined as "expressed:".
#' @slot not_expressed character. A list of genes defined as "not expressed:".
#' @slot gene_rules vector of GeneRules-class. A list of genes defined under
#'  specific rules using "expressed below:", "expressed above:", or
#'  "expressed between:".
#' @slot meta data.frame of meta data rules specified in marker file.
#' @slot parenttype character. The name of the parent type - specified by
#'  "subtype of:".
#' @slot references character. A list of references included in the definition.
#'
#' @name cell_rules
#' @rdname cell_rules
#' @aliases cell_rules-class
#' @exportClass cell_rules
setClass("cell_rules", representation(name = "character",
                                      gene_names = "character",
                                      expressed = "character",
                                      not_expressed = "character",
                                      gene_rules = "vector",
                                      meta = "data.frame",
                                      parenttype = "character",
                                      references = "character"))

setGeneric("collect_genes", signature = "x", function(x)
  standardGeneric("collect_genes"))

setMethod("collect_genes", c("x" = "cell_rules"), function(x) {
  n1 <- as.character(c(x@expressed, x@not_expressed))
  n2 <- as.character(unlist(lapply(x@gene_rules, function(y) y@gene_name)))
  x@gene_names <- unique(c(n1, n2))
  x
})

collect_gene_names <- function(cellont_list) {
  genes <- lapply(cellont_list, function(x) {
    if(length(collect_genes(x)@gene_names) != 0) {
      pt <- ifelse(identical(x@parenttype, character(0)), "root", x@parenttype)
      data.frame(genes = collect_genes(x)@gene_names, parent = pt,
                 cell_type = x@name)
    }
  })
  all <- do.call("rbind",genes)
  return(all)
}

# Part3. Methods for garnett classifier
# --------------------------------------------------------------------
new_garnett_classifier <- function()
{
  garc <- new( "garnett_classifier",
               classification_tree = igraph::graph.empty())

  root_node_id <- "root"

  garc@classification_tree <- garc@classification_tree +
    igraph::vertex(root_node_id,
                   classify_func=list(function(x) {rep(TRUE, ncol(x))}),
                   model = NULL)

  return(garc)
}

add_cell_type <- function(classifier,
                          cell_type_name,
                          classify_func,
                          parent_cell_type_name="root")
{
  if (cell_type_name %in% igraph::V(classifier@classification_tree)$name){
    stop(paste("Error: cell type",cell_type_name, "already exists."))
  }

  classifier@classification_tree <- classifier@classification_tree +
    igraph::vertex(cell_type_name, classify_func=list(classify_func),
                   model=NULL)

  classifier@classification_tree <- classifier@classification_tree +
    igraph::edge(parent_cell_type_name, cell_type_name)
  return (classifier)
}

make_predictions <- function(cds,
                             classifier,
                             curr_node,
                             rank_prob_ratio,
                             cores = 1,
                             s) {
  cvfit <- igraph::V(classifier@classification_tree)[curr_node]$model[[1]]

  predictions <- tryCatch({
    if(is.null(cvfit)) {
      child_cell_types <- igraph::V(classifier@classification_tree)[
        suppressWarnings(outnei(curr_node)) ]$name
      predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                            ncol=length(child_cell_types),
                            dimnames=list(row.names(pData(cds)),
                                          child_cell_types))
      predictions <- split(predictions, rep(1:ncol(predictions),
                                            each = nrow(predictions)))
      names(predictions) <- child_cell_types
      predictions
    } else {
      candidate_model_genes <- cvfit$glmnet.fit$beta[[1]]@Dimnames[[1]]
      good_genes <- intersect(row.names(exprs(cds)),
                              candidate_model_genes)
      if (length(good_genes) == 0) stop(paste("None of the model genes are in",
                                              "your CDS object. Did you",
                                              "specify the correct",
                                              "cds_gene_id_type and the",
                                              "correct db?"))
      x <- Matrix::t(exprs(cds[intersect(row.names(exprs(cds)),
                                         candidate_model_genes),])) #slow

      extra <- as(matrix(0, nrow = nrow(x),
                         ncol = length(setdiff(candidate_model_genes,
                                               colnames(x)))), "sparseMatrix")
      row.names(extra) <- row.names(x)
      colnames(extra) <- setdiff(candidate_model_genes, colnames(x))

      x <- cbind(x, extra)
      x <- x[,candidate_model_genes]

      # predict probabilities using fitted model
      nonz <- Matrix::rowSums(do.call(cbind,
                                      glmnet:::coef.glmnet(cvfit,
                                                           s="lambda.min")))
      nonz <- nonz[2:length(nonz)]
      nonz <- names(nonz[nonz != 0])

      if (sum(!nonz %in% row.names(exprs(cds))) > 0) {
        warning(paste("The following genes used in the classifier are not",
                      "present in the input CDS. Interpret with caution.",
                      nonz[!nonz %in% row.names(exprs(cds))]))
      }

      temp <- stats::predict(cvfit, #slow
                             newx = x,
                             s = s,
                             type = "response")
      temp[is.nan(temp)] <- 0
      prediction_probs <- as.matrix(as.data.frame(temp))

      # normalize probabilities by dividing by max
      prediction_probs <- prediction_probs/Biobase::rowMax(prediction_probs)

      prediction_probs[is.nan(prediction_probs)] <- 0

      # find the odds ratio of top prob over second best
      prediction_probs <- apply(prediction_probs, 1, function(x) {
        m <- names(which.max(x))
        s <- sort(x, decreasing = T)
        c(cell_type = m, odds_ratio = s[1]/s[2])
      })

      prediction_probs <- as.data.frame(t(prediction_probs))
      prediction_probs$cell_name <- row.names(prediction_probs)
      names(prediction_probs) <- c("cell_type", "odds_ratio", "cell_name")
      prediction_probs$odds_ratio <-
        as.numeric(as.character(prediction_probs$odds_ratio))

      # odds ratio has to be larger than rank_prob_ratio
      assignments <- prediction_probs[prediction_probs$odds_ratio >
                                        rank_prob_ratio,]

      # odds ratio also must be larger than expected by random guess
      # (1/number of cell types)
      random_guess_thresh <- 1.0 / length(cvfit$glmnet.fit$beta)
      assignments <- assignments[assignments$odds_ratio > random_guess_thresh,]

      not_assigned <- row.names(pData(cds))[ !row.names(pData(cds)) %in%
                                               assignments$cell_name]
      if(length(not_assigned) > 0) {
        assignments <- rbind(assignments,
                             data.frame(cell_name = not_assigned,
                                        cell_type = NA, odds_ratio = NA))
      }

      assignments$cell_type <- stringr::str_replace_all(assignments$cell_type,
                                                        "\\.1",
                                                        "")

      # reformat predictions
      predictions <- reshape2::dcast(assignments, cell_name ~ cell_type,
                                     value.var = "odds_ratio")
      predictions <- predictions[!is.na(predictions$cell_name),]
      row.names(predictions) <- predictions$cell_name

      if (ncol(predictions) > 2){
        predictions <- predictions[,setdiff(colnames(predictions), "NA")]
        predictions <- predictions[,-1, drop=FALSE]
        predictions <- predictions[rownames(pData(cds)),,drop=FALSE]
        predictions <- as.matrix(predictions)
        predictions[is.na(predictions)] <- FALSE
        predictions[predictions != 0] <- TRUE
        cell_type_names <- colnames(predictions)

        predictions <- split(predictions, rep(1:ncol(predictions),
                                              each = nrow(predictions)))
        names(predictions) <- cell_type_names

      } else {
        cell_type_names <- names(cvfit$glmnet.fit$beta)
        one_type <- names(predictions)[2]
        if (one_type == "NA") {
          names(predictions)[2] <- "Unknown"
          one_type <- "Unknown"
        }
        predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                              ncol=length(cell_type_names),
                              dimnames=list(row.names(pData(cds)),
                                            cell_type_names))
        predictions[,one_type] <- TRUE

        predictions <- split(predictions, rep(1:ncol(predictions),
                                              each = nrow(predictions)))
        names(predictions) <- cell_type_names
      }
      predictions
    }

  },
  #warning = function(w) print(w),
  error = function(e) {
    if (e$message == paste("None of the model genes are in your CDS object.",
                           "Did you specify the correct cds_gene_id_type and",
                           "the correct db?"))
      stop(e)
    print (e)
    cell_type_names <- names(cvfit$glmnet.fit$beta)
    predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                          ncol=length(cell_type_names),
                          dimnames=list(row.names(pData(cds)),
                                        cell_type_names))
    predictions <- split(predictions, rep(1:ncol(predictions),
                                          each = nrow(predictions)))
    names(predictions) <- cell_type_names
    predictions
  })

  for (i in 1:length(predictions)){
    p <- as(as(predictions[[i]], "sparseVector"), "sparseMatrix")
    row.names(p) <- row.names(pData(cds))
    predictions[[i]] <- p
  }

  return(predictions)
}
# --------------------------------------------------------------------


# Part4. Methods for get training sample
# --------------------------------------------------------------------
get_training_sample <- function(cds,
                                orig_cds,
                                classifier,
                                tf_idf,
                                gene_table,
                                curr_node,
                                parse_list,
                                name_order,
                                max_training_samples,
                                num_unknown,
                                back_cutoff,
                                training_cutoff,
                                marker_scores,
                                return_initial_assign) {

  ##### Find type assignment from expressed/not expressed #####

  child_cell_types <- igraph::V(classifier@classification_tree)[
    suppressWarnings(outnei(curr_node)) ]$name
  parent <- igraph::V(classifier@classification_tree)[curr_node]$name
  if (length(child_cell_types) > 0) {
    if (length(intersect(child_cell_types,
                         colnames(marker_scores))) == 0) {
      return(NULL)
    }
    assigns <- assign_type(marker_scores[,intersect(child_cell_types,
                                                    colnames(marker_scores)),
                                         drop=FALSE],
                           training_cutoff, return_initial_assign)
    if(return_initial_assign) {
      return(assigns)
    }
  }

  assigns <- as.data.frame(assigns)

  names(assigns) <- "assigns"
  pData(orig_cds)$assigns <- assigns[row.names(pData(orig_cds)),"assigns"]
  pData(orig_cds)$assigns <- as.character(pData(orig_cds)$assigns)
  pData(orig_cds)$assigns[is.na(pData(orig_cds)$assigns)] <- "None"

  ##### Exclude possibles using other definitions #####
  child_rules <- list()
  for (child in child_cell_types) {
    cell_class_func <-
      igraph::V(classifier@classification_tree)[child]$classify_func[[1]]

    parent <- environment(cell_class_func)
    if (is.null(parent))
      parent <- emptyenv()
    e1 <- new.env(parent=parent)

    Biobase::multiassign(names(pData(orig_cds)), pData(orig_cds), envir=e1)
    environment(cell_class_func) <- e1

    type_res <- cell_class_func(exprs(orig_cds))
    if (length(type_res)!= ncol(orig_cds)){
      message(paste("Error: classification function for",
                    igraph::V(classifier@classification_tree)[child]$name,
                    "returned a malformed result."))
      stop()
    }

    type_res <- as(as(type_res,"sparseVector"), "sparseMatrix")
    row.names(type_res) <- row.names(pData(orig_cds))
    colnames(type_res) <- child
    child_rules[[ child ]] <- type_res
  }


  # Assign cell types by classifier rules and check for conflicts
  ctf_cell_type <- rep("Unknown", length(child_rules[[1]]))
  names(ctf_cell_type) <- row.names(child_rules[[1]])
  for (child in child_cell_types) {
    ctf_cell_type[Matrix::which(as.matrix(child_rules[child][[1]]))] <- child
  }

  # Find training sample cells and downsample if necessary
  good_cell_type <- ctf_cell_type[ctf_cell_type %in% child_cell_types]
  obs_counts <- table(good_cell_type)
  training_sample <- c()
  for(i in names(obs_counts)){
    num_obs_for_type_i <- min(max_training_samples, obs_counts[i])
    obs_for_type_i <- sample(which(good_cell_type == i), num_obs_for_type_i)
    training_sample <- append(training_sample, obs_for_type_i)
  }
  training_sample <- good_cell_type[training_sample]

  # Find outgroup samples and downsample if necessary
  outgroup_samples <- !ctf_cell_type %in% child_cell_types

  if (length(outgroup_samples) > 0){
    outgroup_samples <- ctf_cell_type[outgroup_samples]

    out_group_cds <- cds[,names(outgroup_samples)]
    out_group_cds <- out_group_cds[,sample(row.names(pData(out_group_cds)),
                                           min(nrow(pData(out_group_cds)),
                                               num_unknown * 10),
                                           replace = F)]

    out_group_cds <- get_communities(out_group_cds)

    per_clust <-
      floor(num_unknown/length(unique(pData(out_group_cds)$louv_cluster)))

    outg <- lapply(unique(pData(out_group_cds)$louv_cluster), function(x) {
      sub <- pData(out_group_cds)[pData(out_group_cds)$louv_cluster == x,]
      sample(row.names(sub), min(nrow(sub), per_clust))
    })

    outg <- unlist(outg)
    if(length(outg) < min(length(outgroup_samples), num_unknown)) {
      to_get <- min(length(outgroup_samples), num_unknown)
      outgroup_samples <- outgroup_samples[!names(outgroup_samples) %in% outg]
      outg <- c(outg, sample(names(outgroup_samples), (to_get - length(outg))))
    }
    outgroup <- rep("Unknown", length(outg))
    names(outgroup) <- outg
    training_sample <- append(training_sample, outgroup)
  }

  training_sample <- factor(training_sample)
  training_sample <- droplevels(training_sample)

  return(training_sample)
}

assign_type <- function(total_vals,
                        training_cutoff,
                        return_initial_assign) {
  cutoffs <- apply(total_vals, 2, stats::quantile, probs = training_cutoff)

  q <- as.data.frame(t(t(total_vals) > cutoffs))
  q$total <- rowSums(q)
  q$assign <- sapply(q$total, function(x) {
    if (x == 0) {
      out <- "None"
    } else if (x > 1) {
      out <- "Ambiguous"
    } else {
      out <- "Assign"
    }
    out
  })
  q$cell <- row.names(total_vals)
  if (return_initial_assign) {
    q$total <- NULL
    types <- colnames(total_vals)
    q$Assignment <- "None"
    for(type in types) {
      q$Assignment[q[,type]] <- type
    }
    q$Assignment[q$assign == "Ambiguous"] <- "Ambiguous"
    q$cell <- NULL
    q$assign <- NULL
    return(q)
  }

  out <- q[q$assign != "Assign",]$assign
  names(out) <- q[q$assign != "Assign",]$cell
  sub <- q[q$assign == "Assign",]
  sub$total <- NULL
  sub$assign <- NULL
  sub <- reshape2::melt(sub, id.vars = "cell")

  sub <- sub[sub$value,]
  out2 <- as.character(sub$variable)
  names(out2) <- sub$cell
  out <- c(out, out2)
  out <- out[out != "None"]
  out
}

aggregate_positive_markers <- function(cell_type,
                                       tf_idf,
                                       gene_table,
                                       back_cutoff,
                                       agg = TRUE) {
  gene_list <- cell_type@expressed
  bad_genes <- gene_table[!gene_table$in_cds,]$orig_fgenes
  gene_list <- gene_list[!gene_list %in% bad_genes]
  gene_list <- gene_table$fgenes[match(gene_list,gene_table$orig_fgenes)]

  rel <- intersect(unlist(gene_list), colnames(tf_idf))
  rel <- unique(rel)
  rel_genes <- tf_idf[,rel, drop=FALSE]
  if(length(rel) == 0) return(NULL)

  if(length(unlist(gene_list)) == 1) {
    background_cutoff <- back_cutoff * stats::quantile(rel_genes[,1],
                                                       prob = .95)
    rel_genes[rel_genes < background_cutoff] <- 0
    agg <- rel_genes
  } else {
    background_cutoff <- back_cutoff * apply(rel_genes, 2,
                                             stats::quantile,
                                             prob = .95)

    temp <- Matrix::t(rel_genes)
    temp[temp < background_cutoff] <- 0
    rel_genes <- Matrix::t(temp)

    if(agg) {
      agg <-  (Matrix::rowSums(rel_genes))
    } else {
      agg <- rel_genes
    }
  }
  agg
}

get_negative_markers <- function(cell_type,
                                 tf_idf,
                                 gene_table,
                                 back_cutoff) {
  not_list <- cell_type@not_expressed
  if (length(not_list) != 0) {
    bad_genes <- gene_table[!gene_table$in_cds,]$orig_fgenes
    not_list <- not_list[!not_list %in% bad_genes]
    not_list <- gene_table$fgenes[match(not_list,gene_table$orig_fgenes)]
    rel <- intersect(unlist(not_list), colnames(tf_idf))
    rel <- unique(rel)
    rel_genes <- tf_idf[,rel, drop=FALSE]
    if(length(rel) == 0) return("")

    if(length(unlist(not_list)) == 1) {
      background_cutoff <- back_cutoff * stats::quantile(rel_genes[,1],
                                                         prob = .95)
      bad_cells <- names(rel_genes[rel_genes[,1] > background_cutoff,])
    } else {
      background_cutoff <- back_cutoff * apply(rel_genes, 2,
                                               stats::quantile,
                                               prob = .95)

      temp <- Matrix::t(rel_genes)
      temp <- temp > background_cutoff
      temp <- Matrix::t(temp)
      bad_cells <- row.names(temp[Matrix::rowSums(temp) != 0,])
    }
  } else {
    bad_cells <- ""
  }
  return(bad_cells)
}

#' This function takes single-cell expression data in the form of a CDS object
#' and a cell type definition file (marker file) and trains a multinomial
#' classifier to assign cell types. The resulting \code{garnett_classifier}
#' object can be used to classify the cells in the same dataset, or future
#' datasets from similar tissues/samples.
#'
#' @param cds Input CDS object.
#' @param marker_file A character path to the marker file to define cell types.
#'  See details and documentation for \code{\link{Parser}} by running
#'  \code{?Parser}for more information.
#' @param db Bioconductor AnnotationDb-class package for converting gene IDs.
#'  For example, for humans use org.Hs.eg.db. See available packages at
#'  \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor}.
#'  If your organism does not have an AnnotationDb-class database available,
#'  you can specify "none", however then Garnett will not check/convert gene
#'  IDs, so your CDS and marker file must have the same gene ID type.
#' @param cds_gene_id_type The type of gene ID used in the CDS. Should be one
#'  of the values in \code{columns(db)}. Default is "ENSEMBL". Ignored if
#'  db = "none".
#' @param marker_file_gene_id_type The type of gene ID used in the marker file.
#'  Should be one of the values in \code{columns(db)}. Default is "SYMBOL".
#'  Ignored if db = "none".
#' @param min_observations An integer. The minimum number of representative
#'  cells per cell type required to include the cell type in the predictive
#'  model. Default is 8.
#' @param max_training_samples An integer. The maximum number of representative
#'  cells per cell type to be included in the model training. Decreasing this
#'  number increases speed, but may hurt performance of the model. Default is
#'  500.
#' @param num_unknown An integer. The number of unknown type cells to use as an
#'  outgroup during classification. Default is 500.
#' @param propogate_markers Logical. Should markers from child nodes of a cell
#'  type be used in finding representatives of the parent type? Should
#'  generally be \code{TRUE}.
#' @param cores An integer. The number of cores to use for computation.
#' @param lambdas \code{NULL} or a numeric vector. Allows the user to pass
#'  their own lambda values to \code{\link[glmnet]{cv.glmnet}}. If \code{NULL},
#'  preset lambda values are used.
#' @param classifier_gene_id_type The type of gene ID that will be used in the
#'  classifier. If possible for your organism, this should be "ENSEMBL", which
#'  is the default. Ignored if db = "none".
#' @param return_initial_assign Logical indicating whether an initial
#'  assignment data frame for the root level should be returned instead of a
#'  classifier. This can be useful while choosing/debugging markers. Please
#'  note that this means that a classifier will not be built, so you will not
#'  be able to move on to the next steps of the workflow until you rerun the
#'  functionwith \code{return_initial_assign = FALSE}. Default is \code{FALSE}.
#'
#' @details This function has three major parts: 1) parsing the marker file 2)
#'  choosing cell representatives and 3) training the classifier. Details on
#'  each of these steps is below:
#'
#'  Parsing the marker file: the first step of this function is to parse the
#'  provided marker file. The marker file is a representation of the cell types
#'  expected in the data and known characteristics about them. Information
#'  about marker file syntax is available in the documentation for the
#'  \code{\link{Parser}} function, and on the
#'  \href{https://cole-trapnell-lab.github.io/garnett}{Garnett website}.
#'
#'  Choosing cell representatives: after parsing the marker file, this function
#'  identifies cells that fit the parameters specified in the file for each cell
#'  type. Depending on how marker genes and other cell type definition
#'  information are specified, expression data is normalized and expression
#'  cutoffs are defined automatically. In addition to the cell types in the
#'  marker file, an outgroup of diverse cells is also chosen.
#'
#'  Training the classifier: lastly, this function trains a multinomial GLMnet
#'  classifier on the chosen representative cells.
#'
#'  Because cell types can be defined hierarchically (i.e. cell types can be
#'  subtypes of other cell types), steps 2 and 3 above are performed iteratively
#'  over all internal nodes in the tree representation of cell types.
#'
#'  See the
#'  \href{https://cole-trapnell-lab.github.io/garnett}{Garnett website} and the
#'  accompanying paper for further details.
#'
#' @export
#'
#' @import garnett
#' @import Biobase
select_fine_samples <- function(cds,
                                marker_file,
                                db,
                                cds_gene_id_type = "ENSEMBL",
                                marker_file_gene_id_type = "SYMBOL",
                                cutoff = .75,
                                min_observations = 8,
                                max_training_samples = 500,
                                num_unknown = 500,
                                propogate_markers = TRUE,
                                cores = 1,
                                lambdas = NULL,
                                classifier_gene_id_type = "ENSEMBL",
                                return_initial_assign = TRUE) {
  # sink("E:/scCancer2/vignettes/garnett-output.txt", append = TRUE, split = FALSE)
  ##### Check inputs #####
  assertthat::assert_that(is(cds, "CellDataSet"))
  assertthat::assert_that(assertthat::has_name(pData(cds), "Size_Factor"),
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling train_cell_classifier"))
  assertthat::assert_that(sum(is.na(pData(cds)$Size_Factor)) == 0,
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling train_cell_classifier"))
  assertthat::assert_that(is.character(marker_file))
  assertthat::is.readable(marker_file)
  if (is(db, "character") && db == "none") {
    cds_gene_id_type <- 'custom'
    classifier_gene_id_type <- 'custom'
    marker_file_gene_id_type <- 'custom'
  } else {
    assertthat::assert_that(is(db, "OrgDb"),
                            msg = paste0("db must be an 'AnnotationDb' object ",
                                         "or 'none' see ",
                                         "http://bioconductor.org/packages/",
                                         "3.8/data/annotation/ for available"))
    assertthat::assert_that(is.character(cds_gene_id_type))
    assertthat::assert_that(is.character(marker_file_gene_id_type))
    assertthat::assert_that(cds_gene_id_type %in% AnnotationDbi::keytypes(db),
                            msg = paste("cds_gene_id_type must be one of",
                                        "keytypes(db)"))
    assertthat::assert_that(classifier_gene_id_type %in% AnnotationDbi::keytypes(db),
                            msg = paste("classifier_gene_id_type must be one of",
                                        "keytypes(db)"))
    assertthat::assert_that(marker_file_gene_id_type %in%
                              AnnotationDbi::keytypes(db),
                            msg = paste("marker_file_gene_id_type must be one of",
                                        "keytypes(db)"))
  }
  assertthat::is.count(num_unknown)
  assertthat::is.count(cores)
  assertthat::assert_that(is.logical(propogate_markers))
  if (!is.null(lambdas)) {
    assertthat::assert_that(is.numeric(lambdas))
  }

  ##### Set internal parameters #####
  rel_gene_quantile <- .9 # exclusion criterion for genes expressed at greater
  # than rel_gene_quantile in all training cell subsets
  back_cutoff <- 0.25 # percent of 95th percentile of expression that marks the
  # cutoff between "expressed" and "not expressed"
  perc_cells <- 0.05 # percent of training cells a gene is expressed to be
  # included in glmnet training
  training_cutoff <- cutoff # percentile of marker score required for training
  # assignment

  ##### Normalize and rename CDS #####
  if (!is(exprs(cds), "dgCMatrix")) {
    sf <- pData(cds)$Size_Factor
    pd <- new("AnnotatedDataFrame", data = pData(cds))
    fd <- new("AnnotatedDataFrame", data = fData(cds))
    cds <- suppressWarnings(newCellDataSet(as(exprs(cds), "dgCMatrix"),
                                           phenoData = pd,
                                           featureData = fd))
    pData(cds)$Size_Factor <- sf
  }

  # pData(cds)$num_genes_expressed <- Matrix::colSums(as(exprs(cds), "lgCMatrix"))
  pData(cds)$num_genes_expressed <- Matrix::colSums(as(exprs(cds), "lMatrix"))
  cell_totals <-  Matrix::colSums(exprs(cds))
  sf <- pData(cds)$Size_Factor

  pd <- new("AnnotatedDataFrame", data = pData(cds))
  fd <- new("AnnotatedDataFrame", data = fData(cds))
  temp <- exprs(cds)
  temp@x <- temp@x / rep.int(pData(cds)$Size_Factor, diff(temp@p))
  norm_cds <- suppressWarnings(newCellDataSet(temp,
                                              phenoData = pd, featureData = fd))
  orig_cds <- cds
  if(cds_gene_id_type != classifier_gene_id_type)  {
    norm_cds <- cds_to_other_id(norm_cds, db=db, cds_gene_id_type,
                                classifier_gene_id_type)
    orig_cds <- cds_to_other_id(cds, db=db, cds_gene_id_type,
                                classifier_gene_id_type)
  }
  pData(norm_cds)$Size_Factor <- sf


  ##### Parse Marker File #####
  file_str = paste0(readChar(marker_file, file.info(marker_file)$size),"\n")

  parse_list <- parse_input(file_str)
  orig_name_order <- unlist(parse_list[["name_order"]])
  rm("name_order", envir=parse_list)

  # Check and order subtypes
  ranks <- lapply(orig_name_order, function(i) parse_list[[i]]@parenttype)
  names(ranks) <- orig_name_order

  if(length(unlist(unique(ranks[which(!ranks %in% names(ranks) & lengths(ranks) != 0L)])) != 0)) {
    stop(paste("Subtype", unlist(unique(ranks[which(!ranks %in% names(ranks) & lengths(ranks) != 0L)])), "is not defined in marker file."))
  }

  if(any(names(ranks) == ranks)) {
    bad <- ranks[names(ranks) == ranks]
    stop(paste0("'", bad,
                "' cannot be a subtype of itself. Please modify marker file."))
  }

  name_order <- names(ranks[lengths(ranks) == 0L])
  ranks <- ranks[!names(ranks) %in% name_order]
  while(length(ranks) != 0) {
    name_order <- c(name_order, names(ranks)[ranks %in% name_order])
    ranks <- ranks[!names(ranks) %in% name_order]
  }


  if(is.null(parse_list)) stop("Parse failed!")
  message(paste("There are", length(parse_list), "cell type definitions"))

  # Check gene names and keywords
  gene_table <- make_name_map(parse_list,
                              as.character(row.names(fData(norm_cds))),
                              classifier_gene_id_type,
                              marker_file_gene_id_type,
                              db)

  ##### Make garnett_classifier #####
  classifier <- new_garnett_classifier()
  classifier@gene_id_type <- classifier_gene_id_type
  if(is(db, "character") && db == "none") classifier@gene_id_type <- "custom"

  for(i in name_order) {
    # check meta data exists
    if (nrow(parse_list[[i]]@meta) != 0) {
      if (!all(parse_list[[i]]@meta$name %in% colnames(pData(norm_cds)))) {
        bad_meta <- parse_list[[i]]@meta$name[!parse_list[[i]]@meta$name %in%
                                                colnames(pData(norm_cds))]
        stop(paste0("Cell type '", parse_list[[i]]@name,
                    "' has a meta data specification '", bad_meta ,
                    "' that's not in the pData table."))
      }
    }
    logic_list <- assemble_logic(parse_list[[i]], gene_table)
    classifier <- add_cell_rule(parse_list[[i]], classifier, logic_list)
  }

  classifier@cell_totals <- exp(mean(log(cell_totals)))/
    stats::median(pData(norm_cds)$num_genes_expressed)

  ##### Create transformed marker table #####
  if(propogate_markers) {
    root <- propogate_func(curr_node = "root", parse_list, classifier)
  }

  tf_idf <- tfidf(norm_cds) #slow


  ### Aggregate markers ###
  marker_scores <- data.frame(cell = row.names(tf_idf))

  for (i in name_order) {
    agg <- aggregate_positive_markers(parse_list[[i]], tf_idf,
                                      gene_table, back_cutoff)
    bad_cells <- get_negative_markers(parse_list[[i]], tf_idf,
                                      gene_table, back_cutoff)
    if(is.null(agg))  {
      warning (paste("Cell type", i, "has no genes that are expressed",
                     "and will be skipped"))
    } else {
      agg[names(agg) %in% bad_cells] <- 0
      marker_scores <- cbind(marker_scores, as.matrix(agg))
      colnames(marker_scores)[ncol(marker_scores)] <- parse_list[[i]]@name
    }
  }
  # return(marker_scores)
  ##### Sample Selection #####
  for (v in igraph::V(classifier@classification_tree)){
    child_cell_types <- igraph::V(classifier@classification_tree)[
      suppressWarnings(outnei(v))]$name

    if(length(child_cell_types) > 0) {
      ### Get CDS subset for training ###
      message(igraph::V(classifier@classification_tree) [ v ]$name)
      if(igraph::V(classifier@classification_tree) [ v ]$name == "root") {
        cds_sub <- norm_cds
        orig_sub <- orig_cds
      } else {
        # loosely classify to subset
        new_assign <-
          make_predictions(norm_cds,
                           classifier,
                           igraph::V(classifier@classification_tree)[
                             suppressWarnings(innei(v))]$name,
                           rank_prob_ratio = 1.1,
                           s = "lambda.min")
        if(!igraph::V(classifier@classification_tree)[v]$name %in%
           names(new_assign)) {
          message(paste0("No cells classified as ",
                         igraph::V(classifier@classification_tree) [ v ]$name,
                         ". No subclassification"))
          next
        }
        good_cells <-
          as.matrix(new_assign[
            igraph::V(classifier@classification_tree)[v]$name][[1]])
        good_cells <- names(good_cells[good_cells[,1] != 0,])
        if(length(good_cells) == 0) {
          message(paste0("No cells classified as ",
                         igraph::V(classifier@classification_tree) [ v ]$name,
                         ". No subclassification"))
          next
        }
        cds_sub <- norm_cds[,good_cells]
        orig_sub <- orig_cds[,good_cells]
      }

      ### Get training sample ###
      sample <- get_training_sample(cds = cds_sub,
                                    orig_cds = orig_sub,
                                    classifier,
                                    tf_idf,
                                    gene_table,
                                    v,
                                    parse_list,
                                    name_order,
                                    max_training_samples,
                                    num_unknown,
                                    back_cutoff,
                                    training_cutoff,
                                    marker_scores,
                                    return_initial_assign)

      if(return_initial_assign) {
        # sink()
        return(sample)
      }
    }
  }
  return(training_samples)
}

parse_input <- function(file_str,
                        debug = F) {
  # Parse input_file
  lexer  <- rly::lex(Lexer, debug=debug)
  parser <- rly::yacc(Parser, debug=debug)
  parse_list <- parser$parse(file_str, lexer)

  parse_list
}

make_name_map <- function(parse_list,
                          possible_genes,
                          cds_gene_id_type,
                          marker_file_gene_id_type,
                          db) {
  gene_start <- collect_gene_names(parse_list)
  gene_table <- data.frame(fgenes = gene_start[,1], parent = gene_start[,2])
  gene_table$parent <- as.character(gene_table$parent)
  gene_table$fgenes <- as.character(gene_table$fgenes)
  gene_table$orig_fgenes <- gene_table$fgenes
  if(cds_gene_id_type != marker_file_gene_id_type) {
    gene_table$fgenes <- convert_gene_ids(gene_table$orig_fgenes,
                                          db,
                                          marker_file_gene_id_type,
                                          cds_gene_id_type)
    bad_convert <- sum(is.na(gene_table$fgenes))
    if (bad_convert > 0) warning(paste(bad_convert,
                                       "genes could not be converted from",
                                       marker_file_gene_id_type,
                                       "to", cds_gene_id_type, "These genes are",
                                       "listed below:", paste0(gene_table$orig_genes[
                                         is.na(gene_table$fgenes)],
                                         collapse="\n")))
  } else {
    gene_table$cds <- gene_table$fgenes
  }

  if(cds_gene_id_type == "ENSEMBL" | marker_file_gene_id_type == "ENSEMBL") {
    gene_table$cds <- NULL
    possibles <- data.frame(cds = possible_genes,
                            ensembl = as.character(
                              stringr::str_split_fixed(possible_genes,
                                                       "\\.",
                                                       2)[,1]))
    gene_table <- merge(gene_table, possibles, all.x=T,
                        by.x="fgenes", by.y="ensembl")
    gene_table$fgenes <- gene_table$cds
  } else {
    gene_table$cds <- gene_table$fgenes
  }

  gene_table$in_cds <- gene_table$f %in% possible_genes
  gene_table$in_cds[is.na(gene_table$in_cds)] <- FALSE

  bad_genes <- gene_table$orig_fgenes[!gene_table$in_cds]
  if (length(bad_genes) > 0) warning(strwrap("The following genes from
                                             the cell type definition file are
                                             not present in the cell dataset.
                                             Please check these genes for
                                             errors. Cell type determination
                                             will continue, ignoring these
                                             genes."), "\n",
                                     paste0(bad_genes, collapse="\n"))

  gene_table$fgenes <- as.character(gene_table$fgenes)
  gene_table$cds <- as.character(gene_table$cds)

  gene_table
}

add_cell_rule <- function(cell_type,
                          classifier,
                          logic_list) {
  # Set parenttype to root if no parent
  if (length(cell_type@parenttype) == 0) {
    cell_type@parenttype <- "root"
  }
  # subtype of
  if (length(cell_type@parenttype) > 1) stop("only 1 parenttype allowed")
  parent_type <- as.character(cell_type@parenttype)

  # references
  if (length(cell_type@references) > 0) {
    if (length(classifier@references) == 0) {
      classifier@references <- list()
    }

    classifier@references <- c(classifier@references,
                               list(cell_type@references))
    names(classifier@references)[length(classifier@references)] <-
      cell_type@name
  }

  if (length(logic_list) == 0) {
    warning (paste("Cell type", cell_type@name,
                   "has no valid rules and will be skipped"))
    classifier <- add_cell_type(classifier, cell_type@name,
                                classify_func = function(x) {rep(FALSE, ncol(x))},
                                parent_type)
    return(classifier)
  }
  logic <- paste(unlist(logic_list), collapse = ' & ')

  tryCatch(
    if(nchar(logic) == 0) {
      classifier <- add_cell_type(classifier, cell_type@name,
                                  classify_func = function(x) {FALSE},
                                  parent_type)

    } else {
      classifier <- add_cell_type(classifier, cell_type@name,
                                  classify_func = function(x) {
                                    eval(parse(text = logic))
                                  },
                                  parent_type)
    },
    error = function(e) {
      msg <- paste("Cell type rule generation failed on the",
                   "cell definition for ", cell_type@name, ".\nError: ",e)
      stop(msg)
    }
  )
  return(classifier)
}

assemble_logic <- function(cell_type,
                           gene_table) {

  logic = ""
  logic_list = list()
  bad_genes <- gene_table[!gene_table$in_cds,]$orig_fgenes

  # expressed/not expressed
  logic_list <- lapply(cell_type@gene_rules, function(rule) {
    log_piece <- ""
    if (!rule@gene_name %in% bad_genes) {
      paste0("(x['",
             gene_table$fgenes[match(rule@gene_name, gene_table$orig_fgenes)],
             "',] > ", rule@lower,
             ") & (x['",
             gene_table$fgenes[match(rule@gene_name, gene_table$orig_fgenes)],
             "',] < ",
             rule@upper,
             ")")
    }
  })
  if (length(cell_type@expressed) > 0 | length(cell_type@not_expressed) > 0) {
    logic_list <- list(logic_list, paste0("assigns == '", cell_type@name, "'"))
  }

  if(length(logic_list) == 0) warning(paste("Cell type", cell_type@name,
                                            "has no valid expression rules."))

  # meta data
  if (nrow(cell_type@meta) > 0) {
    mlogic <- plyr::dlply(cell_type@meta, plyr::.(name), function(x) {
      if(nrow(x) == 1){
        out <- paste0(x["name"], " %in% c('", x[,"spec"][1],"')")
      } else {
        out <- paste0(x[,"name"][1], " %in% c('", paste(x[,"spec"],
                                                        collapse = "', '"),
                      "')")
      }
      out
    })
    logic_list <- c(logic_list, unname(mlogic))
  }
  logic_list <- logic_list[!is.na(logic_list)]
  logic_list
}

propogate_func <- function(curr_node,
                           parse_list,
                           classifier) {
  children <- igraph::V(classifier@classification_tree)[
    suppressWarnings(outnei(curr_node))]$name

  if(length(children) == 0) {
    return(parse_list[[curr_node]]@expressed)
  } else {
    child_genes <- c()
    if (curr_node != "root") {
      child_genes <- parse_list[[curr_node]]@expressed
    }
    for(child in children) {
      child_genes <- union(child_genes,
                           propogate_func(child, parse_list, classifier))
    }
    if(curr_node != "root") {
      parse_list[[curr_node]]@expressed <- child_genes
    }
    return(child_genes)
  }
}

tfidf <- function(input_cds) {
  ncounts <- exprs(input_cds)
  ncounts <- ncounts[Matrix::rowSums(ncounts) != 0,]
  nfreqs <- ncounts
  nfreqs@x <- ncounts@x / rep.int(Matrix::colSums(ncounts), diff(ncounts@p))
  tf_idf_counts <- nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts > 0))
  Matrix::t(tf_idf_counts)
}
