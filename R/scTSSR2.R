#' use scTSSR2 to impute dropout values in scRNA-seq data
#'
#' @param X.count Raw read count matrix. The rows correspond to genes and the columns correspond to cells.
#' Can be sparse.
#'
#' @param k.gene A hyper-parameter that controls the sparsity level of the estimated coefficient matrices, A1 and A2.
#' Default is k_{gene} = min(100, m/30).
#'
#' @param k.cell A hyper-parameter that controls the sparsity level of the estimated coefficient matrices, B1 and B2.
#' Default is k_{cell} = min(100, n/30).
#'
#' @param W A weight matrix with element W_{gc} denotes the non-dropout probability of the expression level of gene g in cell c.
#' Default is W_{gc}=X_{gc}/max(X_{gc}).
#'
#'@param lambda Ridge penalty parameter. Default is 256.
#'
#' @param percent The expression count matrix is preprocessed by filtering out the genes
#' expressed in at most percent*\eqn{100\%} of the cells. Default is 0.05.
#'
#' @param ncores Number of cores to use. Default is 1.
#'
#' @param MAX.ITER Maximum iteration of the external circulation of scTSSR2. Default is 4.
#'
#' @param ABSTOL Absolute tolerance of the external circulation. Default is 1e-3.
#'
#' @param learning.rate A hyper-parameter that controls the speed of adjusting the weights of the network with respect to the loss gradient.
#' Default is 0.0001.
#'
#' @param epochs The number of the entire training set going through the entire network. Default is 100.
#'
#' @param verbose Whether to output the value of metrics at the end of each epoch. Default is TRUE.
#'
#' @param estimates.only logical item. A logical flag to determine whether to output only imputed estimates.
#' Default is FALSE.
#'
#'
#' @return If `estimates.only = TRUE', then a matrix of scTSSR2 estimates.
#'
#' If `estimates.only = FALSE', a list with the following components
#'
#' \item{\code{estimate}}{Recovered (normalized) expression.}
#'
#' \item{\code{se}}{Standard error of estimates.}
#'
#' \item{\code{info}}{Information about dataset.}
#'
#'
#' The \code{info} element is a list with the following components:
#'
#' \item{\code{size.factor}}{Size factor used for normalization.}
#'
#' \item{\code{pred.time}}{Time taken to generate predictions.}
#'
#' \item{\code{posterior.time}}{Time taken to compute the posterior distribution.}
#'
#' \item{\code{total.time}}{Total time for scTSSR2 estimation.}
#'
#'
#' @export
#'
#'
#' @import keras
#'
#' @import tensorflow
#'
#' @import SAVER
#'
#' @author Ke Jin, \email{kej13@mails.ccnu.edu.cn}
#'
#' @examples
#'
#' data("baron")
#'
#' baron_imputation_result = scTSSR2(baron$count.samp)
#'
scTSSR2 <- function(X.count, k.gene = NULL, k.cell = NULL, W = NULL, lambda = 256,

                   percent = 0.05, ncores  = 1, MAX.ITER = 4, ABSTOL = 1e-3,

                   learning.rate = 0.0001, epochs = 100, verbose = TRUE, estimates.only = FALSE){



  # preprocessing and log-normalization

  message("Starting preprocessing and log-normalization ...")

  X <- log.normalization(X.count, percent, preprocess.only = FALSE)

  X.count <- log.normalization(X.count, percent, preprocess.only = TRUE)

  message("Done!")


  # compute weights

  if (is.null(W)) {

    W <- X/max(X)

  }



  if (is.null(k.gene)) {

    k.gene <- round(min(100, nrow(X.count)/30))

  }

  if (is.null(k.cell)) {

    k.cell <- round(min(100, ncol(X.count)/30))

  }



  X.count <- clean.data(X.count)

  ngenes <- nrow(X.count)

  ncells <- ncol(X.count)

  gene.names <- rownames(X.count)

  cell.names <- colnames(X.count)


  # assign size factor

  sf.out <- calc.size.factor(X.count)

  sf <- sf.out[[1]]

  scale.sf <- sf.out[[2]]



  result <- list()

  result$estimate <- matrix(NA, ngenes, ncells, dimnames = list(gene.names, cell.names))

  result$A1 <- matrix(NA, ngenes, k.gene)

  result$A2 <- matrix(NA, k.gene, ngenes)

  result$B1 <- matrix(NA, ncells, k.cell)

  result$B2 <- matrix(NA, k.cell, ncells)





  if (!estimates.only) {

    result$se <- matrix(NA, ngenes, ncells, dimnames = list(gene.names, cell.names))

  } else {

    result$se <- NA

  }

  result$info <- c(list(0), list(0), list(0), list(0))

  names(result$info) <- c("size.factor", "pred.time", "posterior.time", "total.time")

  result$info$size.factor <- scale.sf*sf




  message("Imputation starts ...")

  message("Calculating the prior mean for each gene in each cell ...")

  message("iter", ' objective')

  pred.st <- Sys.time()

  k <- 1

  history.objval <- c()

  zeros = matrix(0, ngenes, ncells)

  while (k <= MAX.ITER){


    batch.size = nrow(X)


    if (k==1){

      Z <- X

    } else {

      Z <- A1%*%A2%*%X

    }


    # using keras with SGD algorithm to calculate B or A

    weight_cells <- keras_weighted_ridge(Y=X, X=Z, k = k.cell, W=W, lambda = lambda, epochs = epochs, batch_size = batch.size,

                                         learning_rate = learning.rate, verbose = verbose)


    B1 <- weight_cells[[1]]

    B2 <- weight_cells[[2]]



    Z <- t(X%*%B1%*%B2)


    batch.size = ncol(X)


    weight_genes <- keras_weighted_ridge(Y=t(X), X=Z, k = k.gene, W=t(W), lambda = lambda, epochs = epochs, batch_size = batch.size,

                                         learning_rate = learning.rate, verbose = verbose)


    A2 <- t(weight_genes[[1]])
    A1 <- t(weight_genes[[2]])



    history.objval[k] <- objective(W=W, X=X, lambda = lambda, A1=A1, A2=A2, B1=B1, B2=B2)


    message(k, '    ', round(history.objval[k], 3))

    if (k > 1 && abs(history.objval[k] - history.objval[k-1]) < ABSTOL) break

    k <- k + 1

  }

  Xhat <- pmax(A1%*%A2%*%X%*%B1%*%B2, zeros)




  pred.time <- Sys.time() - pred.st

  message("Calculating the posterior means with ", ncores, " worker(s)")

  posterior.st <- Sys.time()

  out <- SAVER::saver(X.count, ncores = ncores, mu =  exp(as.matrix(Xhat)) )

  posterior.time <- Sys.time() - posterior.st

  total.time <- Sys.time() - pred.st

  result$estimate <- out$estimate

  result$se <- out$se

  result$info$pred.time <- pred.time

  result$info$posterior.time <- posterior.time

  result$info$total.time <- total.time

  result$A1 <- A1

  result$A2 <- A2

  result$B1 <- B1

  result$B2 <- B2



  message("Done!")

  message("Finish time: ", Sys.time())

  message("Total time: ", format(result$info$total.time))

  if (!estimates.only) {

    class(result) <- "scTSSR2"

    result

  } else {

    result$estimate

  }

}
