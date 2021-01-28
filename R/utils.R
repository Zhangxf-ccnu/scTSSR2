clean.data <- function(x) {

  if (!(grepl("matrix", class(x), ignore.case = TRUE))) {

    x <- Matrix::Matrix(as.matrix(x))

    message("Converting x to matrix.")

    if (!is.numeric(x)) {

      warning("Make sure x is numeric.")

    }

  }

  np <- dim(x)

  size <- as.numeric(np[1])*as.numeric(np[2])

  if(size > 2^31-1){

    inds <- split(1:np[2], ceiling(1:np[2]/1000))

    for(i in 1:length(inds)){

      x[, inds[[i]]][x[, inds[[i]]] < 0.001] <- 0

    }

  } else {

    x[x < 0.001] <- 0

  }

  if (is.null(np) | (np[2] <= 1))

    stop("x should be a matrix with 2 or more columns")

  if (min(Matrix::colSums(x)) == 0) {

    nzerocells <- sum(Matrix::colSums(x) == 0)

    x <- x[, Matrix::colSums(x) != 0]

    message("Removing ", nzerocells, " cell(s) with zero expression.")

  }

  if (is.null(rownames(x))) {

    rownames(x) <- 1:np[1]

  }

  x

}





calc.size.factor <- function(x) {

  sf <- Matrix::colSums(x)/mean(Matrix::colSums(x))

  scale.sf <- 1

  list(unname(sf), scale.sf)

}







regularizer_define <- function(weight_matrix, lambda = 1.0){
  lambda * k_sum(k_square(weight_matrix))
}






loss_define <- function(W, y_true, y_pred){

  W = keras_array(W,'float32')

  0.5 * k_sum(W*k_square(y_true - y_pred))
}





objective <- function(W=W, X=X, lambda = lambda, A1=A1, A2=A2, B1=B1, B2=B2){

  objval <- 0.5*sum(W*(X - A1%*%A2%*%X%*%B1%*%B2)^2)+lambda*norm(A1, 'F')^2 + lambda*norm(A2, 'F')^2+
    lambda*norm(B1, 'F')^2 + lambda*norm(B2, 'F')^2

  return(objval)

}

