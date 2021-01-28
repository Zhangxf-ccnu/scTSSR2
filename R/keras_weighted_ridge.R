keras_weighted_ridge <- function(Y=Y, X=X, k=k, W=W, epochs = 100, batch_size = 128,lambda=1.0,
                                   learning_rate = 0.0001, lasso_threshold = 0, verbose = TRUE){

  p <- ncol(X)
  n <- nrow(X)

  X <- array_reshape(X, dim = c(nrow(X), ncol(X)))
  Y <- array_reshape(Y, dim = c(nrow(Y), ncol(Y)))

  W <- array_reshape(W, dim = c(nrow(W), ncol(W)))



  network <- keras::keras_model_sequential()

  network %>%
    layer_dense(units = k, input_shape = p, activation = "linear", kernel_regularizer = function(weight_matrix)
      regularizer_define(weight_matrix = weight_matrix, lambda = lambda), use_bias = FALSE) %>%

    layer_dense(units = p, activation = "linear", kernel_regularizer = function(weight_matrix)
      regularizer_define(weight_matrix = weight_matrix, lambda = lambda), use_bias = FALSE)



  network %>% keras::compile(loss = function(y_true, y_pred) loss_define(W, y_true, y_pred),
                             optimizer = optimizer_rmsprop(lr = learning_rate), metrics = list("mean_squared_error"))



  history <- network %>% keras::fit(x = X, y = Y, epochs = epochs, batch_size = batch_size, validation_split = 0,
                                    verbose = verbose)



  weights = keras::get_weights(network)



  return(weights)

}






