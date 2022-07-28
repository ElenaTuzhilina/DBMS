#' @title Weighted principal curve metric scaling
#'
#' @description WPCMS function calculates the weighted principal curve metric scaling solution for a contact matrix \eqn{Z} (e.g. representing a Hi-C contact matrix after some proper transformation applied) and a spline basis matrix \eqn{H}.
#' The optimal solution is found via minimizing the weighted Frobenius norm \deqn{ \|W * (Z - D^2(X) + \beta) \|^2_F }{|| W x (Z - D^2(X) + \beta) ||} w.r.t. \eqn{\Theta} subject to the smooth curve constraint \eqn{X = H\Theta}.
#' Here \eqn{D(X)} refers to the matrix of pairwise distances.
#' The spatial coordinates of the resulting reconstruction are presented in \eqn{X}.
#'
#' @param Z a square symmetric matrix.
#' @param H a spline basis matrix. By default assumed to have orthogonal columns. If not, orthogonalization should be done via QR decomposition.
#' @param W a weight matrix, same dimension as \code{Z}. By default equal weights \code{W = 1} are assumed.
#' @param Theta an initialization for spline basis coefficient matrix \eqn{Theta}. If \code{Theta = NULL} the default option \code{Theta = matrix(rnorm(ncol(H) * 3), ncol(H), 3)} is used, i.e. a random initialization is considered.
#' @param beta an initialization for intercept \eqn{beta}. If \code{beta = NULL} the default option \code{beta = -min(Z)} is used.
#' @param update_beta logical. If \code{update_beta = TRUE}, then the algorithm finds an optimal intercept value. If \code{update_beta = FALSE}, then the intercept is considered to be fixed and set to \code{beta}.
#' @param eps a positive convergence tolerance rate. Default value is \code{eps = 1e-6}.
#' @param maxiter an integer giving the maximal number of iterations. Default value is \code{maxiter = 100}.
#' @param verbose logical. If \code{verbose = TRUE}, the WPCMS loss after each iteration is printed.
#' @return A list containing the WPCMS problem solution:
#' \itemize{
#'   \item \code{Theta} -- the matrix of spline parameters.
#'   \item \code{X} -- the resulting conformation reconstruction.
#'   \item \code{beta} -- the resulting intercept value.
#'   \item \code{loss} -- the resulting value of the WPCMS loss.
#'   \item \code{iter} -- the total number of iterations.
#'   \item \code{plot} -- the list of WPCMS plots: 'loss' corresponds to loss vs. iteration plot, 'intercept' represents beta vs. iteration plot.
#' }
#'
#' @examples
#' data(C)
#' 
#' #transform contact counts to distances
#' Z = 1/(C+1)
#' 
#' #create spline basis matrix
#' H = splines::bs(1:ncol(C), df = 5)
#' 
#' #orthogonalize H using QR decomposition
#' H = qr.Q(qr(H))
#' 
#' #run WPCMS with equal weights; optimize intercept 
#' WPCMS(Z, H)$X
#' 
#' #run WPCMS with random weights; fixed intercept
#' W = matrix(runif(length(Z)), nrow(Z), ncol(Z))
#' WPCMS(Z, H, W, beta = 1, update_beta = FALSE)$X
#'
#' @export WPCMS


WPCMS = function(Z, H, W = matrix(1, nrow(Z), ncol(Z)), Theta = NULL, beta = NULL,
                update_beta = TRUE, eps = 1e-6, maxiter = 100, verbose = FALSE){
  #Initialize
  if(is.null(Theta)) Theta = matrix(rnorm(ncol(H) * 3), ncol(H), 3)
  if(is.null(beta)) beta = -min(Z)
  X = H %*% Theta
  X = scale(X, scale = FALSE, center = TRUE)
  obj = loss_WPCMS(X, Z, W, beta)
  
  delta = Inf
  iter = 0
  
  info = c(iter, 0, beta, obj, delta)
  info_names = c("iter", "rate", "beta", "objective", "delta")
  if(verbose){
    cat("\n\n--------------------[WPCMS loop]--------------------\n")
    paste("\n ", paste0(info_names, ": ", info))
  } 
  
  #Iterate
  while(delta > eps && iter < maxiter){
    iter = iter + 1
    obj0 = obj
    
    #Line search
    find = line_search_WPCMS(X, Z, W, H, beta, obj0)
    
    #Update solution
    if(find$rate > 0) Theta = find$Theta
    X = H %*% Theta
    X = scale(X, scale = FALSE, center = TRUE)
    
    #Update beta
    if(update_beta) beta = update_param_WPCMS(X, Z, W)
    
    #Update loss
    obj = loss_WPCMS(X, Z, W, beta)
    delta = abs((obj0 - obj)/obj0)
    info = rbind(info, c(iter, find$rate, beta, obj, delta))
    
    if(verbose) cat("\n ", paste(info_names, ":", info[iter + 1,]))
  }
  info = data.frame(info)
  names(info) = info_names
  info$df = ncol(H)
  return(list(Theta = Theta, X = X, beta = beta, info = info, iter = iter, loss = obj))
}

loss_WPCMS = function(X, Z, W, beta){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  obj = W * (Z + beta - D^2)^2
  return(mean(obj))
}

deriv_WPCMS = function(X, Z, W, beta){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  G = W * (Z + beta - D^2)
  G_plus = diag(rowSums(G))
  return(G - G_plus)
}

line_search_WPCMS = function(X, Z, W, H, beta, obj0){
  S = X %*% t(X)
  G = deriv_WPCMS(X, Z, W, beta)
  
  rate = 2
  obj = Inf
  rank = 0
  
  while(obj0 < obj || rank < 3){
    rate = rate / 2
    pcms = PCMS(S - rate * G, H)
    obj = loss_WPCMS(H %*% pcms$Theta, Z, W, beta)
    rank = pcms$rank
    if(rate < 1e-20) return(list(Theta = NA, rate = 0))
  }
  return(list(Theta = pcms$Theta, rate = rate))
}

update_param_WPCMS = function(X, Z, W){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  beta = -sum(W * (Z - D^2))/sum(W)
  return(beta)
}


