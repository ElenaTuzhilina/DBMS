#' @title Principal curve metric scaling
#'
#' @description PCMS function calculates the principal curve metric scaling solution for some similarity matrix \eqn{Z} and a spline basis matrix \eqn{H}.
#' The optimal solution is found via minimizing the Frobenius norm \deqn{\|Z - S(X)\|^2_F}{||Z - S(X)||^2_F} w.r.t. \eqn{\Theta} subject to the smooth curve constraint \eqn{X = H\Theta}. 
#' Here \deqn{S(X) = XX^T}{S(X) = XX'} refers to the matrix of inner products.
#' The spatial coordinates of the resulting reconstruction are presented in \eqn{X}.
#'
#' @param Z a square symmetric matrix.
#' @param H a spline basis matrix, assumed to have orthonormal columns. If not orthonotmal, orthogonalization should be done via QR decomposition prior to running this function.
#' @return A list containing the PCMS problem solution:
#' \itemize{
#'   \item \code{Theta} -- the matrix of spline parameters.
#'   \item \code{X} -- the resulting conformation reconstruction.
#'   \item \code{loss} -- the resulting value of the PCMS objective.
#'   \item \code{rank} -- the resulting rank of the PCMS objective (can be \eqn{\leq 3}).
#' }
#'
#' @examples
#' data(C)
#' 
#' #transform contact counts to similarities
#' D = 1/(C+1)
#' Z = -D^2/2
#' Z = scale(Z, scale = FALSE, center = TRUE)
#' Z = t(scale(t(Z), scale = FALSE, center = TRUE))
#' 
#' #create spline basis matrix
#' H = splines::bs(1:ncol(C), df = 5)
#' 
#' #orthogonalize H using QR decomposition
#' H = qr.Q(qr(H))
#' 
#' #compute the PCMS solution
#' PCMS(Z, H)
#'
#' @export PCMS


PCMS = function(Z, H){
  ED3 = rARPACK::eigs_sym(t(H) %*% Z %*% H, 3, which = "LA")
  U3 = ED3$vectors
  d3 = ED3$values
  rank = sum(d3 > 1e-12)
  d3 = pmax(d3, 0)
  Theta = U3 %*% diag(sqrt(d3))
  X = H %*% Theta
  return(list(Theta = Theta, X = X, loss = loss_PCMS(X, Z), rank = rank))
}

loss_PCMS = function(X, Z){
  obj = (Z - X %*% t(X))^2
  return(mean(obj))
}