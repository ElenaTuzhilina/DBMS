#' @title Distribution-based metric scaling
#'
#' @description DBMS function calculates the Distribution-based metric scaling solution for a contact matrix \eqn{C} and a spline basis matrix \eqn{H}.
#' The optimal solution is found via minimizing the negative log-likelihood for a chosen probabilistic model \eqn{C~f(X, \Omega)} w.r.t. 3D conformation \eqn{X} and nuisance parameters \eqn{\Omega} subject to the smooth curve constraint \eqn{X = H\Theta}.
#' The solution can be calculated via iterating the second order approximation of the objective (outer epoch loop) and applying WPCMS to optimize the obtained quadratic approximation (inner iteration loop).
#' The spatial coordinates of the resulting reconstruction are presented in \eqn{X}.
#' There are four models for contact counts available (however, one can implement their own model):
#' \itemize{
#'   \item Poisson model \eqn{C~Pois(X, \beta)};
#'   \item Hurdle Poisson model \eqn{C~HPois(X, \beta, p)}, where \eqn{p} is the mixture weight;
#'   \item zero-inflated Poisson model \eqn{C~ZIPois(X, \beta, p)}, where \eqn{p} is the probability of extra zeros;
#'   \item negative binomial \eqn{C~NB(X, \beta, r)}, where \eqn{r} is the overdispersion parameter.
#' }
#' For each model the mean parameter \eqn{\Lambda} is linked to the 3D spatial structure via the log-link \deqn{\log(\Lambda) = -D^2(X) + \beta.}{log(\Lambda) = -D^2(X) + \beta.}
#' Here \eqn{\beta} is the intercept and \eqn{D(X)} refers to the matrix of pairwise distances between the genomic loci.
#' 
#' @param C a square symmetric matrix representing a contact matrix.
#' @param H a spline basis matrix, assumed to have orthonormal columns. If not orthonotmal, orthogonalization should be done via QR decomposition prior to running this function.
#' @param method a probabilistic model used for the contact counts. By default, \code{method = Pois}; other options available are \code{HPois}, \code{ZIPois} and \code{NB}.
#' @param Theta initialization for spline basis coefficient matrix \eqn{\Theta}. By default, \code{Theta = NULL}, so random initialization is done for this parameter.
#' @param param a list containing initialization for the nuisance parameters. If \code{param = NULL}, then the default option is used: \itemize{\item for all methods \code{beta = log(mean(C[C>0]))}, i.e. logarithm of the average of non-zero counts; \item for the \code{HPois} and \code{ZIPois} methods \code{p = mean(C == 0)}, i.e. the proportion of zero counts; \item for the \code{NB} method \code{r = 1}.}
#' @param update_param logical. If \code{update_param = TRUE}, then the algorithm finds optimal nuisance parameter values using \code{param} values as initialization. If \code{update_param = FALSE}, then the parameters are considered to be fixed at \code{param} values.
#' @param eps_wpcms,eps_dbms positive convergence tolerance rates for the inner and outer loops. Default values are \code{eps_wpcms = 1e-6} and \code{eps_dbms = 1e-6}.
#' @param maxiter,maxepoch integers corresponding to maximum number of iterations for the inner loop and maximum number of epochs for the outer loop. Default values are \code{maxiter = 100} and \code{maxepoch = 100}.
#' @param verbose_wpcms,verbose_dbms logical. If \code{verbose_wpcms = TRUE}, the WPCMS loss after each iteration in the inner loop is printed. If \code{verbose_dbms = TRUE} the DBMS loss after each epoch in the outer loop is printed. Default values are \code{verbose_wpcms = FALSE} and \code{verbose_dbms = FALSE}.
#' @return A list containing the DBMS problem solution:
#' \itemize{
#'   \item \code{Theta} -- the matrix of spline parameters.
#'   \item \code{X} -- the resulting conformation reconstruction.
#'   \item \code{param} -- the resulting list of the nuisance parameters.
#'   \item \code{info} -- table with detailed convergence information that can be used ofr plotting.
#'   \item \code{epoch} -- the total number of epochs.
#'   \item \code{iter_total} -- the total number of iterations.
#'   \item \code{loss} -- the resulting value of the DBMS objective.
#' }
#'
#' @examples
#' data(C)
#'
#' #create spline basis matrix
#' H = splines::bs(1:ncol(C), df = 5)
#' 
#' #orthogonalize H using QR decomposition
#' H = qr.Q(qr(H))
#' 
#' #run the PoisMS approach; optimize beta; print the epoch and iteration convergence progress
#' DBMS(C, H, verbose_wpcms = TRUE, verbose_dbms = TRUE)
#' 
#' #run the HPoiMS approach; fix beta and p to log of the average contact count and 0.5, respectively
#' param = list(beta = log(mean(C)), p = 0.5)
#' DBMS(C, H, method = HPois, param = param, update_param = FALSE)
#' 
#' #run the ZIPoisMS approach; optimize beta and p initializing them to param; print the outer loop convergence progress
#' DBMS(C, H, method = HPois, param = param, verbose_dbms = TRUE)
#' 
#' #run the NBMS approach; fix beta and r
#' param = list(beta = log(mean(C)), r = 1)
#' DBMS(C, H, method = NB, param = param, update_param = FALSE)
#' 
#' @export DBMS


DBMS = function(C, H, Theta = NULL, method = Pois, param = NULL, update_param = TRUE, eps_wpcms = 1e-6, maxiter = 100, verbose_wpcms = FALSE, eps_dbms = 1e-6, maxepoch = 100, verbose_dbms = FALSE){
  #Initialize
  if(is.null(Theta)) Theta = matrix(rnorm(ncol(H) * 3), ncol(H), 3)
  method_name = deparse(substitute(method))
  if(is.null(param)) param = list()
  if(is.null(param$beta)) param$beta = log(mean(C[C>0]))
  if(method_name %in% c("HPois" , "ZIPois") & is.null(param$p)) param$p = mean(C == 0)
  if(method_name == "NB" & is.null(param$r)) param$r = 1
  
  X = H %*% Theta
  X = scale(X, scale = FALSE, center = TRUE)
  obj = method$loss(X, C, param)
  
  delta = Inf
  deltaX = Inf
  iter_total = 0
  epoch = 0
  
  info = c(epoch, 0, iter_total, unlist(param), obj, delta, deltaX)
  info_names = c("epoch", "rate", "iter_total", names(param), "loss", "delta", "deltaX")
  if(verbose_dbms){
    cat("\n\n====================[DBMS loop]====================\n")
    cat("\n  ", paste0(info_names, ": ", info))
  }  
  
  #Iterate
  while(deltaX > eps_dbms && epoch < maxepoch){
    epoch = epoch + 1
    obj0 = obj
    X0 = X
    
    #SOA
    soa = SOA(X, method$deriv(X, C, param, "D2"))
    W = soa$W
    Z = soa$Z
    
    #WPCMS (beta = 0)
    wpcms = WPCMS(Z, H, W, Theta, 0, FALSE, eps_wpcms, maxiter, verbose_wpcms)
    iter_total = iter_total + wpcms$iter
    
    #Line search
    find = line_search(X, wpcms$X, C, H, param, obj, method$loss)
    
    #Update solution
    if(find$rate > 0) Theta = find$Theta
    X = H %*% Theta
    X = scale(X, scale = FALSE, center = TRUE)
    
    #Update parameters
    if(update_param) param = method$update_param(X, C, param, method, eps_dbms, verbose_dbms)
    
    #Update loss
    obj = method$loss(X, C, param)
    delta = abs((obj0 - obj)/obj0)
    deltaX = vegan::procrustes(X0, X, scale = TRUE, symmetric = TRUE)$ss
    info = rbind(info, c(epoch, find$rate, iter_total, unlist(param), obj, delta, deltaX))
     
    if(verbose_dbms){
      cat("\n====================[DBMS loop]====================\n")
      cat("\n  ", paste0(info_names, ": ", info[epoch + 1,]))
    }  
  }
  info = data.frame(info)
  names(info) = info_names
  info$df = ncol(H)
  return(list(Theta = Theta, X = X, param = param, info = info, epoch = epoch, iter_total = iter_total, loss = obj))
}

############################################

line_search = function(X0, X, C, H, param, obj0, loss){
  S0 = X0 %*% t(X0)
  S = X %*% t(X)
  
  rate = 2
  obj = Inf
  rank = 0
  
  while(obj0 < obj || rank < 3){
    rate = rate / 2
    pcms = PCMS((1 - rate) * S0 + rate * S, H)
    obj = loss(H %*% pcms$Theta, C, param)
    rank = pcms$rank
    if(rate < 1e-20) return(list(Theta = NA, rate = 0))
  }
  return(list(Theta = pcms$Theta, rate = rate))
}

############################################

SOA = function(X, deriv){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  G = deriv$G
  H = deriv$H
  W = H
  diag(W) = 0
  W = W/max(W)
  Z = D^2 - G/H
  Z[W == 0] = 0
  return(list(Z = Z, W = W))
}

############################################

newton = function(param, method, eps, verbose){
  obj = method$loss(param)
  delta = Inf
  while(delta > eps){
    der = method$deriv(param)
    rate = 1
    while(!method$cond(param - rate * der$G/der$H)) rate = rate/2
    param = param - rate * der$G/der$H
    obj0 = obj
    obj = method$loss(param)
    delta = abs((obj0 - obj)/obj0)
    if(verbose) cat("\n  param:", param, "rate:", rate, "loss:", obj, "delta:", delta)
  }
  cat("\n")
  return(param)
}

sgd = function(param, method, eps, verbose){
  obj = method$loss(param)
  delta = Inf
  while(delta > eps){
    der = method$deriv(param)
    rate = 1
    while(!method$cond(param - rate * der$G)) rate = rate/2
    param = param - rate * der$G
    obj0 = obj
    obj = method$loss(param)
    delta = abs((obj0 - obj)/obj0)
    if(verbose) cat("\n  param:", param, "rate:", rate, "loss:", obj, "delta:", delta)
  }
  cat("\n")
  return(param)
}

############################################

Pois = list(
  loss = function(X, C, param){
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    obj = L - C * logL
    return(mean(obj))
  },
  
  deriv = function(X, C, param, var){
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    H = L
    G = C - L
    if(var == "D2") return(list(G = G, H = H))
  },
  
  update_param = function(X, C, param, method, eps, verbose){
    if(verbose)  cat("\n\n--------------------[update beta]--------------------\n")
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    param$beta = log(sum(C)/sum(exp(-D^2)))
    return(param)
  }
)

HPois = list(
  loss = function(X, C, param){
    beta = param$beta
    p = param$p
    
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    
    obj = matrix(0, nrow(C), ncol(C))
    obj[C == 0] = -log(p)
    obj[C != 0] = -log(1 - p) + (L - C * logL + log(1 - exp(-L)))[C != 0]
    
    return(mean(obj))
  },
  
  deriv = function(X, C, param, var){
    beta = param$beta
    
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    M = L/(1 - exp(-L))
    
    G = H = matrix(0, nrow(C), ncol(C))
    G[C != 0] = (C - M)[C != 0]
    H[C != 0] =  (M * (1 - M * exp(-L)))[C != 0]
    
    if(var == "D2") return(list(G = G, H = H))
    if(var == "beta") return(list(G = -sum(G), H = sum(H)))
  },
  
  update_param = function(X, C, param, method, eps, verbose){
    if(verbose)  cat("\n\n--------------------[update beta]--------------------\n")
    method_newton = list(loss = function(beta) method$loss(X, C, list(p = param$p, beta = beta)),
                         deriv = function(beta) method$deriv(X, C, list(p = param$p, beta = beta), "beta"),
                         cond = function(beta) TRUE)
    param$beta = newton(param$beta, method_newton, eps, verbose)
    
    if(verbose)  cat("\n--------------------[update p]--------------------\n")
    param$p = mean(C == 0)
    return(param)
  }
)


ZIPois = list(
  loss = function(X, C, param){
    p = param$p
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    a = p /(1 - p)
    obj = matrix(0, nrow(C), ncol(C))
    obj[C == 0] = -log(a + exp(-L))[C == 0] - log(1 - p)
    obj[C != 0] = (L - C * logL)[C != 0] - log(1 - p)
    return(mean(obj))
  },
  deriv = function(X, C, param, var){
    p = param$p
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    expL = exp(L)
    a = p/(1 - p)  
    
    G = H = matrix(0, nrow(C), ncol(C))
    
    if(var %in% c("D2", "beta")){
      G[C == 0] = (-L / (a * expL + 1))[C == 0]
      G[C != 0] = (C - L)[C != 0]
      H[C == 0] =  (L * (a * expL * (1 - L) + 1)/(a * expL + 1)^2)[C == 0]
      H[C != 0] = L[C != 0]
      if(var == "D2") return(list(G = G, H = H))
      else return(list(G = -sum(G), H = sum(H)))
    }
    if(var == "p"){
      G[C == 0] = (-(expL - 1) / (p * (expL - 1) + 1))[C == 0]
      G[C != 0] = 1/(1 - p)
      H[C == 0] =  ((expL - 1)^2 / (p * (expL - 1) + 1)^2)[C == 0]
      H[C != 0] = 1/(1 - p)^2
      return(list(G = sum(G), H = sum(H)))
    }
  },
  update_param = function(X, C, param, method, eps, verbose){
    if(verbose)  cat("\n\n--------------------[update beta]--------------------\n")
    method_newton = list(loss = function(beta) method$loss(X, C, list(p = param$p, beta = beta)),
                         deriv = function(beta) method$deriv(X, C, list(p = param$p, beta = beta), "beta"), 
                         cond = function(beta) TRUE)
    param$beta = newton(param$beta, method_newton, eps, verbose)
    
    if(verbose)  cat("\n--------------------[update p]--------------------\n")
    method_newton = list(loss = function(p) method$loss(X, C, list(p = p, beta = param$beta)),
                         deriv = function(p) method$deriv(X, C, list(p = p, beta = param$beta), "p"), 
                         cond = function(p) p > 0 & p < 1)
    param$p = newton(param$p, method_newton, eps, verbose)
    return(param)
  }
)

NB = list(
  loss = function(X, C, param){
    r = param$r
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    obj =  lgamma(r) - lgamma(C + r) - r * log(r)  + (C + r) * log(r + L) - C * logL
    return(mean(obj))
  },
  
  deriv = function(X, C, param, var){
    r = param$r
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    if(var %in% c("D2", "beta")){
      G = r * (C - L) / (L + r)
      H = r * L * (C + r)/(L + r)^2
      if(var == "D2") return(list(G = G, H = H))
      if(var == "beta") return(list(G = -sum(G), H = sum(H)))
    }
    if(var == "r"){
      G = digamma(r) - digamma(C + r) - log(r) + log(L + r) + (C - L)/(L + r)
      H = trigamma(r) - trigamma(C + r) - 1/r + 1/(L + r) - (C - L)/(L + r)^2
      return(list(G = sum(G), H = sum(H)))
    }
  },
  
  update_param = function(X, C, param, method, eps, verbose){
    method_newton = list(loss = function(beta) method$loss(X, C, list(r = param$r, beta = beta)),
                         deriv = function(beta) method$deriv(X, C, list(r = param$r, beta = beta), "beta"),
                         cond = function(beta) TRUE)
    if(verbose)  cat("\n\n--------------------[update beta]--------------------\n")
    param$beta = newton(param$beta, method_newton, eps, verbose)
    
    method_newton = list(loss = function(r) method$loss(X, C, list(r = r, beta = param$beta)),
                         deriv = function(r) method$deriv(X, C, list(r = r, beta = param$beta), "r"), 
                         cond = function(r) r > 0)
    if(verbose)  cat("\n--------------------[update r]--------------------\n")
    param$r = sgd(param$r, method_newton, eps, verbose)
    return(param)
  }
)


