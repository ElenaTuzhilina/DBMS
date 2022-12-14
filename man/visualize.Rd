% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{visualize}
\alias{visualize}
\title{Reconstruction vizualization}
\usage{
visualize(X, type = "projection", index = 1:nrow(X), title = NULL)
}
\arguments{
\item{X}{a matrix representing spatial coordinates of resulting chromatin reconstruction.}

\item{type}{the type of plot returned. Set \code{type = 'projection'} and \code{type = '3D'} to output the projection and 3D model of chromatin conformation reconstruction, respectively.}

\item{index}{points where spline basis is evaluated; each corresponds to a particular genomic loci.}

\item{title}{optional, adds title to the plot. Default value \code{title = NULL}.}
}
\value{
Reconstruction plots.
}
\description{
This function allows to plot reconstructed chromatin conformation \emph{X} and corresponding contact matrix approximation \emph{XX'}.
}
\examples{
data(C)

#transform contact counts to distances
Z = 1/(C+1)

#create spline basis matrix
H = splines::bs(1:ncol(C), df = 5)

#orthogonalize H using QR decomposition
H = qr.Q(qr(H))

#run WPCMS with equal weights 
wpcms = WPCMS(Z, H)

#plot projection of reconstructed chromatin conformation
visualize(wpcms$X, type = 'projection')

#plot 3D model of reconstructed chromatin conformation
visualize(wpcms$X, type = '3D')

}
