% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10-grafos-TSP.R
\encoding{UTF-8}
\name{compute_lower_bound_1tree}
\alias{compute_lower_bound_1tree}
\title{Computing the 1-tree lower bound for a TSP instance}
\usage{
compute_lower_bound_1tree(d, n, degree = FALSE)
}
\arguments{
\item{d}{Distance matrix.}

\item{n}{Number of vertices of the TSP complete graph.}

\item{degree}{Boolean: Should the routine return the degree
seguence of the internal minimum spanning tree? Defaults to
FALSE.}
}
\value{
The 1-tree lower bound -- A scalar if the optional
argument "degree" is FALSE.  Otherwise, a list with the
previous 1-tree lower bound in the $bound component and the
degree sequence of the internal minimum spanning tree in the
$degree component.
}
\description{
It computes the 1-tree lower bound for an optimum tour for a TSP instance.
}
\details{
It computes the 1-tree lower bound for an optimum tour for a TSP
instance from vertex 1.  Internally, it creates the graph Kn-{v1}
and invokes \link{mst} from package \link{igraph} to compute the minimum
weight spanning tree.  If optional argument "degree" is TRUE, it
returns the degree seguence of this internal minimum spanning
tree, which is very convenient when embedding this routine in the
Held-Karp lower bound estimation routine.
}
\examples{
m <- 10   # Generate some points in the plane
z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
n <- nrow(z)
d <- compute_distance_matrix(z)
b <- build_tour_2tree(d, n)
b$distance           # Distance 57.868
bi <- improve_tour_2opt(d, n, b$tour)
bi$distance          # Distance 48 (optimum)
compute_lower_bound_1tree(d,n)       # 45

## Random points
set.seed(1)
n <- 25
z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
d <- compute_distance_matrix(z)
compute_lower_bound_1tree(d,n)       # 31.4477
bn <- build_tour_nn_best(d,n)
b3 <- improve_tour_3opt(d,n,bn$tour)
b3$distance                          # 35.081

}
\seealso{
\link{improve_tour_2opt} tour improving using 2-opt,
\link{improve_tour_3opt} tour improving using 3-opt,
\link{compute_lower_bound_HK} for Held-Karp lower bound estimates.
}
\author{
Cesar Asensio
}
