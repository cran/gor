% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10-grafos-TSP.R
\encoding{UTF-8}
\name{improve_tour_3opt}
\alias{improve_tour_3opt}
\title{Tour improving for a TSP using the 3-opt heuristic}
\usage{
improve_tour_3opt(d, n, C)
}
\arguments{
\item{d}{Distance matrix of the TSP.}

\item{n}{Number of vertices of the TSP complete graph.}

\item{C}{Starting tour to be improved.}
}
\value{
A list with two components: $tour contains a permutation
of the 1:n sequence representing the tour constructed by the
algorithm, $distance contains the value of the distance
covered by the tour.
}
\description{
3-opt heuristic tour-improving algorithm for the Traveling
Salesperson Problem
}
\details{
It applies the 3-opt algorithm to a starting tour of a TSP
instance until no further improvement can be found.  The tour
thus improved is a 3-opt local minimum.

The 3-opt algorithm consists of applying all possible
3-interchanges on the starting tour.  A 3-interchange removes
three non-indicent edges from the tour, leaving three pieces,
and combine them to form a new tour by interchanging the
endpoints in all possible ways and gluing them together by
adding the missing edges.
}
\examples{
## Regular example with obvious solution (minimum distance 32)
m <- 6   # Generate some points in the plane
z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
n <- nrow(z)
d <- compute_distance_matrix(z)
b <- build_tour_2tree(d, n)
b$distance    # Distance 38.43328
bi <- improve_tour_3opt(d, n, b$tour)
bi$distance   # Distance 32 (optimum)
plot_tour(z,b)
plot_tour(z,bi)

## Random points
set.seed(1)
n <- 15
z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
d <- compute_distance_matrix(z)
b <- build_tour_2tree(d, n)
b$distance    # Distance 45.788
bi <- improve_tour_3opt(d, n, b$tour)
bi$distance   # Distance 32.48669
plot_tour(z,b)
plot_tour(z,bi)

}
\seealso{
\link{improve_tour_2opt} improves a tour using the 2-opt
algorithm, \link{build_tour_nn_best} nearest neighbor heuristic,
\link{build_tour_2tree} double-tree heuristic,
\link{compute_tour_distance} computes tour distances,
\link{compute_distance_matrix} computes a distance matrix,
\link{plot_tour} plots a tour.
}
\author{
Cesar Asensio
}
