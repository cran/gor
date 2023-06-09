% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11-grafos-NPdif.R
\encoding{UTF-8}
\name{build_cover_greedy}
\alias{build_cover_greedy}
\title{Greedy algorithm for vertex cover in a graph}
\usage{
build_cover_greedy(G)
}
\arguments{
\item{G}{Graph}
}
\value{
A list with two components: $set contains the cover, $size
contains the number of vertices of the cover.
}
\description{
This routine uses a greedy algorithm to build a cover selecting
the highest degree vertex first and removing its incident
edges.
}
\details{
This algorithm builds a vertex cover since no edge remains to be
covered when it returns.  However, it is no guaranteed that
the cover found by this algorithm has minimum cardinality.
}
\examples{
## Example with known cover
K25 <- make_full_graph(25)   # Cover of size 24
X0 <- build_cover_greedy(K25)
X0$size  # 24
plot_cover(X0, K25)
plot_cover(list(set = c(1,2), size = 2), K25)

## Vertex-cover of a random graph
set.seed(1)
n <- 25
g <- sample_gnp(n, p=0.25)
X1 <- build_cover_greedy(g)
X1$size   # 17
plot_cover(X1, g)

}
\references{
Korte, Vygen \emph{Combinatorial Optimization. Theory
and Algorithms.}
}
\seealso{
\link{is_cover} checks if a vertex subset is a vertex cover,
\link{build_cover_approx} builds a cover using a 2-approximation
algorithm, \link{improve_cover_flip} improves a cover using local
search, \link{search_cover_random} looks for a random cover of
fixed size, \link{search_cover_ants} looks for a random cover using
a version of the ant-colony optimization heuristic,
\link{find_cover_BB} finds covers using a branch-and-bound
technique, \link{plot_cover} plots a cover.
}
\author{
Cesar Asensio
}
