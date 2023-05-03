##' Greedy algorithm for coloring the vertices of a graph
##'
##' "Colors" are integers from 1 to the order of the graph to be
##'     colored.  The greedy strategy assigns to each vertex v the least
##'     color not assigned to the neighbors of v.
##' 
##' @title Greedy coloring of a graph
##' @param g Graph to be colored
##' @param ord Specified vertex ordering or NULL if natural vertex
##'     ordering is preferred
##' @param ran Choose random vertex ordering; it defaults to FALSE.
##'     It is ignored if ord is non-NULL
##' @return Vertex colors in a vector, with component i being the
##'     (integer) color of vertex i. 
##' @author Cesar Asensio
##' @examples
##' library(igraph)
##' g <- make_graph("Petersen")
##' cg <- color_graph_greedy(g)
##' plot(g, vertex.color = rainbow(6)[cg])
##' max(cg)         # = 3: Number of colors used by the algorithm
##' sum(g[cg == 1, cg == 1])  # = 0: Color classes are stable sets
##' 
##' g <- make_graph("Dodecahedron")
##' cg <- color_graph_greedy(g)
##' plot(g, vertex.color = rainbow(6)[cg])
##' max(cg)   # = 4: Number of colors used by the algorithm
##' sum(g[cg == 1, cg == 1])  # = 0: Color classes are stable sets
##'
##' ## However, the dodecahedron has a 3-coloring:
##' cdod <- rep(1, 20)
##' cdod[c(1,3,7,12,13,20)] <- 2
##' cdod[c(5,8,9,11,17,19)] <- 3
##' plot(g, vertex.color = rainbow(6)[cdod])
##' sum(g[cdod == 1, cdod == 1]) # = 0
##' sum(g[cdod == 2, cdod == 2]) # = 0
##' sum(g[cdod == 3, cdod == 3]) # = 0
##'
##' ## Some vertex orderings can use less colors:
##' cg <- color_graph_greedy(g, ord = 20:1)
##' plot(g, vertex.color = rainbow(6)[cg])
##' max(cg)         # = 3: Number of colors used by the algorithm
##' 
##' @encoding UTF-8
##' @md 
##' @export 
color_graph_greedy <- function(g, ord = NULL, ran = FALSE) {
    n <- gorder(g)
    C <- 1:n                              # Available colors
    if (is.null(ord)) {
	if (ran) {
	    V <- sample(1:n, n)           # Use random vertex ordering
	} else {
	    V <- 1:n                      # Use natural vertex ordering
	}
    } else {
	V <- ord                          # Use specified vertex ordering
    }
    col <- rep(NA,n);                     # All vertices uncolored
    for (v in V) {
	N <- neighbors(g, v)              # Neighbors of vertex v
	col[v] <- min(setdiff(C, col[N])) # Least unused color
    }
    col                                   # Final colors
}
