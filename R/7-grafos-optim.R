##' Dijkstra's algorithm finding the sorthest paths from a root vertex
##' to the remaining vertices of a graph using a spanning tree
##'
##' An implementation of Dijkstra's algorithm.
##' @title Dijkstra' algorithm for shortest paths
##' @param g An igraph Graph
##' @param d Weights (lengths) of the edges of \eqn{g}
##' @param r Starting vertex --- root of the output tree
##' @return A list with components: $tree, which is a sequence of
##'     pairs of vertices parent-son; $distances, which is a
##'     \eqn{2\times n} matrix with distances from the root vertex to
##'     the remaining vertices, and $parents, which contains the
##'     parent of each vertex in the tree, except for the root which
##'     has no parent, so its entry is NA.
##' @author Cesar Asensio
##' @seealso [shortest_paths] in the [igraph] package.
##' @examples
##' library(igraph)
##' g <- make_graph("Frucht")
##' n <- gorder(g)
##' set.seed(1);
##' d <- matrix(round(runif(n^2, min = 1, max = 10)), nrow = n)  # Distances
##' d <- d + t(d); for (i in 1:n) { d[i,i] <- 0 }          # Distance matrix
##' Td <- dijk(g, d, r = 1)
##' Td$distances
##' Td$parents
##' gTd <- make_graph(Td$tree, n = gorder(g))   # igraph tree
##' Eg <- as_edgelist(g)
##' dl <- c()   # We convert the matrix in a list:
##' for (e in 1:nrow(Eg)) { dl <- c(dl, d[Eg[e,1], Eg[e,2]]) }
##' z <- layout_with_kk(g)
##' plot(g, layout = z, edge.label = dl)
##' plot(gTd, layout = z, edge.color = "red3", add = TRUE)
##' 
##' @encoding UTF-8
##' @md 
##' @export 
dijk <- function(g, d, r = 1) { # Dijkstra's algorithm
    n <- gorder(g)
    Vg <- V(g)
    U <- c()                # Will contain vertices added to the tree
    L <- rbind(1:n,         # Vertex labels...
	       rep(Inf,n))  # ..and infinite initial distances...
    L[2,r] <- 0             # ...except for the root, which is 0
    p <- rep(NA,n)          # Predecessors
    while (length(U) < n) {            # When |U| = n we will stop
	Uc <- setdiff(Vg,U)            # Non-added vertices
	k <- which.min(L[2,Uc])        # Take the minimum distance vertex...
	u <- L[1,Uc][k]                # ...which is this...
	U <- c(U,u)                    # ...and add it to the tree
	UciNu <- intersect(Uc,neighbors(g,u))  # Non-added vertices
	for (w in UciNu) {             # Update L:
	    if (L[2,w] > L[2,u] + d[u,w]) {
		L[2,w] <- L[2,u] + d[u,w]
		p[w] <- u
	    }
	}
    }
    T <- c()              # Tree
    for (v in 1:n) {
	if (v != r) { T <- c(T,p[v],v) }
    }
    list(tree = T, distances = L, parents = p)
}
