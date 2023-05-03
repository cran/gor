##' Computation of the depth-first tree search in an undirected graph.
##'
##' Starting from a root vertex, the tree is grown by adding neighbors
##'     of the last vertex added to the tree.  In this way, the tree
##'     has many levels and few branches and leaves.  When the tree
##'     cannot grow further, it backtracks to previously added
##'     vertices with neighbors outside the tree, adding them until
##'     none is left.
##' 
##' @title Depth-first search tree
##' @param g Graph
##' @param r Root: Starting vertex growing the tree.
##' @return A directed spanning subgraph of g containing the edges of
##'     the DFS tree.
##' @author Cesar Asensio
##' @seealso [bfs_tree] breadth-first search tree; [bfs] and [dfs] in
##'     the igraph package.
##' @examples
##' g <- make_graph("Frucht")
##' T <- dfs_tree(g, 5)  # Root at v = 5
##' z <- layout_with_gem(g)
##' plot(g, layout = z, main = "Depth-first search tree")
##' plot(T, layout = z, add = TRUE, edge.color = "cyan4", edge.width = 2)
##' plot(T, layout = layout_as_tree(T))
##' 
##' @encoding UTF-8
##' @md 
##' @export 
dfs_tree <- function(g, r) {
    Q <- R <- c(r)
    T <- c()
    while(length(Q) > 0) {
	v <- Q[1]
	N <- setdiff(neighbors(g, v), R)
	if (length(N)==0) {
	    Q <- setdiff(Q, v)
	    next
	} else { w <- N[1] }
	Q <- c(w, Q)
	R <- c(R, w)
	T <- c(T, v, w)
    }
    make_graph(T)
}

##' Computation of the breadth-first tree search in an undirected graph.
##'
##' Starting from a root vertex, the tree is grown by adding neighbors
##'     of the first vertex added to the tree until no more neighbors
##'     are left; then it passes to another vertex with neighbors
##'     outside the tree.  In this way, the tree has few levels and
##'     many branches and leaves.
##' 
##' @title Breadth-first search tree
##' @param g Graph
##' @param r Root: Starting vertex growing the tree.
##' @return A directed spanning subgraph of g containing the edges of
##'     the BFS tree.
##' @author Cesar Asensio
##' @examples 
##' g <- make_graph("Frucht")
##' T <- bfs_tree(g, 2)  # Root at v = 2
##' z <- layout_with_gem(g)
##' plot(g, layout = z, main = "Breadth-first search tree")
##' plot(T, layout = z, add = TRUE, edge.color = "cyan4", edge.width = 2)
##' plot(T, layout = layout_as_tree(T))
##' 
##' @encoding UTF-8
##' @md 
##' @export 
bfs_tree <- function(g, r) {
    Q <- R <- c(r)
    T <- c()
    while(length(Q) > 0) {
	v <- Q[length(Q)]
	N <- setdiff(neighbors(g, v), R)
	if (length(N) == 0) {
	    Q <- setdiff(Q, v)
	    next
	} else { w <- N[1] }
	Q <- c(w, Q)
	R <- c(R, w)
	T <- c(T, v, w)
    }
    make_graph(T)
}
