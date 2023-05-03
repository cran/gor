##' Generation of a system of fundamental cycles in a connected graph
##'      with respect of a given spanning tree.
##'
##' The routine loops through the edges of the graph outside the
##'     spanning tree (there are |E| - |V| + 1 of them); in each step,
##'     it adds an edge to the tree, thus closing a cycle, which has
##'     some "hair" in it in the form of dangling vertices.  Then all
##'     those dangling vertices are removed from the cycle (the "hair"
##'     is "shaven").
##' 
##' @title Generate fundamental cycles in a connected graph 
##' @param eT Spanning tree of the graph in edgelist representation,
##'     see [as_edgelist].
##' @param eG Graph in edgelist representation, see [as_edgelist].
##' @return A matrix with the fundamental cycles in its rows, in edge
##'     vector representation, that is, a binary vector with 1 if the
##'     edge belongs to the cycle and 0 otherwise.  This interpretation
##'     of the edge vectors of each fundamental cycle refers to the
##'     edgelist of the graph given in eG.
##' @author Cesar Asensio
##' @seealso [shave_cycle] shaves hairy cycles, [apply_incidence_map]
##'     applies the incidence map of a graph to an edge vector.
##' @examples
##' g <- make_graph("Dodecahedron")
##' n <- gorder(g)
##' b <- bfs(g, 1, father = TRUE)                 # BFS tree
##' T <- make_graph(rbind(b$father[2:n], 2:n), n) # Tree as igraph graph
##' eT <- as_edgelist(T)
##' eG <- as_edgelist(g)
##' C <- generate_fundamental_cycles(eT, eG)      # Fundamental cycles
##' mu <- gsize(g) - gorder(g) + 1                # Cyclomatic number
##' z <- layout_with_gem(g)
##' for (i in 1:mu) {                             # Cycle drawing
##'     c1 <- make_graph(t(eG[which(C[i,] == 1),]) , dir = FALSE)
##'     plot(g, layout = z)
##'     plot(c1, layout = z, add = TRUE, edge.color = "cyan4",
##'          edge.lty = "dashed", edge.width = 3)
##'     title(paste0("Cycle ", i, " of ", mu))
##'     #Sys.sleep(1) # Adjust time to see the cycles
##' }
##' 
##' @encoding UTF-8
##' @md 
##' @export 
generate_fundamental_cycles <- function(eT, eG) {
    qT <- nrow(eT)
    qG <- nrow(eG)
    v <- rep(0, qG)
    for (i in 1:qT) {
	inc <- which((eG[,1] == eT[i,1] & eG[,2] == eT[i,2]) |
		     (eG[,1] == eT[i,2] & eG[,2] == eT[i,1]))
	v[inc] <- 1
    }  # v: Edge-space vector associated to the tree given by eT
    mu <- qG - sum(v)   # Cyclomatic number of G
    C <- matrix(0, nrow = mu, ncol = qG)
    js <- which(v==0)
    for (i in 1:mu) {
	cy <- v
	cy[js[i]] <- 1
	cy <- shave_cycle(cy, eG)
	C[i,] <- cy
    }
    C
}

##' Removing dangling vertices of a cycle obtained by adding a single
##'     edge to a spanning tree.
##'
##' When generating a fundamental cycle in a graph, addition of a
##'     single edge to a spanning tree gives a "hairy" cycle, that is,
##'     a single cycle with some dangling branches of the tree.  This
##'     routine removes iteratively all leaves from this "hairy" tree
##'     until only a 2-regular, connected cycle remains, which is a
##'     fundamental cycle of the graph with respect the given spanning
##'     tree.
##' 
##' @title Shaving a hairy cycle
##' @param v Edge vector of the hairy cycle
##' @param eG Graph given as edgelist, see [as_edgelist]
##' @return Edge vector of the shaven cycle, to be interpreted with
##'     respect to the edgelist eG.
##' @author Cesar Asensio
##' @seealso [generate_fundamental_cycles] generates the edge vectors
##'     of a system of fundamental cycles of a graph,
##'     [apply_incidence_map] applies the incidence map of a graph to
##'     an edge vector.
##' @examples
##' ## It is used as a subroutine in [generate_fundamental_cycles].
##' 
##' @encoding UTF-8
##' @md 
##' @export 
shave_cycle <- function(v, eG) {
    d <- apply_incidence_map(eG, v)
    h <- which(d==1)   # Leaves: vertices of degree one
    va <- v
    while (length(h) > 0) {
	for (i in 1:length(h)) {
	    eh <- which(eG[,1]==h[i] | eG[,2]==h[i])  # Leaf edges
	    va[eh] <- 0
	}
	d <- apply_incidence_map(eG,va)  # Recomputing degrees...
	h <- which(d==1)                 # ...and leaves
    }
    va  # When there is no more leaves a shaven cycle remains!
}

##' Apply incidence map of a graph to an edge vector.  It uses the
##'     edgelist of the graph instead of the incidence matrix.
##'
##' The incidence map is the linear transformation from the edge
##'     vector space to the vertex vector space of a graph associating
##'     to each edge its incident vertices.  It is customarily
##'     represented by the incidence matrix, which is a very large
##'     matrix for large graphs; for this reason it not efficient to
##'     use directly the incidence matrix.  This function uses the
##'     edgelist of the graph as returned by the [as_edgelist]
##'     function to compute the result of the incidence map on an edge
##'     vector, which is interpreted with respect to the same
##'     edgelist.
##' 
##' @title Apply incidence map of a graph to an edge vector
##' @param eG Graph in edgelist representation, see [as_edgelist].
##' @param v Edge vector to which the incidence map will be applied.
##' @return A vertex vector, having the degree of each vertex in the
##'     subgraph specified by the edge vector.
##' @author Cesar Asensio
##' @seealso [shave_cycle], for shaving hairy cycles, which makes use
##'     of this routine, and [generate_fundamental_cycles], using the
##'     former.
##' @examples 
##' g <- make_graph("Dodecahedron")
##' eG <- as_edgelist(g)
##' set.seed(1)
##' v <- sample(0:1, gsize(g), replace = TRUE) # Random edge vector
##' apply_incidence_map(eG, v) # 1 1 0 1 2 0 1 1 3 2 0 1 1 1 1 1 0 0 1 2
##' ## Plotting the associated subgraph
##' h <- make_graph(t(eG[v==1,]))
##' z <- layout_with_gem(g)
##' plot(g, layout = z)
##' plot(h, layout = z, add = TRUE, edge.color = "red3", edge.width = 3)
##' 
##' @encoding UTF-8
##' @md 
##' @export 
apply_incidence_map <- function(eG, v) {
    chi <- eG[v==1,]  # Edges transported by the edge vector "v"
    n <- max(eG)      # Graph order
    deg <- rep(0,n)   # Degree vector to be filled
    for (i in 1:n) {  # Computing degrees a vertex at a time
	deg[i] <- sum(chi==i)
    }
    deg
}
