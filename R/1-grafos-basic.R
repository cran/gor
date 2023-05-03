##' @title Constructing an Eulerian Cycle
##'
##' @description It finds an \strong{Eulerian cycle} using a
##'     \eqn{O(q)} algorithm where \eqn{q} is the number of edges of
##'     the graph.
##'
##' @details
##' Recursive algorithm for undirected graphs.
##' The input graph should be Eulerian and connected.
##'
##' @section Disclaimer:
##' This function is part of the subject "Graphs and Network
##'     Optimization".  It is designed for teaching purposes only, and
##'     not for production.
##'
##' * It is introduced in the "Connectivity" section.
##'
##' * It is used as a subroutine in the "TSP" section.
##'
##' @references
##' \emph{Korte, Vygen: Combinatorial Optimization (Springer)} sec 2.4. 
##' 
##' @seealso [build_tour_2tree] double-tree algorithm.
##' 
##' @param G Eulerian and connected graph 
##' @param v Any vertex as starting point of the cycle
##' @return A two-element list: the $walk component is a
##'     \eqn{q\times2} matrix edgelist, and the $graph component is
##'     the input graph with no edges; it is used in intermediate
##'     steps, when the function calls itself.
##' @author Cesar Asensio (2021)
##' @examples
##' library(igraph)
##' find_euler(make_full_graph(5), 1)$walk   # Walk of length 10
##' 
##' @md
##' @encoding UTF-8
##' @export
find_euler <- function(G,v) {
    q <- gsize(G)
    W <- matrix(0,nrow=q,ncol=2)       # (1)
    i <- 1
    x <- v                             # (1)
    while(TRUE) {
	dx <- incident(G,x)
	if (length(dx)==0) {           # (2)
	    kp1 <- which(W[,1]==0)[1]
	    if (is.na(kp1)) { k <- q } else { k <- kp1-1 }
	    if (k==0) {
		return(list(walk=c(),graph=G))
	    }
	    for (j in 2:k) {           # (4)
		ce <- find_euler(G,W[j,1])
		G <- ce$graph
		assign(paste0("W",j),ce$walk)
		if (gsize(G)==0) break
	    }
	    break
	} else {                       # (2)
	    e <- dx[1]
	    e2 <- ends(G,e)
	    if (e2[1]==x) {            # (3)
		x <- e2[2]
		W[i,] <- e2
	    } else {
		x <- e2[1]
		W[i,] <- rev(e2)
	    }
	    i <- i+1
	    G <- delete_edges(G,e)
	}
    }
    Wt <- W[1,]                        # (5)
    for (j in 2:k) {
	if (exists(paste0("W",j))) {
	    Wt <- rbind(Wt,get(paste0("W",j)),deparse.level=0)
	}
	Wt <- rbind(Wt,W[j,],deparse.level=0)
    }
    list(walk=Wt,graph=G)
}

## Test:
## find_euler(make_full_graph(5),1)$walk
## rbind(c(1, 2, 3, 4, 5, 3, 1, 4, 2, 5),
##       c(2, 3, 4, 5, 3, 1, 4, 2, 5, 1))
