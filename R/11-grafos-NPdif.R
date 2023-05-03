##' Plot of a vertex cover in a graph.
##'
##' It plots a graph, then superimposes a vertex cover in a different
##'     color.  It also draws the covered edges, to help in detecting
##'     non-covers by inspection.
##' 
##' @title Vertex cover plotting
##' @param X Cover to be plotted; an output list returned by some
##'     cover-building function, see below.
##' @param G Graph on which to superimpose the cover.
##' @return This function is called for its side effect of plotting.
##' @importFrom graphics legend
##' @importFrom graphics text
##' @author Cesar Asensio
##' @seealso [is_cover] checks if a vertex subset is a vertex cover,
##'     [build_cover_greedy] builds a cover using a greedy heuristic,
##'     [build_cover_approx] builds a cover using a 2-approximation
##'     algorithm, [improve_cover_flip] improves a cover using local
##'     search, [search_cover_random] looks for a random cover of
##'     fixed size, [search_cover_ants] looks for a random cover using
##'     a version of the ant-colony optimization heuristic,
##'     [find_cover_BB] finds covers using a branch-and-bound
##'     technique.
##' @examples 
##' set.seed(1)
##' g <- sample_gnp(25, p=0.25)   # Random graph
##' X1 <- build_cover_greedy(g)
##' plot_cover(X1, g)
##'
##' st <- 1:5   # Not a vertex cover
##' plot_cover(list(set = st, size = length(st)), g)  # See covered edges
##' 
##' @encoding UTF-8
##' @md 
##' @export 
plot_cover <- function(X, G) {
    n <- gorder(G)
    z <- layout_with_gem(G)
    plot(G, layout=z)
    eG <- as_edgelist(G)
    gcov <- make_empty_graph(n, directed=FALSE)
    for (v in X$set) {
	gcov <- gcov %>% add_edges(t(eG[as.numeric(incident(G,v)),]))
    }
    gcov <- simplify(gcov)
    vcol <- rep("pink", n)
    lcol <- rep("black", n)
    vcol[X$set] <- "red4"
    lcol[X$set] <- "white"
    plot(gcov, layout = z, add = TRUE, edge.color = "red3",
	 edge.width = 2, edge.lty = "dashed",
	 vertex.color = vcol, vertex.label.color = lcol)
    legend("topleft", legend = c("Vertex set", "Covered edges"),
	   lwd = c(NA,2), pch = c(21, NA), col = c("black", "red3"),
	   pt.bg = c("red4", NA), pt.cex = c(2.5, NA),
	   lty = c("solid", "dashed"))
    text(0, 1.2,
	 labels = paste0("Size = ", X$size, ", Uncovered edges = ",
			 gsize(G) - gsize(gcov)))
}


##' This routine uses a greedy algorithm to build a cover selecting
##'     the highest degree vertex first and removing its incident
##'     edges.
##'
##' This algorithm builds a vertex cover since no edge remains to be
##'     covered when it returns.  However, it is no guaranteed that
##'     the cover found by this algorithm has minimum cardinality.
##'
##' @title Greedy algorithm for vertex cover in a graph
##' @param G Graph
##' @return A list with two components: $set contains the cover, $size
##'     contains the number of vertices of the cover.
##' @author Cesar Asensio
##' @references Korte, Vygen \emph{Combinatorial Optimization. Theory
##'     and Algorithms.}
##' @seealso [is_cover] checks if a vertex subset is a vertex cover,
##'     [build_cover_approx] builds a cover using a 2-approximation
##'     algorithm, [improve_cover_flip] improves a cover using local
##'     search, [search_cover_random] looks for a random cover of
##'     fixed size, [search_cover_ants] looks for a random cover using
##'     a version of the ant-colony optimization heuristic,
##'     [find_cover_BB] finds covers using a branch-and-bound
##'     technique, [plot_cover] plots a cover.
##' @examples 
##' ## Example with known cover
##' K25 <- make_full_graph(25)   # Cover of size 24
##' X0 <- build_cover_greedy(K25)
##' X0$size  # 24
##' plot_cover(X0, K25)
##' plot_cover(list(set = c(1,2), size = 2), K25)
##'
##' ## Vertex-cover of a random graph
##' set.seed(1)
##' n <- 25
##' g <- sample_gnp(n, p=0.25)
##' X1 <- build_cover_greedy(g)
##' X1$size   # 17
##' plot_cover(X1, g)
##'
##' @encoding UTF-8
##' @md 
##' @export 
build_cover_greedy <- function(G) {
    n <- gorder(G)
    X <- rep(FALSE, n)
    h <- G
    while (gsize(h) > 0) {
	dv <- degree(h)
	v <- which.max(dv)
	h[v,] <- FALSE
	X[v] <- TRUE
    }
    list(set = as.numeric(V(G)[X]), size = sum(X))
}

##' Gavril's 2-approximation algorithm to build a vertex cover.
##'
##' This algorithm computes a maximal matching and takes the ends of
##'     the edges in the matching as a vertex cover.  No edge is
##'     uncovered by this vertex subset, or the matching would not be
##'     maximal; therefore, the vertex set thus found is indeed a
##'     vertex cover.
##'
##' Since no vertex can be incident to two edges of a matching M, at
##'     least |M| vertices are needed to cover the edges of the
##'     matching; thus, any vertex cover X should satisfy |X| >= |M|.
##'     Moreover, the vertices incident to the matching are always a
##'     vertex cover, which implies that, if X* is a vertex cover of
##'     minimum sise, |X*| <= 2|M|.
##' 
##' @title 2-approximation algorithm for vertex cover
##' @param G Graph
##' @return A list with two components: $set contains the cover, $size
##'     contains the number of vertices of the cover.
##' @author Cesar Asensio
##' @seealso [is_cover] checks if a vertex subset is a vertex cover,
##'     [build_cover_greedy] builds a cover using a greedy heuristic,
##'     [improve_cover_flip] improves a cover using local
##'     search, [search_cover_random] looks for a random cover of
##'     fixed size, [search_cover_ants] looks for a random cover using
##'     a version of the ant-colony optimization heuristic,
##'     [find_cover_BB] finds covers using a branch-and-bound
##'     technique, [plot_cover] plots a cover.
##' @references Korte, Vygen \emph{Combinatorial Optimization. Theory
##'     and Algorithms.}
##' @examples 
##' ## Example with known vertex cover
##' K25 <- make_full_graph(25)   # Cover of size 24
##' X0 <- build_cover_approx(K25)
##' X0$size  # 24
##' plot_cover(X0, K25)
##'
##' ## Vertex-cover of a random graph
##' set.seed(1)
##' n <- 25
##' g <- sample_gnp(n, p=0.25)
##' X2 <- build_cover_approx(g)
##' X2$size   # 20
##' plot_cover(X2, g)
##'
##' @encoding UTF-8
##' @md 
##' @export 
build_cover_approx <- function(G) {
    n <- gorder(G)
    q <- gsize(G)
    X <- rep(FALSE, n)
    h <- G
    eG <- as_edgelist(G)
    M <- rep(FALSE, q)
    for (i in 1:q) {
	e <- eG[i,]
	X[e] <- TRUE
	M[i] <- TRUE
	h[e,] <- FALSE
	if (gsize(h) == 0) break
    }
    list(set = as.numeric(V(G)[X]), size = sum(X))    
}

##' Local search to improve a cover by using "neighboring" vertex
##'     subsets differing in just one element from the initial subset.
##' 
##' Given some cover specified by a vertex subset X in a graph, this
##'     routine scans the neighboring subsets obtained from X by
##'     removing a vertex from X looking for a smaller cover.  If such
##'     a cover is found, it replaces the starting cover and the
##'     search starts again.  This iterative procedure continues until
##'     no smaller cover can be found.  Of course, the resulting cover
##'     is only a local minimum.
##'
##' @title Improving a cover with local search
##' @param G A graph
##' @param X A cover list with components $set, $size as returned by
##'     routines [build_cover_greedy] or [build_cover_approx].  X
##'     represents the cover to be improved
##' @return A list with two components: $set contains the subset of
##'     V(g) representing the cover and $size contains the number of
##'     vertices of the cover.
##' @author Cesar Asensio
##' @seealso [is_cover] checks if a vertex subset is a vertex cover,
##'     [build_cover_greedy] builds a cover using a greedy heuristic,
##'     [build_cover_approx] builds a cover using a 2-approximation
##'     algorithm, [search_cover_random] looks for a random cover of
##'     fixed size, [search_cover_ants] looks for a random cover using
##'     a version of the ant-colony optimization heuristic,
##'     [find_cover_BB] finds covers using a branch-and-bound
##'     technique, [plot_cover] plots a cover.
##' @examples 
##' set.seed(1)
##' n <- 25
##' g <- sample_gnp(n, p=0.25)  # Random graph
##'
##' X1 <- build_cover_greedy(g)
##' X1$size    # 17
##' plot_cover(X1, g)
##'
##' X2 <- build_cover_approx(g)
##' X2$size    # 20
##' plot_cover(X2, g)
##'
##' X3 <- improve_cover_flip(g, X1)
##' X3$size    # 17 : Not improved
##' plot_cover(X3,g)
##'
##' X4 <- improve_cover_flip(g, X2)
##' X4$size    # 19 : It is improved by a single vertex
##' plot_cover(X4,g)
##'
##' # Vertex subsets or n-1 elements are always vertex covers:
##' for (i in 1:25) {
##'    X3 <- improve_cover_flip(g, list(set = setdiff(1:25,i), size = 24))
##'    print(X3$size)
##' } # 19 18 18 18 18 18 17 20 19 17 17 18 18 18 17 19 20 19 19 17 19 19 19 19 19
##'  
##' @encoding UTF-8
##' @md 
##' @export 
improve_cover_flip <- function(G, X) {
    n <- gorder(G)
    vG <- 1:n
    eG <- as_edgelist(G)
    Xmin <- Xb <- Xv <- vG %in% X$set
    sX <- sXmin <- sXmin.old <- X$size
    while(TRUE) {
	for (i in vG) {
	    Xv <- Xb
	    if (Xv[i]) {
		Xv[i] <- FALSE
	    } else {
		next
	    }
	    if (is_cover(vG[Xv], eG)) {
		Xmin <- Xv
		sXmin <- sum(Xmin)
	    }
	}
	if (sXmin < sXmin.old) {
	    sX <- sXmin.old <- sXmin
	    Xb <- Xmin
	} else {
	    break
	}
    }
    list(set = vG[Xmin], size = sXmin)
}

##' Check if some vertex subset of a graph covers all its edges.
##'
##' The routine simply goes through the edge list of the graph to see
##'      if both ends of each edge are inside the vertex subset to be
##'      checked.  When an edge with both ends ouside X is
##'      encountered, the routine returns FALSE; otherwise, it returns
##'      TRUE.
##' 
##' @title Check vertex cover
##' @param X Vertex subset to check.
##' @param eG Edgelist of the graph as returned by [as_edgelist]
##' @return Boolean: TRUE if X is a vertex cover of the graph
##'     represented by eG, FALSE otherwise.
##' @author Cesar Asensio
##' @seealso [build_cover_greedy] builds a cover using a greedy heuristic,
##'     [build_cover_approx] builds a cover using a 2-approximation
##'     algorithm, [improve_cover_flip] improves a cover using local
##'     search, [search_cover_random] looks for a random cover of
##'     fixed size, [search_cover_ants] looks for a random cover using
##'     a version of the ant-colony optimization heuristic,
##'     [find_cover_BB] finds covers using a branch-and-bound
##'     technique, [plot_cover] plots a cover.
##' @examples
##' set.seed(1)
##' n <- 25
##' g <- sample_gnp(n, p=0.25)  # Random graph
##' eg <- as_edgelist(g)
##'
##' X1 <- build_cover_greedy(g)
##' is_cover(X1$set, eg)   # TRUE
##' is_cover(c(1:10),eg)   # FALSE
##' plot_cover(list(set = 1:10, size = 10), g)  # See uncovered edges
##' 
##' @encoding UTF-8
##' @md 
##' @export 
is_cover <- function(X, eG) {
    q <- nrow(eG)
    verdict <- TRUE
    for (i in 1:q) {
	e <- eG[i,]
	if (sum(e %in% X) == 0) {
	    verdict <- FALSE
	    break
	}
    }
    verdict
}


##' Random algorithm for vertex-cover.
##'
##' It builds N random vertex sets by inserting elements with
##'     probability p, and it verifies if the subset so chosen is a
##'     vertex cover by running [is_cover] on it.  It is very
##'     difficult to find a good vertex cover in this way, so this
##'     algorithm is very inefficient and it finds no specially good
##'     covers.
##'
##' Currently, this function is *not* exported.  The random sampling
##'     performed by [search_cover_random] is faster and more
##'     efficient.
##' 
##' @title Random vertex covers
##' @param G Graph.
##' @param N Number of random vertex set to try.
##' @param p Probability of each element to be selected.
##' @return A list with four components: $set contains the subset of
##'     V(g) representing the cover and $size contains the number of
##'     vertices of the cover; $found is the number of vertex covers
##'     found and $failed is the number of generated subset that were
##'     not vertex covers. 
##' @author Cesar Asensio
##' @examples
##' n <- 25
##' g <- sample_gnp(n, p=0.25)  # Random graph
##' X5 <- build_cover_random(g,10000,p=0.65)
##' X5$size            # 19
##' plot_cover(X5, g)
##' X6 <- improve_cover_flip(g, X5)   # Improved : 17
##' plot_cover(X6, g)
##' 
##' @encoding UTF-8
##' @md 
##' @export 
build_cover_random <- function(G, N, p = 0.75) {
    n <- gorder(G)
    vG <- 1:n
    eG <- as_edgelist(G)
    Xmin <- X <- rep(TRUE, n)
    sXmin <- n
    sX <- num <- numF <- 0
    h <- rep(0,n)
    for (i in 1:N) {
	X <- sample(c(TRUE,FALSE), n, replace = TRUE, prob = c(p,1-p))
	if (is_cover(vG[X], eG)) {
	    num <- num + 1
	    sX <- sum(X)
	    h[sX] <- h[sX] + 1
	    if (sX < sXmin) {
		Xmin <- X
		sXmin <- sX
	    }
	} else { numF <- numF + 1}
    }
    list(set = vG[Xmin], size = sXmin, found = num, failed = numF,
	 hist = h)
}


##' Random algorithm for vertex-cover.
##'
##' This routine performs N iterations of the simple procedure of
##'     selecting a random sample of size k on the vertex set of a
##'     graph and check if it is a vertex cover, counting successes
##'     and failures.  The last cover found is returned, or an empty
##'     set if none is found.
##' 
##' @title Random vertex covers
##' @param G Graph.
##' @param N Number of random vertex set to try.
##' @param k Cardinality of the random vertex sets generated by the
##'     algorithm.
##' @param alpha Exponent of the probability distribution from which
##'     vertices are drawn: P(v) ~ d(v)^alpha.
##' @return A list with four components: $set contains the subset of
##'     V(g) representing the cover and $size contains the number of
##'     vertices of the cover (it coincides with k).  $found is the
##'     number of vertex covers found and $failed is the number of
##'     generated subset that were not vertex covers.
##' @author Cesar Asensio
##' @seealso [is_cover] checks if a vertex subset is a vertex cover,
##'     [build_cover_greedy] builds a cover using a greedy heuristic,
##'     [build_cover_approx] builds a cover using a 2-approximation
##'     algorithm, [improve_cover_flip] improves a cover using local
##'     search, [search_cover_ants] looks for a random cover using a
##'     version of the ant-colony optimization heuristic,
##'     [find_cover_BB] finds covers using a branch-and-bound
##'     technique, [plot_cover] plots a cover.
##' @examples 
##' set.seed(1)
##' n <- 25
##' g <- sample_gnp(n, p=0.25)  # Random graph
##' X7 <- search_cover_random(g, 10000, 17, alpha = 3)
##' plot_cover(X7, g)
##' X7$found    # 21 (of 10000) covers of size 17
##'
##' ## Looking for a cover of size 16...
##' X8 <- search_cover_random(g, 10000, 16, alpha = 3)  # ...we don't find any!
##' plot_cover(X8, g)   # All edges uncovered
##' X8$found    # 0
##' 
##' @encoding UTF-8
##' @md 
##' @export 
search_cover_random <- function(G, N, k, alpha = 1) {
    n <- gorder(G)
    vG <- 1:n
    eG <- as_edgelist(G)
    p <- degree(G)^alpha
    Xmin <- X <- c()
    num <- 0
    for (i in 1:N) {
	X <- sample(vG, k, prob = p)
	if (is_cover(X, eG)) {
	    Xmin <- X
	    num <- num + 1
	}
    }
    list(set = Xmin, size = length(Xmin), found = num, failed = N - num)
}

##' Ant colony optimization (ACO) heuristic algorithm to search for a
##'     vertex cover of small size in a graph.  ACO is a random
##'     algorithm; as such, it yields better results when longer
##'     searches are run.  To guess the adequate parameter values
##'     resulting in better performance in particular instances
##'     requires some experimentation, since no universal values of
##'     the parameters seem to be appropriate to all examples.
##'
##' ACO is an optimization paradigm that tries to replicate the
##'     behavior of a colony of ants when looking for food.  Ants
##'     leave after them a soft pheromone trail to help others follow
##'     the path just in case some food has been found.  Pheromones
##'     evaporate, but following again the trail reinforces it, making
##'     it easier to find and follow.  Thus, a team of ants search a
##'     vertex cover in a graph, leaving a pheromone trail on the
##'     chosen vertices.  At each step, each ant decides the next
##'     vertex to add based on the pheromone level and on the degree
##'     of the remaining vertices, according to the formula \eqn{P(v)
##'     ~ phi(v)^alpha*exp(beta*d(v))}, where \eqn{phi(v)} is the
##'     pheromone level, \eqn{d(v)} is the degree of the vertex
##'     \eqn{v}, and \eqn{alpha}, \eqn{beta} are two exponents to
##'     broad or sharpen the probability distribution.  After each
##'     vertex has been added to the subset, its incident adges are
##'     removed, following a randomized version of the greedy
##'     heuristic.  In a single iteration, each ant builds a vertex
##'     cover, and the best of them is recorded.  Then the pheromone
##'     level of the vertices of the best cover are enhanced, and the
##'     remaining pheromones begin to evaporate.
##' 
##' Default parameter values have been chosen in order to find the
##'     optimum in the examples considered below.  However, it cannot
##'     be guarateed that this is the best choice for all cases.  Keep
##'     in mind that no polynomial time exact algorithm can exist for
##'     the VCP, and thus harder instances will require to fine-tune
##'     the parameters.  In any case, no guarantee of optimality of
##'     covers found by this method can be given, so they might be
##'     improved further by other methods.
##' 
##' @title Ant colony optimization algorithm for Vertex-Cover
##' @param g Graph.
##' @param K Number of ants per iteration.
##' @param N Number of iterations.
##' @param alpha Exponent of the pheronome index, see details.
##' @param beta Exponent of the vertex degree, see details.
##' @param dt Pheromone increment.
##' @param rho Pheromone evaporation rate.
##' @param verb Boolean; if TRUE (default) it echoes to the console
##'     the routine progress .
##' @return A list with three components: $set contains the subset of
##'     V(g) representing the cover and $size contains the number of
##'     vertices of the cover; $found is the number of vertex covers
##'     found in subsequent iterations (often they are repeated, that
##'     is, different explorations may find the same vertex cover).
##' @author Cesar Asensio
##' @seealso [is_cover] checks if a vertex subset is a vertex cover,
##'     [build_cover_greedy] builds a cover using a greedy heuristic,
##'     [build_cover_approx] builds a cover using a 2-approximation
##'     algorithm, [improve_cover_flip] improves a cover using local
##'     search, [search_cover_random] looks for a random cover of
##'     fixed size, [find_cover_BB] finds covers using a
##'     branch-and-bound technique, [plot_cover] plots covers.
##' @examples
##' set.seed(1)
##' g <- sample_gnp(25, p=0.25)                  # Random graph
##'
##' X6 <- search_cover_ants(g, K = 20, N = 10)
##' plot_cover(X6, g)
##' X6$found
##' 
##' @encoding UTF-8
##' @md 
##' @export 
search_cover_ants <- function(g, K, N, alpha = 2, beta = 2, dt = 1,
			      rho = 0.1, verb=TRUE) {
    n <- gorder(g)
    Q <- gsize(g)
    V <- as.numeric(V(g))
    dg <- degree(g)
    pher <- rep(1, n)
    Co.min <- matrix(0, nrow = N, ncol = n)
    Nu.min <- rep(0, N)
    zeroM <- matrix(0, nrow = K, ncol = n)
    zeroV <- rep(0, K)
    if (verb) {
	cat("\nAnt-Colony algorithm for VERTEX-COVER\n\n   ",
	    N, " iterations with ", K, " ants ",
	    "(alpha = ", alpha, " / beta = ", beta,
	    " / rho = ", rho,")\n",
	    "   Graph instance: \"", g$name,"\" (n = ", n,
	    ", q = ", Q, ")\n",
	    "   Iteration: ", sep="")
    }
    for (i in 1:N) {       # Iteration of explorations
	if (verb) { cat(i, ",", sep="") }
	covers <- zeroM
	numins <- zeroV
	for (j in 1:K) {   # Send K explorer ants
	    h <- g
	    q <- Q
	    cov <- rep(0,n)
	    while (q>0) {  # Randomly build cover
		dh <- degree(h)[cov == 0]
		p <- exp(beta*dh)*pher[cov == 0]^alpha
		vp <- sample(V[cov == 0], size = 1, prob = p)
		cov[vp] <- 1
		h[vp,] <- FALSE
		q <- gsize(h)
	    }
	    covers[j,] <- cov       # Keep cover...
	    numins[j] <- sum(cov)   # ...its size...
	}
	min <- which.min(numins)    # ...and select the smaller one
	Nu.min[i] <-  numins[min]
	Co.min[i,] <- covers[min,]
	pher <- (1-rho)*pher             # Pheromone evaporation
	for (k in 1:n) {                 # Pheromone enhancement
	    pher[covers[min,] == 1] <- pher[covers[min,] == 1] + rho*dt
	}
    }
    nu.min <- min(Nu.min)
    co.min <- c()
    for (k in 1:nrow(Co.min)) {
	if (Nu.min[k] == nu.min) {
	    co.min <- rbind(co.min, V[Co.min[k,] == 1])
	}
    }
    found <- nrow(co.min)
    if (verb) {
	cat("DONE\n   Found cover with ", nu.min,
	    " vertices:\n", "    (", sep="")
	cat(co.min[1,], sep=",")
	cat(") [1 of ", found, "]\n\n",sep="")
    }
    list(set = co.min[1,], size = nu.min, found = found)
}

##' This routine performs a version of the Branch-and-Bound algorithm
##'     for the VCP.  It is an exact algorithm with exponential
##'     worst-case running time; therefore, it can be run only with a
##'     small number of vertices.
##'
##' The algorithm traverses a binary tree in which each bifurcation
##'     represents if a vertex is included in or excluded from a
##'     partial cover.  The leaves of the tree represent vertex
##'     subsets; the algorithm checks if at some point the partial
##'     cover cannot become a full cover because of too many uncovered
##'     edges with too few remaining vertices to decide.  In this way,
##'     the exponential complexity is somewhat reduced.  Furthermore,
##'     the vertices are considered in decreasing degree order, as in
##'     the greedy algorithm, so that some cover is found in the early
##'     steps of the algorithm and thus a good upper bound on the
##'     solution can be used to exclude more subsets from being
##'     explored.  The full algorithm has been extracted from the
##'     reference below.
##'
##' In this routine, the binary tree search is implemented by
##'     recursive calls (that is, a dynamic programming algorithm).
##'     Although the worst case time complexity is exponential (recall
##'     that the Minimum Vertex Cover Problem is NP-hard), the
##'     approach is fairly efficient for a branch-and-bound technique.
##'
##' The tree node in which the algorithm is when it is called (by the
##'     user or by itself) is encoded in a sequence of vertex marks.
##'     Marks come in three flavors: "C" is assigned to "Covered"
##'     vertices, that is, already included in the partial cover.  "U"
##'     is asigned to "Uncovered" vertices, that is, those excluded
##'     from the partial cover.  Mark "F" is assigned to "Free"
##'     vertices, those not considered yet by the algorithm; one of
##'     them is considered in the actual function call, and the
##'     subtree under this vertex is explored before returning.  This
##'     mark sequence starts and ends with all vertices marked "F",
##'     and is used only by the algorithm, which modifies and passes
##'     it on to succesive calls to itself.
##'
##' When the verb argument is TRUE, the routine echoes to the console
##'     the newly found cover only if it is better than the last.
##'     This report includes the size, actual cover and number call of
##'     the routine.
##' 
##' The routine can drop the best cover found so far in a file so that
##'     the user can stop the run afterwards; this technique might be
##'     useful when the full run takes too much time to complete.
##' 
##' @title Branch-and-Bound algorithm for the Vertex-Cover problem
##' @param g Graph.
##' @param verb Boolean: Should echo each newly found cover to the
##'     console? Defaults to TRUE.
##' @param save.best.result Boolean: Should the algorithm save the
##'     result of the algorithm in a file?  It defaults to FALSE.
##'     When save.best.result = TRUE, a file is created with the
##'     variable "Xbest" being the best result achieved by the
##'     algorithm before its termination. 
##' @param filename.best.result Name of the file created when
##'     save.best.result = TRUE.  It defaults to
##'     "best_result_find_cover_BB.Rdata".
##' @param nu Size of the best cover currently found.
##' @param X Partial cover.
##' @param Xmin Best cover found so far.
##' @param marks Mark sequence storing the current state of the
##'     algorithm, see details.
##' @param call Number of recursive calls performed by the algorithm.
##' @return A list with five components: $set contains the subset of
##'     V(g) representing the cover and $size contains the number of
##'     vertices of the cover. Component $call is the number of calls
##'     the algorithm did on itself.  The remaining components are
##'     used to transfer the state of the algorithm in the search
##'     three from one call to the next; they are $partial, the
##'     partially constructed cover, and $marks, a sequence encoding
##'     the tree node in which the algorithm is when this function is
##'     called, see details.
##' @author Cesar Asensio
##' @seealso [is_cover] checks if a vertex subset is a vertex cover,
##'     [build_cover_greedy] builds a cover using a greedy heuristic,
##'     [build_cover_approx] builds a cover using a 2-approximation
##'     algorithm, [improve_cover_flip] improves a cover using local
##'     search, [search_cover_random] looks for a random cover of
##'     fixed size, [search_cover_ants] looks for a random cover using
##'     a version of the ant-colony optimization heuristic,
##'     [plot_cover] plots a cover.
##' @examples
##' set.seed(1)
##' g <- sample_gnp(25, p=0.25)        # Random graph
##' X7 <- find_cover_BB(g)
##' X7$size                            # Exact result: 16
##' X7$call                            # 108 recursive calls!
##' plot_cover(X7, g)
##'
##' ## Saving best result in a file (useful if the algorithm takes too
##' ## long and should be interrupted by the user)
##' ## It uses tempdir() to store the file
##' ## The variable is called "Xbest"
##' find_cover_BB(g, save.best.result = TRUE, filename.best.result = "BestResult_BB.Rdata")
##' 
##' @encoding UTF-8
##' @md 
##' @export 
find_cover_BB <- function(g, verb = TRUE,
			  save.best.result = FALSE,
			  filename.best.result = "best_result_find_cover_BB.Rdata",
			  nu    = gorder(g),
			  X     = c(),
			  Xmin  = c(),
			  marks = rep("F", gorder(g)),
			  call  = 0) {

    if (call == 0) {
	cat("Running branch-and-bound on a VCP with n =",
	    gorder(g), "vertices...")
    }

    call <- call + 1 # Call count
    J <- call

    ## Are there some remaining edges in the graph?  Otherwise,
    ## the graph is covered!
    q <- gsize(g)
    if (q == 0) {
	if (length(X) < nu) {
	    nu <- length(X)
	    Xmin <- X
	    if(verb) {
		cat("\nCover found:\n",
		    "  size = ", nu, "\n",
		    " cover = ", paste(X, collapse = ","), "\n",
		    "  call = ",  J, "\n",
		    sep="")
	    }
	    if (save.best.result)  {
		assign("Xbest", value = list(set = Xmin, size = nu))
		dirsep <- if(Sys.info()["sysname"] == "Windows") { "\\" } else { "/" }
		tfile <- paste(tempdir(), filename.best.result, sep = dirsep)
		save(Xbest, file = tfile)
		cat("\nSaved best result found in variable \"Xbest\" in file ", tfile, "\n", sep = "")
	    }
	}
	## if (plot) { plot_cover(list(set = X, size = nu), gini) }
	return(list(size = nu, set = Xmin, call = J, partial = X, marks = marks))
    }

    ## We compute bounds:
    ## F = Number of vertices to reach "nu" in X
    ## D = Upper bound of the number of uncovered edges
    F <- nu - length(X)
    ds <- degree(g)
    D <- sum_g(ds[marks == "F"], F)

    ## If uncovered edges outnumber the bound, covering in not
    ## possible anymore
    if (D < q) {
	return(list(size = nu, set = Xmin, call = J, partial = X, marks = marks)) 
    }

    ## We choose the "F" (free) vertex with maximum degree, as in the
    ## greedy heuristic
    V <- as.numeric(V(g))
    k <- which.max(ds[marks == "F"])
    v0 <- V[marks == "F"][k]

    ## We mark the chosen vertex as "C" (covered)  [subtree "v0 covered"]
    marks[v0] <- "C"
    X <- c(X,v0)   # ...end we add it to the subset X

    ## We remove the edges incident to v0
    ## Before, we record its neighbors to restore them later
    Nv0 <- neighbors(g, v0)
    g[v0,] <- FALSE

    ## Entering subtree "v0 covered" with updated "g" and "X"
    s <- find_cover_BB(g, verb, save.best.result, filename.best.result, nu, X, Xmin, marks, J)
    J <- s$call
    nu <- s$size
    X <- s$partial
    Xmin <- s$set
    marks <- s$marks

    ## At this point, subtree "v0 covered" has been explored; we
    ## restore the previously removed edges, and we remove v0 from X
    for (v in Nv0) { g <- g + edge(v0,v) }
    X <- X[X != v0]

    ## We check the second bound: Is F > d(v0)?  If it is, we mark
    ## with "C" all neighbors; if it is not, marking would lead to a
    ## suboptimal cover, so we just free v0.
    if (F > length(Nv0)) { # Subtree "v0 uncovered"
	marks[v0] <- "U"
	aris <- c()
	for (v in Nv0) {
	    marks[v] <- "C"
	    X <- c(X, v)
	    Nv <- neighbors(g, v)
	    for (u in Nv) { # We take note of edges to erase...
		aris <- c(aris, v, u)
	    }
	    g[v,] <- FALSE             # ...and we remove them
	}

	## Entering subtree "v0 uncovered"
	s <- find_cover_BB(g, verb, save.best.result, filename.best.result, nu, X, Xmin, marks, J)
	J <- s$call
	nu <- s$size
	X <- s$partial
	Xmin <- s$set
	marks <- s$marks

	## We free neighbors and restore the previously removed edges
	for (v in Nv0) {
	    X <- X[X!=v]
	    marks[v] <- "F"
	}
	g <- g + edges(aris)
    }

    ## Once checked both subtrees, we free v0 and we are done!
    marks[v0] <- "F"
    if (sum(marks == "F") == gorder(g)) cat("done\n")
    return(list(size = nu, set = Xmin, call = J, partial = X, marks = marks))
}

##' Sum of the higher terms of a list
##'
##' It sorts the list L in decreasing order and returns the sum of the
##'     first m components of the ordered list.
##' 
##' @title Sum of the higher terms of a list
##' @param L A numeric sequence
##' @param m A scalar
##' @return The sum of the m higher terms of the list L
##' @author Cesar Asensio
##' @examples 
##' sum_g(1:10, 3)  # 8 + 9 + 10 = 27
##' 
##' @encoding UTF-8
##' @md 
##' @export 
sum_g <- function(L, m) { sum(sort(L, decreasing = TRUE)[1:m]) }

##' Random cut generation on a graph.  This function generates a
##'     hopefully large cut on a graph by randomly selecting vertices;
##'     it does not attempt to maximize the cut size or weigth, so it
##'     is intended to be used as part of some smarter strategy.
##'
##' It selects a random subset of the vertex set of the graph,
##'     computing the associated cut, its size and its weigth,
##'     provided by the user as a weight matrix.  If the weight
##'     argument w is NA, the weights are taken as 1.
##' 
##' @title Random cut generation on a graph
##' @param G Graph
##' @param w Weight matrix (defaults to NA).  It should be zero for
##'     those edges not in G
##' @return A list with four components: $set contains the subset of
##'     V(g) representing the cut, $size contains the number of edges
##'     of the cut, $weight contains the weight of the cut (which
##'     coincides with $size if w is NA) and $cut contains the edges
##'     of the cut, joining vertices inside $set with vertices outside
##'     $set.
##' @importFrom stats runif
##' @author Cesar Asensio
##' @seealso [build_cut_greedy] builds a cut using a greedy algorithm,
##'     [compute_cut_weight] computes cut size, weight and edges,
##'     [improve_cut_flip] uses local search to improve a cut obtained
##'     by other methods, [plot_cut] plots a cut.
##' @examples 
##' ## Example with known maximum cut
##' K10 <- make_full_graph(10)   # Max cut of size 25
##' c0 <- build_cut_random(K10)
##' c0$size  # Different results: 24, 21, ...
##' plot_cut(c0, K10)
##'
##' ## Max-cut of a random graph
##' set.seed(1)
##' n <- 25
##' g <- sample_gnp(n, p=0.25)
##' c1 <- build_cut_random(g)   # Repeat as you like
##' c1$size  # Different results: 43, 34, 39, 46, 44, 48...
##' plot_cut(c1, g)
##'
##' @encoding UTF-8
##' @md 
##' @export 
build_cut_random <- function(G, w=NA) {
    n <- gorder(G)
    S <- as.numeric(V(G)[runif(n) < 0.5])
    eG <- as_edgelist(G)
    compute_cut_weight(S, n, eG, w, return.cut = TRUE)
}


##' Plot of a cut in a graph.
##'
##' It plots a graph, then superimposes a cut, drawing the associated
##'     vertex set in a different color.
##' 
##' @title Cut plotting
##' @param K Cut to be plotted; an output list returned by some
##'     cut-building function, see below.
##' @param G Graph on which to superimpose the cut.
##' @return This function is called for its side effect of plotting.
##' @importFrom graphics legend
##' @importFrom graphics text
##' @author Cesar Asensio
##' @seealso [build_cut_random] builds a random cut,
##'     [build_cut_greedy] builds a cut using a greedy algorithm,
##'     [improve_cut_flip] uses local search to improve a cut obtained
##'     by other methods, [compute_cut_weight] computes cut size,
##'     weight and edges.
##' @examples 
##' K10 <- make_full_graph(10)   # Max cut of size 25
##' c0 <- build_cut_random(K10)
##' plot_cut(c0, K10)
##'
##' @encoding UTF-8
##' @md 
##' @export 
plot_cut <- function(K, G) {
    n <- gorder(G)
    z <- layout_with_gem(G)
    plot(G, layout=z)
    gcut <- make_graph(t(K$cut), dir=FALSE)
    vcol <- rep("pink", n)
    lcol <- rep("black", n)
    vcol[K$set] <- "red4"
    lcol[K$set] <- "white"
    plot(gcut, layout=z, add=TRUE, edge.color="red3", edge.width=2,
	 vertex.color=vcol, vertex.label.color=lcol)
    legend(-1.2, 1.2, legend = c("Set", "Cut"),
	   lwd = c(NA,2), pch = c(21, NA), col = c("black", "red3"),
	   pt.bg = c("red4", NA), pt.cex = c(2.5, NA))
    text(0, 1.2,
	 labels = paste0("Size = ", K$size, ", Weight = ", K$weight))
}

##' Compute cut weight and size from its associated vertex set.  It
##'     can also return the edges in the cut.
##'
##' In a graph, a cut K is defined by means of a vertex subset S as
##'     the edges joining vertices inside S with vertices outside S.
##'     This routine computes these edges and their associated weight.
##' 
##' @title Compute cut weight and size
##' @param S Subset of the vertex set of the graph.
##' @param n Size of the graph.
##' @param eG Edgelist of the graph as returned by [as_edgelist], that
##'     is, a matrix with q rows and 2 columns.  Note that this is the
##'     graph format used by the routine.
##' @param w Weight matrix or NA if all edge weights are 1.  It should
##'     be zero for those edges not in G
##' @param return.cut Boolean.  Should the routine return the edges in
##'     the cut? It defaults to FALSE.  When TRUE, the routine also
##'     returns the input subset S, for easier cut plotting with
##'     [plot_cut].
##' @return A list with two components: $size is the number of edges
##'     in the cut, $weight is the weight of the cut, that is, the sum
##'     of the weights of the edges in the cut.  If w=NA these two
##'     numbers coincide.  When return.cut is TRUE, there are two
##'     additional components of the list: $cut, which contains the
##'     edges in the cut as rows of a two-column matrix, and $set,
##'     which contains the input set, as a convenience for plotting
##'     with [plot_cut].
##' @author Cesar Asensio
##' @seealso [build_cut_random] builds a random cut,
##'     [build_cut_greedy] builds a cut using a greedy algorithm,
##'     [improve_cut_flip] uses local search to improve a cut obtained
##'     by other methods, [plot_cut] plots a cut.
##' @examples 
##' K10 <- make_full_graph(10)
##' S <- c(1,4,7)
##' compute_cut_weight(S, gorder(K10), as_edgelist(K10))
##' cS <- compute_cut_weight(S, gorder(K10), as_edgelist(K10),
##'     return.cut = TRUE)
##' plot_cut(cS, K10)
##' 
##' @encoding UTF-8
##' @md 
##' @export 
compute_cut_weight <- function(S, n, eG, w=NA, return.cut = FALSE) {
    M <- 2*as.numeric(1:n %in% S) - 1
    q <- nrow(eG)
    k <- rep(FALSE, q)
    val <- 0
    for (e in 1:q) {
	u <- eG[e, 1]
	v <- eG[e, 2]
	if (prod(M[c(u,v)]) == -1) {
	    k[e] <- TRUE
	    if (is.na(w[1])) {
		val <- val + 1
	    } else {
		val <- val + w[u, v]
	    }
	}
    }
    if (return.cut) {
	list(set = S, size = sum(k), weight = val, cut = eG[k,])
    } else {
	list(size = sum(k), weight = val)
    }
}

##' This routine uses a greedy algorithm to build a cut with large
##'     weight.  This is a 2-approximation algorithm, which means that
##'     the weight of the cut returned by this algorithm is larger
##'     than half the maximum possible cut weight for a given graph.
##'
##' The algorithm builds a vertex subset S a step a a time.  It starts
##'     with S = c(v1), and with vertices v1 and v2 marked.  Then it
##'     iterates from vertex v3 to vn checking if the weight of the
##'     edges joining vi with marked vertices belonging to S is less
##'     than the weight of the edges joining vi with marked vertices
##'     not belonging to S.  If the former weight is less than the
##'     latter, then vi is adjoined to S.  At the end of each
##'     iteration, vertex vi is marked.  When all vertices are marked
##'     the algorithm ends and S is already built.
##'
##' @title Greedy algorithm aimed to build a large weight cut in a
##'     graph
##' @param G Graph
##' @param w Weight matrix (defaults to NA).  It should be zero for
##'     those edges not in G
##' @return A list with four components: $set contains the subset of
##'     V(g) representing the cut, $size contains the number of edges
##'     of the cut, $weight contains the weight of the cut (which
##'     coincides with $size if w is NA) and $cut contains the edges
##'     of the cut, joining vertices inside $set with vertices outside
##'     $set.
##' @author Cesar Asensio
##' @references Korte, Vygen \emph{Combinatorial Optimization. Theory
##'     and Algorithms.}
##' @seealso [build_cut_random] builds a random cut,
##'     [improve_cut_flip] uses local search to improve a cut obtained
##'     by other methods, [compute_cut_weight] computes cut size,
##'     weight and edges, [plot_cut] plots a cut.
##' @examples 
##' ## Example with known maximum cut
##' K10 <- make_full_graph(10)   # Max cut of size 25
##' c0 <- build_cut_greedy(K10)
##' c0$size  # 25
##' plot_cut(c0, K10)
##'
##' ## Max-cut of a random graph
##' set.seed(1)
##' n <- 25
##' g <- sample_gnp(n, p=0.25)
##' c2 <- build_cut_greedy(g)
##' c2$size   # 59
##' plot_cut(c2, g)
##'
##' @encoding UTF-8
##' @md 
##' @export 
build_cut_greedy <- function(G, w=NA) {
    n <- gorder(G)
    S <- c(1)
    R <- c(1,2)
    if (is.na(w[1])) {
	w <- as_adjacency_matrix(G)
    }
    for (i in 3:n) {
	RiS <- intersect(R,S)
	RsS <- setdiff(R,S)
	if (sum(w[i, RiS]) < sum(w[i, RsS])) {
	    S <- c(S, i)
	}
	R <- c(R, i)
    }
    eG <- as_edgelist(G)
    compute_cut_weight(S, n, eG, w, return.cut = TRUE)
}

##' Local search to improve a cut by using "neighboring" vertex
##'     subsets differing in just one element from the initial subset.
##' 
##' Given some cut specified by a vertex subset S in a graph, this
##'     routine scans the neighboring subsets obtained from S by
##'     adding/removing a vertex from S looking for a larger cut.  If
##'     such a cut is found, it replaces the starting cut and the
##'     search starts again.  This iterative procedure continues until
##'     no larger cut can be found.  Of course, the resulting cut is
##'     only a local maximum.
##'
##' @title Improving a cut with local search
##' @param G A graph
##' @param K A cut list with components $set, $size, $weight and $cut
##'     as returned by routines [build_cut_greedy], [build_cut_random]
##'     or [compute_cut_weight].  Only the $set and $weight components
##'     are used.  K represents the cut to be improved
##' @param w Weight matrix (defaults to NA).  It should be zero for
##'     those edges not in G
##' @param return.cut Boolean.  Should the routine return the cut?  It
##'     is passed on to [compute_cut_weight] on return.  It defaults
##'     to TRUE
##' @return A list with four components: $set contains the subset of
##'     V(g) representing the cut, $size contains the number of edges
##'     of the cut, $weight contains the weight of the cut (which
##'     coincides with $size if w is NA) and $cut contains the edges
##'     of the cut, joining vertices inside $set with vertices outside
##'     $set.  When return.cut is FALSE, components $set and $cut are
##'     omitted.
##' @author Cesar Asensio
##' @seealso [build_cut_random] builds a random cut,
##'     [build_cut_greedy] builds a cut using a greedy algorithm,
##'     [compute_cut_weight] computes cut size, weight and edges,
##'     [plot_cut] plots a cut.
##' @examples 
##' set.seed(1)
##' n <- 25
##' g <- sample_gnp(n, p=0.25)  # Random graph
##'
##' c1 <- build_cut_random(g)
##' c1$size    # 44
##' plot_cut(c1, g)
##'
##' c2 <- build_cut_greedy(g)
##' c2$size    # 59
##' plot_cut(c2, g)
##'
##' c3 <- improve_cut_flip(g, c1)
##' c3$size    # 65
##' plot_cut(c3,g)
##'
##' c4 <- improve_cut_flip(g, c2)
##' c4$size    # 60
##' plot_cut(c4,g)
##'
##' @encoding UTF-8
##' @md 
##' @export 
improve_cut_flip <- function(G, K, w=NA, return.cut = TRUE) {
    n <- gorder(G)
    vG <- 1:n
    eG <- as_edgelist(G)
    Smax <- S <- Sv <- vG %in% K$set
    wK <- wKmax <- wKmax.old <- K$weight
    while(TRUE) {
	for (i in vG) {
	    Sv <- S
	    Sv[i] <- !Sv[i]
	    cSv <- compute_cut_weight(vG[Sv], n, eG, w)
	    wK <- cSv$weight
	    if (wK > wKmax) {
		Smax <- Sv
		wKmax <- wK
	    }
	}
	if (wKmax > wKmax.old) {
	    wKmax.old <- wKmax
	    S <- Smax
	} else {
	    break
	}
    }
    compute_cut_weight(vG[Smax], n, eG, w, return.cut = return.cut)
}

##' Crossover sequence operation for use in the genetic cut-search algorithm.
##'
##' This operation takes two sequences of the same lenght "n" and
##'     splits them in two at a crossover point between 1 and "n-1".
##'     Then it produces two "offsprings" by interchanging the pieces
##'     and gluing them together.
##'
##' The crossover point can be specified in argument cpoint.  By
##'     providing NA (the default), cpoint is chosen randomly.
##'
##' Note that this crossover operation is the "classic" crossover
##'     included in the original genetic algorithm, and it is adequate
##'     when applied to binary sequences.  However, when applied to
##'     permutations, the result of this function can have repeated
##'     elements; hence, it is not adequate for the TSP.
##' 
##' @title Crossover of sequences
##' @param s1 Sequence
##' @param s2 Sequence of the same lenght as s1
##' @param cpoint Crossover point, an integer between 1 and
##'     length(s1)-1.  Defaults to NA, in which case it will be
##'     randomly chosen
##' @return A two-row matrix.  Rows are the offsprings produced by the
##'     crossover
##' @author Cesar Asensio
##' @seealso [search_cut_genetic] genetic algorithm cut-search,
##'     [mutate_binary_sequence] binary sequence mutation
##' @examples
##' set.seed(1)
##' s1 <- sample(0:1, 10, replace = TRUE)
##' s2 <- sample(0:1, 10, replace = TRUE)
##' crossover_sequences(s1, s2)
##' 
##' set.seed(1)
##' s1 <- sample(1:10, 10)
##' s2 <- sample(1:10, 10)
##' crossover_sequences(s1, s2, cpoint = 5)
##' 
##' @encoding UTF-8
##' @md 
##' @export 
crossover_sequences <- function(s1, s2, cpoint = NA) {
    n <- length(s1)
    if (is.na(cpoint)) cpoint <- sample(1:(n-1),1)
    s1.1 <- s1[1:cpoint]
    s1.2 <- s1[(cpoint+1):n]
    s2.1 <- s2[1:cpoint]
    s2.2 <- s2[(cpoint+1):n]
    rbind(c(s1.1, s2.2), c(s2.1,s1.2))
}


##' Mutation of binary sequences for use in the genetic algorithm
##'
##' This routine takes a binary sequence and it flips ("mutates") each
##'     bit with a fixed probability.  In the genetic algorithm
##'     context, this operation randomly explores regions of
##'     configuration space which are far away from the starting
##'     point, thus trying to avoid local optima.  The fitting
##'     function values of mutated individuals are generically very
##'     poor, and this behavior is to be expected.  Thus, mutation is
##'     not an optimization procedure per se.
##' 
##' @title Binary sequence mutation
##' @param s Sequence consisting of 0 and 1
##' @param p Mutation probability.  Defaults to 0.1
##' @return A mutated binary sequence 
##' @author Cesar Asensio
##' @importFrom stats runif
##' @seealso [search_cut_genetic] genetic cut-searching algorithm,
##'     [crossover_sequences] crossover operation
##' @examples
##' set.seed(1)
##' s <- sample(0:1, 10, replace = TRUE)  # 0 0 1 1 0 1 1 1 1 0
##' mutate_binary_sequence(s, p = 0.5)    # 1 1 1 0 0 0 1 1 0 0
##' mutate_binary_sequence(s, p = 1)      # 1 1 0 0 1 0 0 0 0 1
##' 
##' @encoding UTF-8
##' @md 
##' @export 
mutate_binary_sequence <- function(s, p = 0.1) {
    n <- length(s)
    r <- runif(n)
    m <- s[r < p]
    s[r < p] <- 1 - m
    s
}

##' Genetic algorithm for Max-Cut.  In addition to crossover and
##'     mutation, which are described below, the algorithm performs
##'     also local search on offsprings and mutants.
##'
##' This algorithm manipulates cuts by means of its associated binary
##'     sequences defined as follows.  Each cut K is defined by its
##'     associated vertex subset S of V(G): K contains all edges
##'     joining vertices inside S with vertices outside S.  If
##'     |V(G)|=n, we can construct a n-bit binary sequence b =
##'     (b1,b2,...,bn) with bi = 1 if vertex vi belongs to S, and 0
##'     otherwise.
##' 
##' The genetic algorithm consists of starting with a cut population,
##'     where each cut is represented by its corresponding binary
##'     sequence defined above, and thus the population is simply a
##'     binary matrix.  This initial cut population can be provided by
##'     the user or can be random.  The initial population can be the
##'     output of a previous run of the genetic algorithm, thus
##'     allowing a chained execution.  Then the routine sequentially
##'     perform over the cuts of the population the
##'     \strong{crossover}, \strong{mutation}, \strong{local search}
##'     and \strong{selection} operations.
##'
##' The \strong{crossover} operation takes two cuts as "parents" and
##'     forms two "offsprings" by cutting and interchanging the binary
##'     sequences of the parents; see [crossover_sequences] for more
##'     information.
##'
##' The \strong{mutation} operation performs a "small" perturbation of
##'     each cut trying to escape from local optima.  It uses a random
##'     flip on each bit of the binary sequence associated with the
##'     cut, see [mutate_binary_sequence] for more information.
##'
##' The \strong{local search} operation takes the cuts found by the
##'     crossover and mutation operations and improves them using some
##'     local search heuristic, in this case [improve_cut_flip].  This
##'     allows this algorithm to approach local maxima faster.
##'
##' The \strong{selection} operation is used when selecting pairs of
##'     parents for crossover and when selecting individuals to form
##'     the population for the next generation.  In both cases, it
##'     uses a probability exponential in the weight with rate
##'     parameter "beta", favouring the better fitted to be selected.
##'     Lower values of beta favours the inclusion of cuts with worse
##'     fitting function values.  When selecting the next population,
##'     the selection uses \emph{elitism}, which is to save the best
##'     fitted individuals to the next generation; this is controlled
##'     with parameter "elite".
##'
##'
##' The usefulness of the crossover and mutation operations stems from
##'     its abitily to escape from the local maxima.  Of course, more
##'     iterations (Ngen) and larger populations (Npop) might improve
##'     the result, but recall that no random algorithm can guarantee
##'     to find the optimum of a given Max-Cut instance.
##'
##' This algorithm calls many times the routines [compute_cut_weight],
##'     [crossover_sequences], [mutate_binary_sequence] and
##'     [improve_cut_flip]; therefore, it is not especially efficient
##'     when called on large problems or with high populations or many
##'     generations.  Please consider chaining the algorithm:  perform
##'     short runs, using the output of a run as the input of the
##'     next.
##'
##' @title Genetic Algorithm for Max-Cut
##' @param G Graph.
##' @param w Weight matrix.
##' @param Npop Population size.
##' @param Ngen Number of generations (iterations of the algorithm).
##' @param pmut Mutation probability.  It defaults to 0.1.
##' @param beta Control parameter of the crossing and selection
##'     probabilities.  It defaults to 1.
##' @param elite Number of better fitted individuals to pass on to the
##'     next generation.  It defaults to 2.
##' @param Pini Initial population.  If it is NA, a random initial
##'     population of Npop individuals is generated.  Otherwise, it
##'     should be a matrix; each row should be an individual (a
##'     permutation of the 1:n sequence) and then Npop is set to the
##'     number of rows of Pini.  This option allows to chain several
##'     runs of the genetic algorithms, which could be needed in the
##'     hardest cases.
##' @param verb Boolean to activate console echo.  It defaults to
##'     TRUE.
##' @param log Boolean to activate the recording of the weights of all
##'     cuts found by the algorithm.  It defaults to FALSE.
##' @return A list with several components: $set contains the subset
##'     of V(g) representing the cut, $size contains the number of
##'     edges of the cut, $weight contains the weight of the cut
##'     (which coincides with $size if w is NA) and $cut contains the
##'     edges of the cut, joining vertices inside $set with vertices
##'     outside $set; $generation contains the generation when the
##'     maximum was found and $population contains the final cut
##'     population.  When log=TRUE, the output includes several lists
##'     of weights of cuts found by the algorithm, separated by
##'     initial cuts, offsprings, mutants, local maxima and selected
##'     cuts.
##' @references Hromkovic \emph{Algorithms for hard problems} (2004),
##'     Hartmann, Weigt, \emph{Phase transitions in combinatorial
##'     optimization problems} (2005).
##' @author Cesar Asensio
##' @seealso [crossover_sequences] performs crossover operation,
##'     [mutate_binary_sequence] performs mutation operation,
##'     [build_cut_random] builds a random cut, [build_cut_greedy]
##'     builds a cut using a greedy algorithm, [improve_cut_flip]
##'     improves a cut by local search, [compute_cut_weight] computes
##'     cut size, weight and edges, [plot_cut] plots a cut.
##' @examples 
##' set.seed(1)
##' n <- 10
##' g <- sample_gnp(n, p=0.5)  # Random graph
##' c5 <- search_cut_genetic(g)
##' plot_cut(c5, g)
##' improve_cut_flip(g, c5)     # It does not improve
##' for (i in 1:5) {            # Weights of final population
##'    s5 <- which(c5$population[i,] == 1)
##'    cs5 <- compute_cut_weight(s5, gorder(g), as_edgelist(g))
##'    print(cs5$weight)
##' }
##'
##' \donttest{
##' ## Longer examples
##' c5 <- search_cut_genetic(g, Npop=10, Ngen=50, log = TRUE)
##' boxplot(c5$Wini, c5$Woff, c5$Wmut, c5$Wvec, c5$Wsel,
##'         names=c("Ini", "Off", "Mut", "Neigh", "Sel"))
##'
##' set.seed(1)
##' n <- 20
##' g <- sample_gnp(n, p=0.25)
##' Wg <- matrix(sample(1:3, n^2, replace=TRUE), nrow=n)
##' Wg <- Wg + t(Wg)
##' A <- as_adjacency_matrix(g)
##' Wg <- Wg * A
##' c6 <- search_cut_genetic(g, Wg, Ngen = 9)   # Size 38, weigth 147
##' plot_cut(c6, g)
##' }
##' 
##' @encoding UTF-8
##' @md 
##' @export 
search_cut_genetic <- function(G, w = NA, Npop = 5, Ngen = 20,
			       pmut = 0.1, beta = 1, elite = 2,
			       Pini = NA, verb = TRUE, log = FALSE) {
    n <- gorder(G)
    vG <- 1:n
    eG <- as_edgelist(G)

    ## Initial population
    if (is.na(Pini)) {
	welc <- "Generated"
	Pini <- matrix(sample(0:1, n*Npop, replace = TRUE),
		       nrow = Npop, ncol = n)
    } else {
	welc <- "Received"
	Npop <- nrow(Pini)
    }

    ## Initialization
    Wmax <- 0
    P <- Pini
    U  <- U0 <- matrix(0, ncol = n, nrow = Npop + Npop + 2*Npop + 4*Npop)
    U[1:Npop,] <- P
    Wu <- Wu0 <- rep(0, 8*Npop)
    for (i in 1:Npop) {
	Wu[i] <- compute_cut_weight(vG[P[i,] == 1], n, eG, w)$weight
    }
    K <- Npop + 1

    if (log) {
	Wini <- Wu[1:Npop]
	Woff <- Wmut <- Wvec <- Wsel <- c()
    }
    if (verb) {
	cat(welc, " initial population (Npop = ", Npop, ")\n",
	    "Evolving Ngen = ", Ngen, " generations:\n", sep="")
    }

    ## Main loop
    for (i in 1:Ngen) {
	if (verb) cat("[Gen = ",i,"] : ", sep="")

	if (verb) cat("Crossover...")                  # Begin crossover 
	Noff <- 0
	Wpar <- Wu[1:Npop]
	pcross <- exp(beta*Wpar/mean(Wpar))
	for (m in 1:Npop) {
	    coup <- sample(1:Npop, 2, prob = pcross)
	    offs <- crossover_sequences(P[coup[1],], P[coup[2],])
	    for (j in 1:2) {
		repe <- FALSE
		off <- offs[j,]
		for (ki in 1:(K-1)) {
		    dis <- sum(abs(off - U[ki,]))
		    eq.test <- (dis == 0)
		    co.test <- (dis == n)
		    if (eq.test | co.test) {
			repe <- TRUE
			break
		    }
		}
		if (!repe) {
		    U[K,] <- off
		    Wo <- compute_cut_weight(vG[off == 1], n, eG, w)
		    Wu[K] <- Wo$weight
		    if (log) { Woff <- c(Woff, Wo$weight) }
		    Noff <- Noff + 1
		    K <- K + 1
		}
		if (Noff == Npop) { break }
	    }
	    if (Noff == Npop) { break }
	}
	if (verb) {
	    if (Noff == 1) { sfin <- "" } else { sfin <- "s" }
	    cat("(", Noff, " offspring", sfin, ") ",sep="")
	}                                               # End crossover

	if (verb) cat("Mutation...")                    # Begin mutation 
	Ncan <- Npop + Noff
	Nmut <- 0
	for (m in 1:Ncan) {
	    mut <- mutate_binary_sequence(U[m, ], pmut)
	    repe <- FALSE
	    for (k in 1:(K-1)) {
		dis <- sum(abs(mut - U[ki,]))
		eq.test <- (dis == 0)
		co.test <- (dis == n)
		if (eq.test | co.test) {
		    repe <- TRUE
		    break
		}
	    }
	    if (!repe) {
		U[K,] <- mut
		Wo <- compute_cut_weight(vG[mut == 1], n, eG, w)
		Wu[K] <- Wo$weight
		if (log) { Wmut <- c(Wmut, Wo$weight) }
		Nmut <- Nmut + 1
		K <- K + 1
	    }
	}
	if (verb) {
	    if (Nmut == 1) { sfin <- "" } else { sfin <- "s" }
	    cat("(", Nmut, " mutant", sfin, ") ",sep="")
	}                                                # End mutation

	if (verb) cat("Local search...")                 # Begin local search
	Ncan <- Npop + Noff + Nmut
	Nvec <- 0
	for (m in 1:Ncan) {
	    Ko <- list(set = vG[U[m, ] == 1], weight = Wu[m])
	    if (verb & !is.na(w[1])) cat("i")
	    vecW <- improve_cut_flip(G, Ko, w)
	    vec <- rep(0,n)
	    vec[vG %in% vecW$set] <- 1
	    repe <- FALSE
	    for (ki in 1:(K-1)) {
		if (verb & !is.na(w[1])) cat("-")
		dis <- sum(abs(vec - U[ki,]))
		eq.test <- (dis == 0)
		co.test <- (dis == n)
		if (eq.test | co.test) {
		    repe <- TRUE
		    break
		}
	    }
	    if (!repe) {
		if (verb & !is.na(w[1])) cat("+")
		Nvec <- Nvec + 1
		U[K,] <- vec
		Wu[K] <- vecW$weight
		if (log) { Wvec <- c(Wvec, vecW$weight) }
		K <- K + 1
	    }
	}
	if (verb) {
	    if (Nvec == 1) { sfin <- "um" } else { sfin <- "a" }
	    cat("(", Nvec, " local maxim", sfin, ") ",sep="")
	}                                               # End local search

	if (verb) cat("Selection...")                    # Begin selection 
	Ncan <- Npop + Noff + Nmut + Nvec
	Wcan <- Wu[1:Ncan]
	Wcs <- sort(Wcan, decreasing = TRUE, index.return = TRUE)
	imax <- Wcs$ix[1]
	if (Wcan[imax] > Wmax) {
	    Wmax <- Wcan[imax]
	    Smax <- vG[U[imax,] == 1]
	    gmax <- i
	}
	P[1:elite,] <- U[Wcs$ix[1:elite],]
	psel <- exp(-beta*(Wmax - Wcan)/mean(Wcan))
	psel[Wcs$ix[1:elite]] <- 0
	isel <- sample(1:Ncan, Npop - elite, prob=psel)
	P[(elite+1):Npop,] <- U[isel,]
	if (verb) cat("done [Wmax = ", Wmax, "]\n",sep="") # End selection

	## Preparing data for next generation
	U <- U0
	U[1:Npop,] <- P
	Wu <- Wu0
	for (i in 1:Npop) {
	    Wo <- compute_cut_weight(vG[P[i,] == 1], n, eG, w)
	    Wu[i] <- Wo$weight
	    if (log) { Wsel <- c(Wsel, Wo$weight)}
	}
	K <- Npop + 1
    }
    Wo <- compute_cut_weight(Smax, n, eG, w, return.cut = TRUE)
    if (log) {
	return(list(set = Smax, size = Wo$size, weight = Wo$weight,
		    cut = Wo$cut, generation = gmax, population = P,
		    Wini = Wini, Woff = Woff, Wmut = Wmut,
		    Wvec = Wvec, Wsel = Wsel))
    } else {
	return(list(set = Smax, size = Wo$size, weight = Wo$weight,
		    cut = Wo$cut, generation = gmax, population = P))
    }
}
