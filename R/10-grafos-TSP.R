##' Nearest neighbor heuristic tour-building algorithm for the
##' Traveling Salesperson Problem
##'
##' Starting from a vertex, the algorithm takes its nearest neighbor
##'     and incorporates it to the tour, repeating until the tour is
##'     complete.  The result is dependent of the initial vertex.
##'     This algorithm is very efficient but its output can be very
##'     far from the minimum.
##' 
##' @title Building a tour for a TSP using the nearest neighbor
##'     heuristic
##' @param d Distance matrix of the TSP.
##' @param n Number of vertices of the TSP complete graph.
##' @param v0 Starting vertex.  Valid values are integers between 1
##'     and n.
##' @return A list with two components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, and $distance contains the value of the distance
##'     covered by the tour.
##' @author Cesar Asensio
##' @seealso [build_tour_nn_best] repeats this algorithm with all
##'     possible starting points, [compute_tour_distance] computes
##'     tour distances, [compute_distance_matrix] computes a distance
##'     matrix, [plot_tour] plots a tour, [build_tour_greedy]
##'     constructs a tour using the greedy heuristic.
##' @examples
##' ## Regular example with obvious solution (minimum distance 48)
##' m <- 10   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' b <- build_tour_nn(d, n, 1)
##' b$distance    # Distance 50
##' plot_tour(z,b)
##' b <- build_tour_nn(d, n, 5)
##' b$distance    # Distance 52.38
##' plot_tour(z,b)
##'
##' ## Random points
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' b <- build_tour_nn(d, n, 1)
##' b$distance    # Distance 46.4088
##' plot_tour(z,b)
##' b <- build_tour_nn(d, n, 9)
##' b$distance    # Distance 36.7417
##' plot_tour(z,b)
##'
##' @encoding UTF-8
##' @md 
##' @export 
build_tour_nn <- function(d, n, v0) {
    h <- c(v0)
    v1 <- v0
    while (length(h) < n) {
	b <- d[v1,]
	b[v1] <- Inf
	b[h] <- rep(Inf, length(h))
	v1 <- which.min(b)
	h <- c(h, v1)
    }
    list(tour = h, distance = compute_tour_distance(h, d))
}

##' Nearest neighbor heuristic tour-building algorithm for the
##' Traveling Salesperson Problem - Better starting point
##'
##' It applies the nearest neighbor heuristic with all possible
##'     starting vertices, retaining the best tour returned by
##'     [build_tour_nn].
##'
##' @title Build a tour for a TSP using the best nearest neighbor
##'     heuristic
##' @param d Distance matrix of the TSP.
##' @param n Number of vertices of the TSP complete graph.
##' @return A list with four components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, $distance contains the value of the distance
##'     covered by the tour, $start contains the better starting
##'     vertex found, and $Lall contains the distances found by
##'     starting from each vertex.
##' @author Cesar Asensio
##' @seealso [build_tour_nn] nearest neighbor heuristic with a single
##'     starting point, [compute_tour_distance] computes tour
##'     distances, [compute_distance_matrix] computes a distance
##'     matrix, [plot_tour] plots a tour.
##' @examples 
##' ## Regular example with obvious solution (minimum distance 48)
##' m <- 10   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' b <- build_tour_nn_best(d, n)
##' b$distance    # Distance 48.6055
##' b$start       # Vertex 12
##' plot_tour(z,b)
##'
##' ## Random points
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' b <- build_tour_nn_best(d, n)
##' b$distance    # Distance 36.075
##' b$start       # Vertex 13
##' plot_tour(z,b)
##'
##' @encoding UTF-8
##' @md 
##' @export 
build_tour_nn_best <- function(d, n) {
    Lmin <- Inf
    vmin <- 0
    hmin <- Lall <- rep(0, n)
    for (v in 1:n) {
	b <- build_tour_nn(d, n, v)
	h <- b$tour
	L <- b$distance
	Lall[v] <- L
	if (L < Lmin) {
	    Lmin <- L
	    vmin <- v
	    hmin <- h
	}
    }
    list(tour = hmin, distance = Lmin, start = vmin, Lall = Lall)
}
##' It computes the distance covered by a tour in a Traveling Salesman Problem
##'
##' This function simply add the distances in a distance matrix
##'     indicated by a vertex sequence defining a tour.  It takes into
##'     account that, in a tour, the last vertex is joined to the
##'     first one by an edge, and adds its distance to the result,
##'     unlike [compute_path_distance].
##' 
##' @title Compute the distance of a TSP tour
##' @param h A tour specified by a vertex sequence
##' @param d Distance matrix to use
##' @return The tour distance \deqn{d(\{v_1,...,v_n\}) = \sum_{j=1}^n
##'     d(v_j,v_{(j mod n) + 1}).}
##' @author Cesar Asensio
##' @seealso [build_tour_nn] nearest neighbor heuristic with a single
##'     starting point, [build_tour_nn_best] repeats the previous
##'     algorithm with all possible starting points,
##'     [compute_distance_matrix] computes a distance matrix,
##'     [compute_path_distance] computes path distances, [plot_tour]
##'     plots a tour.
##' @examples
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' h <- sample(1:n)              # A random tour
##' compute_tour_distance(h, d)   # 114.58
##'
##' @encoding UTF-8
##' @md 
##' @export 
compute_tour_distance <- function(h, d) {
    t <- 0
    m <- length(h)
    for (i in 2:m) {
	t <- t + d[h[i-1], h[i]]
    }
    t <- t + d[h[m], h[1]]
    t
}

##' It computes the distance covered by a path in a Traveling Salesman Problem
##'
##' This function simply add the distances in a distance matrix
##'     indicated by a vertex sequence defining a path.  It takes into
##'     account that, in a path, the last vertex is \strong{not}
##'     joined to the first one by an edge, unlike [compute_tour_distance].
##' 
##' @title Compute the distance of a TSP path
##' @param h A path specified by a vertex sequence
##' @param d Distance matrix to use
##' @return The path distance \deqn{d(\{v_1,...,v_n\}) = \sum_{j=1}^{n-1}
##'     d(v_j,v_{j + 1}).}
##' @author Cesar Asensio
##' @seealso [build_tour_nn] nearest neighbor heuristic with a single
##'     starting point, [build_tour_nn_best] repeats the previous
##'     algorithm with all possible starting points,
##'     [compute_distance_matrix] computes a distance matrix,
##'     [compute_tour_distance] computes tour distances, [plot_tour]
##'     plots a tour.
##' @examples
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' h <- sample(1:n)              # A random tour
##' compute_path_distance(h, d)   # 107.246
##' compute_tour_distance(h, d) - compute_path_distance(h, d) - d[h[1], h[n]]
##'
##' @encoding UTF-8
##' @md 
##' @export 
compute_path_distance <- function(h, d) {
    t <- 0
    m <- length(h)
    for (i in 2:m) {
	t <- t + d[h[i-1], h[i]]
    }
    t
}

##' It computes the distance matrix of a set of \eqn{n}
##'     two-dimensional points given by a \eqn{n\times 2} matrix using
##'     the distance-\eqn{p} with \eqn{p=2} by default.
##'
##' Given a set of \eqn{n} points \eqn{\{z_j\}_{j=1,...,n}}, the
##'     distance matrix is a \eqn{n\times n} symmetric matrix with
##'     matrix elements \deqn{d_{ij} = d(z_i,z_j)} computed using the
##'     distance-\eqn{p} given by \deqn{d_p(x,y) = \left(\sum_i
##'     (x_i-y_i)^p\right)^{\frac{1}{p}}}.
##' 
##' @title \eqn{p}-distance matrix computation
##' @param z A \eqn{n\times 2} matrix with the two-dimensional points
##' @param p The \eqn{p} parameter of the distance-\eqn{p}.  It
##'     defaults to 2.
##' @return The distance-\eqn{p} of the points.
##' @author Cesar Asensio
##' @seealso [compute_p_distance] computes the distance-\eqn{p},
##'     [compute_tour_distance] computes tour distances.  A distance
##'     matrix can also be computed using [dist]. 
##' @examples 
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##'
##' @encoding UTF-8
##' @md 
##' @export 
compute_distance_matrix <- function(z, p = 2) {
    n <- nrow(z)
    d <- matrix(rep(0, n^2), ncol = n)
    for (i in 1:(n-1)) {
	for (j in (i+1):n) {
	    d[i,j] <- d[j,i] <- compute_p_distance(z[i,], z[j,], p)
	}
    }
    d
}

##' It computes the distance-\eqn{p} between two-dimensional points.
##'
##' The distance-\eqn{p} is defined by \deqn{d_p(x,y) = \left(\sum_i
##'     (x_i-y_i)^p\right)^{\frac{1}{p}}}.
##' 
##' @title Distance-p between two-dimensional points
##' @param x A two-dimensional point
##' @param y A two-dimensional point
##' @param p The \eqn{p}-parameter of the distance-\eqn{p}
##' @return The distance-\eqn{p} between points \eqn{x} and \eqn{y}.
##' @author Cesar Asensio
##' @seealso [compute_distance_matrix] computes the distance matrix of
##'     a set of two-dimensional points, [compute_tour_distance]
##'     computes tour distances.
##' @examples 
##' compute_p_distance(c(1,2),c(3,4))      # 2.8284
##' compute_p_distance(c(1,2),c(3,4),p=1)  # 4
##'
##' @encoding UTF-8
##' @md 
##' @export 
compute_p_distance <- function(x, y, p = 2) {
    (abs(x[1] - y[1])^p + abs(x[2] - y[2])^p)^(1/p)
}

##' Plotting tours constructed by tour-building routines for TSP
##'
##' It plots the two-dimensional cities of a TSP and a tour among them
##'     for visualization purposes.  No aesthetically appealing effort
##'     has been invested in this function.
##' 
##' @title TSP tour simple plotting
##' @param z Set of points of a TSP
##' @param h List with $tour and $distance components returned from a
##'     TSP tour building algorithm
##' @param ... Parameters to be passed to [plot]
##' @return This function is called by its side effect.
##' @importFrom graphics lines
##' @importFrom graphics text
##' @author Cesar Asensio
##' @seealso [build_tour_nn] nearest neighbor heuristic with a single
##'     starting point, [build_tour_nn_best] repeats the previous
##'     algorithm with all possible starting points,
##'     [compute_distance_matrix] computes the distance matrix of a
##'     set of two-dimensional points.
##' @examples 
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' b <- build_tour_nn_best(d, n)
##' plot_tour(z,b)
##'
##' @encoding UTF-8
##' @md 
##' @export 
plot_tour <- function(z,h,...) {
    plot(z[,1],z[,2], pch=19, xlim=c(0,max(z[,1])+1),
	 main=paste("Distance =", h$distance),...)
    N <- nrow(z)
    for(i in 1:N) {
	text(z[i,1],z[i,2], bquote(v[.(i)]), pos=4, col="gray60")
    }
    m <- length(h$tour)
    for (i in 2:m) {
	lines(z[h$tour[c(i-1,i)],1], z[h$tour[c(i-1,i)],2],
	      lwd=1.5, lty="dashed", col="red3")
    }
    lines(z[h$tour[c(1,m)],1], z[h$tour[c(1,m)],2],
	  lwd=1.5, lty="dashed", col="red3")    
}

##' Greedy heuristic tour-building algorithm for the Traveling
##' Salesperson Problem
##'
##' The greedy heuristic begins by sorting the edges by increasing
##'     distance.  The tour is constructed by adding an edge under the
##'     condition that the final tour is a connected spanning cycle. 
##'
##' @title Building a tour for a TSP using the greedy heuristic
##' @param d Distance matrix of the TSP.
##' @param n Number of vertices of the TSP complete graph.
##' @return A list with two components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, and $distance contains the value of the distance
##'     covered by the tour.
##' @author Cesar Asensio
##' @seealso [build_tour_nn] uses the nearest heighbor heuristic,
##'     [build_tour_nn_best] repeats the previous algorithm with all
##'     possible starting points, [compute_tour_distance] computes
##'     tour distances, [compute_distance_matrix] computes a distance
##'     matrix, [plot_tour] plots a tour, [build_tour_2tree]
##'     constructs a tour using the double tree 2-factor approximation
##'     algorithm.
##' @examples
##' ## Regular example with obvious solution (minimum distance 48)
##' m <- 10   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' b <- build_tour_greedy(d, n)
##' b$distance    # Distance 50
##' plot_tour(z,b)
##' 
##' ## Random points
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' b <- build_tour_greedy(d, n)
##' b$distance    # Distance 36.075
##' plot_tour(z,b)
##' 
##' @encoding UTF-8
##' @md 
##' @export 
build_tour_greedy <- function(d, n) {
    q <- n*(n-1)/2
    qlist <- matrix(0,nrow=q,ncol=3)
    k <- 0

    ## Transform the distance matrix to a list:
    for (i in 1:(n-1)) {
	for (j in (i+1):n) {
	    k <- k+1
	    qlist[k,] <- c(j,i,d[i,j])
	}
    }

    ## Sort the distance list:
    qls <- sort(qlist[,3],index.return=TRUE)
    qlist <- qlist[qls$ix,]

    ## Form a tour by adding edges from the sorted distance list 
    qacep <- rep(FALSE,q)     # Edges forming the tour
    gras <- cols <- rep(0,n)  # Degrees and colors of vertices
    curcol <- 1
    for (k in 1:q) {
	v <- qlist[k, 1]
	w <- qlist[k, 2]
	dv <- gras[v] + 1
	dw <- gras[w] + 1
	if (dv > 2 | dw > 2) { next }
	if (sum(cols == curcol) > 0) { curcol <- curcol + 1 }
	incor <- FALSE
	if (cols[v] == 0 & cols[w] == 0) {
	    cols[v] <- cols[w] <- curcol
	    incor <- TRUE
	}
	if (cols[v] == 0 & cols[w] > 0)  {
	    cols[v] <- cols[w]
	    incor <- TRUE
	}
	if (cols[v] > 0 & cols[w] == 0)  {
	    cols[w] <- cols[v]
	    incor <- TRUE
	}
	if ((cols[v] > 0 & cols[w] > 0) & cols[v] != cols[w]) {
	    cols[cols == cols[v]] <- cols[w] 
	    incor <- TRUE
	}
	if ((cols[v] > 0 & cols[w] > 0) & cols[v] == cols[w]) {
	    if (sum(qacep) == n - 1) { incor <- TRUE }
	}
	if (incor) {
	    qacep[k] <- TRUE
	    gras[v] <- dv
	    gras[w] <- dw
	}
	if (sum(qacep) == n) { break }
    }

    ## We transform the edge list of the tour into a vertex sequence:
    tour <- qlist[qacep,]
    L <- sum(tour[,3])
    h <- rep(0,n)
    h[1:2] <- tour[1, 1:2]
    v <- h[2]
    tour <- tour[-1,]
    for (i in 3:n) {
	if (is.matrix(tour)) {
	    rw <- which(tour[,1]==v | tour[,2]==v)
	    trw <- tour[rw, 1:2]
	} else {
	    trw <- tour[1:2]
	}
	v <- trw[trw != v]
	h[i] <- v
	tour <- tour[-rw,]
    }

    list(tour = h, distance = L)
}

##' Double-tree heuristic tour-building algorithm for the Traveling
##' Salesperson Problem
##'
##' The \strong{double-tree} heuristic is a 2-factor approximation
##'     algorithm which begins by forming a minimum distance spanning
##'     tree, then it forms the double-tree by doubling each edge of
##'     the spanning tree.  The double tree is Eulerian, so an
##'     Eulerian walk can be computed, which gives a well-defined
##'     order of visiting the cities of the problem, thereby yielding
##'     the tour.
##'
##' In practice, this algorithm performs poorly when compared with
##'     another simple heuristics such as nearest-neighbor or
##'     insertion methods.
##' 
##' @title Double-tree heuristic for TSP
##' @param d Distance matrix defining the TSP instance
##' @param n Number of cities to consider with respect to the distance matrix
##' @param v0 Initial vertex to find the eulerian walk; it defaults to 1.
##' @return A list with two components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, and $distance contains the value of the distance
##'     covered by the tour.
##' @author Cesar Asensio
##' @seealso [build_tour_nn] uses the nearest heighbor heuristic,
##'     [build_tour_nn_best] repeats the previous algorithm with all
##'     possible starting points, [compute_tour_distance] computes
##'     tour distances, [compute_distance_matrix] computes a distance
##'     matrix, [plot_tour] plots a tour, [find_euler] finds an
##'     Eulerian walk.     
##' @examples
##' ## Regular example with obvious solution (minimum distance 48)
##' m <- 10   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' b <- build_tour_2tree(d, n)
##' b$distance    # Distance 57.86
##' plot_tour(z,b)
##' 
##' ## Random points
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' b <- build_tour_2tree(d, n)
##' b$distance    # Distance 48.63
##' plot_tour(z,b)
##' 
##' @encoding UTF-8
##' @md 
##' @export 
build_tour_2tree <- function(d, n, v0 = 1) {
    ## Form a weighted complete graph using input distances
    Kn <- make_full_graph(n) 
    EKn <- as_edgelist(Kn)
    q <- gsize(Kn)
    dl <- rep(0,q)
    for (e in 1:q) dl[e] <- d[EKn[e,1],EKn[e,2]]
    Kn <- set_edge_attr(Kn,"weight",value=dl)

    ## Compute minimum distance spanning tree
    T <- mst(Kn)

    ## Form the double tree
    T2 <- add_edges(T,as.vector(t(as_edgelist(T))))

    ## Find an Eulerian walk on the double tree
    euT2 <- find_euler(T2, v0)

    ## Find the order of visited cities by the Eulerian walk
    vis <- euT2$walk[,1]
    h <- rep(0,n)
    m <- rep(FALSE,n)
    i <- 1
    for (j in 1:(2*n-2)) {
	v <- vis[j]
	if (m[v]) {
	    next
	} else {
	    m[v] <- TRUE
	    h[i] <- v
	    i <- i+1
	}
    }

    list(tour = h, distance = compute_tour_distance(h,d))
}

##' 2-opt heuristic tour-improving algorithm for the Traveling
##' Salesperson Problem
##'
##' It applies the 2-opt algorithm to a starting tour of a TSP
##'     instance until no further improvement can be found.  The tour
##'     thus improved is a 2-opt local minimum.
##'
##' The 2-opt algorithm consists of applying all possible
##'     2-interchanges on the starting tour.  Informally, a
##'     2-interchange is the operation of cutting the tour in two
##'     pieces (by removing two nonincident edges) and gluing the
##'     pieces together to form a new tour by interchanging the
##'     endpoints.
##' 
##' @title Tour improving for a TSP using the 2-opt heuristic
##' @param d Distance matrix of the TSP.
##' @param n Number of vertices of the TSP complete graph.
##' @param C Starting tour to be improved.
##' @return A list with two components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, $distance contains the value of the distance
##'     covered by the tour.
##' @author Cesar Asensio
##' @seealso [improve_tour_3opt] improves a tour using the 3-opt
##'     algorithm, [build_tour_nn_best] nearest neighbor heuristic,
##'     [build_tour_2tree] double-tree heuristic,
##'     [compute_tour_distance] computes tour distances,
##'     [compute_distance_matrix] computes a distance matrix,
##'     [plot_tour] plots a tour.
##' @examples 
##' ## Regular example with obvious solution (minimum distance 48)
##' m <- 10   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' b <- build_tour_2tree(d, n)
##' b$distance    # Distance 57.868
##' bi <- improve_tour_2opt(d, n, b$tour)
##' bi$distance   # Distance 48 (optimum)
##' plot_tour(z,b)
##' plot_tour(z,bi)
##'
##' ## Random points
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' b <- build_tour_2tree(d, n)
##' b$distance    # Distance 48.639
##' bi <- improve_tour_2opt(d, n, b$tour)
##' bi$distance   # Distance 37.351
##' plot_tour(z,b)
##' plot_tour(z,bi)
##'
##' @encoding UTF-8
##' @md 
##' @export 
improve_tour_2opt <- function(d,n,C) {
    Lc <- compute_tour_distance(C,d)
    Ln <- Lmin <- Lmin.old <- Lc
    Cmin <- C
    while(TRUE) {
	for (i in 1:(n-2)) {
	    ## a: pos i, b: pos i+1 (1 <= i <= n-2)
	    if (i == 1) { m <- n-1 } else { m <- n }
	    for (j in (i+2):m) {
		## c: pos j, b: pos j+1 (i+2 <= j <= n+1 mod n)
		if (j < n) {
		    c1 <- C[(i+1):j]
		    c2 <- c(C[(j+1):n],C[1:i])
		} else {
		    c1 <- C[1:i]
		    c2 <- C[(i+1):n]
		}
		Cn <- c(c1,rev(c2))
		Ln <- compute_tour_distance(Cn,d)
		if (Ln < Lmin) {
		    Cmin <- Cn
		    Lmin <- Ln
		}
	    }
	}
	if (Lmin < Lmin.old) {
	    Lmin.old <- Lmin
	    C <- Cmin
	} else {
	    break
	}
    }
    list(tour = Cmin, distance = Lmin)
}

##' 3-opt heuristic tour-improving algorithm for the Traveling
##' Salesperson Problem
##'
##' It applies the 3-opt algorithm to a starting tour of a TSP
##'     instance until no further improvement can be found.  The tour
##'     thus improved is a 3-opt local minimum.
##'
##' The 3-opt algorithm consists of applying all possible
##'     3-interchanges on the starting tour.  A 3-interchange removes
##'     three non-indicent edges from the tour, leaving three pieces,
##'     and combine them to form a new tour by interchanging the
##'     endpoints in all possible ways and gluing them together by
##'     adding the missing edges.
##' 
##' @title Tour improving for a TSP using the 3-opt heuristic
##' @param d Distance matrix of the TSP.
##' @param n Number of vertices of the TSP complete graph.
##' @param C Starting tour to be improved.
##' @return A list with two components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, $distance contains the value of the distance
##'     covered by the tour.
##' @author Cesar Asensio
##' @seealso [improve_tour_2opt] improves a tour using the 2-opt
##'     algorithm, [build_tour_nn_best] nearest neighbor heuristic,
##'     [build_tour_2tree] double-tree heuristic,
##'     [compute_tour_distance] computes tour distances,
##'     [compute_distance_matrix] computes a distance matrix,
##'     [plot_tour] plots a tour.
##' @examples 
##' ## Regular example with obvious solution (minimum distance 32)
##' m <- 6   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' b <- build_tour_2tree(d, n)
##' b$distance    # Distance 38.43328
##' bi <- improve_tour_3opt(d, n, b$tour)
##' bi$distance   # Distance 32 (optimum)
##' plot_tour(z,b)
##' plot_tour(z,bi)
##'
##' ## Random points
##' set.seed(1)
##' n <- 15
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' b <- build_tour_2tree(d, n)
##' b$distance    # Distance 45.788
##' bi <- improve_tour_3opt(d, n, b$tour)
##' bi$distance   # Distance 32.48669
##' plot_tour(z,b)
##' plot_tour(z,bi)
##'
##' @encoding UTF-8
##' @md 
##' @export 
improve_tour_3opt <- function(d,n,C) {
    Lc <- compute_tour_distance(C,d)
    Ln <- Lmin <- Lmin.old <- Lc
    Cmin <- C
    while(TRUE) {
	for (i in 1:(n-4)) {
	    for (j in (i+2):(n-2)) {
		if (j == 1) { m <- n-1 } else { m <- n }
		for (k in (j+2):m) {
		    if (k < n) {
			c1 <- C[(i+1):j]
			c2 <- C[(j+1):k]
			c3 <- c(C[(k+1):n],C[1:i])
		    } else {
			c1 <- C[1:i]
			c2 <- C[(i+1):j]
			c3 <- C[(j+1):n]
		    }
		    Cn1 <- c(rev(c1),c2,c3)
		    Cn2 <- c(c1,rev(c2),c3)
		    Cn3 <- c(c1,c2,rev(c3))
		    Cn4 <- c(c1,c3,c2)
		    Cn5 <- c(c1,rev(c3),c2)
		    Cn6 <- c(c1,c3,rev(c2))
		    Cn7 <- c(c1,rev(c2),rev(c3))
		    for (Cns in c("Cn1", "Cn2", "Cn3", "Cn4",
				  "Cn5", "Cn6", "Cn7")) {
			Cn <- get(Cns)
			Ln <- compute_tour_distance(Cn,d)
			if (Ln < Lmin) {
			    Cmin <- Cn
			    Lmin <- Ln
			}
		    }
		}
	    }
	}
	if (Lmin < Lmin.old) {
	    Lmin.old <- Lmin
	    C <- Cmin
	} else {
	    break
	}
    }
    list(tour = Cmin, distance = Lmin)
}

## The following routine is a very poor version of LK using a simple
## neighborhood: firstly I used transpositions, but the final result
## were not even 2-opt! Then I used fixed 2-opt, but then the routine
## is no better than 2-opt.  This raises the question: Should it be
## explained??
## ------------------------------------------------------------------


##' Lin-Kernighan heuristic tour-improving algorithm for the Traveling
##'     Salesperson Problem using fixed 2-opt instead of variable
##'     k-opt exchanges.
##'
##' It applies a version of the core Lin-Kernighan algorithm to a
##'     starting tour of a TSP instance until no further improvement
##'     can be found.  The tour thus improved is a local minimum.
##'
##' The Lin-Kernighan algorithm implemented here is based on the core
##'     routine described in the reference below.  It is provided here
##'     as an example of a local search routine which can be embedded
##'     in larger search strategies.  However, instead of using
##'     variable k-opt moves to improve the tour, it uses 2-exhanges
##'     only, which is far easier to program.  Tours improved with
##'     this technique are of course 2-opt.
##'
##' The TSP library provides an interface to the Lin-Kernighan
##'     algorithm with all its available improvements in the external
##'     program Concorde, which should be installed separately.
##' 
##' @title Tour improving for a TSP using a poor version of the Lin-Kernighan heuristic
##' @param d Distance matrix of the TSP.
##' @param n Number of vertices of the TSP complete graph.
##' @param C Starting tour to be improved.
##' @param try Number of tries before quitting.
##' @return A list with two components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, $distance contains the value of the distance
##'     covered by the tour.
##' @references Hromkovic \emph{Algorithmics for Hard Problems} (2004)
##' @author Cesar Asensio
##' @seealso [improve_tour_2opt] improves a tour using the 2-opt
##'     algorithm, [improve_tour_3opt] improves a tour using the 3-opt
##'     algorithm, [build_tour_nn_best] nearest neighbor heuristic,
##'     [build_tour_2tree] double-tree heuristic,
##'     [compute_tour_distance] computes tour distances,
##'     [compute_distance_matrix] computes a distance matrix,
##'     [plot_tour] plots a tour.
##' @examples 
##' ## Regular example with obvious solution (minimum distance 48)
##' m <- 10   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' b <- build_tour_2tree(d, n)
##' b$distance    # Distance 57.868
##' bi <- improve_tour_LinKer(d, n, b$tour)
##' bi$distance   # Distance 48 (optimum)
##' plot_tour(z,b)
##' plot_tour(z,bi)
##'
##' ## Random points
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' b <- build_tour_2tree(d, n)
##' b$distance    # Distance 48.639
##' bi <- improve_tour_LinKer(d, n, b$tour)
##' bi$distance   # Distance 37.351 (2-opt)
##' plot_tour(z,b)
##' plot_tour(z,bi)
##'
##' @encoding UTF-8
##' @md 
##' @export 
improve_tour_LinKer <- function(d,n,C,try=5) {
    Lmin <- Lact <- Ltra <- Lnei <- Lmin.old <- compute_tour_distance(C,d)
    Cmin <- Cact <- Ctra <- Cnei <- C
    T <- T0 <- matrix(0,nrow=n,ncol=n)
    tk <- tmin <- c(0,0)
    max1inT <- n*(n-3)/2 # n*(n-1)/2
    num1inT <- 0
    ntry.neg <- ntry.zero <- 0
    K <- 0
    while(TRUE) {
	K <- K+1
	##cat(rep("-",30,),"[K = ",K,"]",rep("-",30),"\n",sep="")
	## for (i in 1:(n-1)) {       # Scan transp neighborhood
	##     for (j in (i+1):n) {
	## 	if (T[i,j]==1) { next }   # Only unmarked transps
	## 	tk <- c(i,j)
	## 	Ctra <- transpose_tour(Cact,tk)
	## 	Ltra <- compute_tour_distance(Ctra,d)
	## 	if (Ltra < Lnei) {   # Retain best in neighborhood
	## 	    tmin <- tk
	## 	    Cnei <- Ctra
	## 	    Lnei <- Ltra
	## 	}
	##     }
	## }
	for (i in 1:(n-2)) {
	    if (i == 1) { m <- n-1 } else { m <- n }
	    for (j in (i+2):m) {
		if (j < n) {
		    c1 <- Cact[(i+1):j]
		    c2 <- c(Cact[(j+1):n],Cact[1:i])
		} else {
		    c1 <- Cact[1:i]
		    c2 <- Cact[(i+1):n]
		}
		Ctra <- c(c1,rev(c2))
		Ltra <- compute_tour_distance(Ctra,d)
		if (Ltra < Lnei) {
		    tmin <- c(i,j)
		    Cnei <- Ctra
		    Lnei <- Ltra
		}
	    }
	}

	if (Lnei > Lact) {            # Should we exit?
	    ntry.neg <- ntry.neg + 1
	    ##cat("Iteraciones  sin encontrar vecinos mejores:", ntry.neg,"\n")
	    if (ntry.neg == try) { break }
	} else {
	    ntry.neg <- 0
	}
	T[tmin[1],tmin[2]] <- 1    # Mark transp already made
	num1inT <- num1inT + 1
	##cat("Ciclo activo: Lact = ", Lact, "\n")
	##cat("Mejor vecino: Lnei = ", Lnei, "\n")
	if (Lnei < Lmin) {        # It is better?
	    ##cat("Mejorado desde Lmin = ", Lmin, " hasta Lmin = " ,Lnei, "\n")
	    Lmin <- Lnei
	    Cmin <- Cnei
	    Cact <- Cmin
	    Lact <- Lmin
	    Lnei <- Inf
	}
	if (num1inT == max1inT & Lmin < Lmin.old) {
	    ##cat("Fin de pasada: distancia: ", Lmin,"\n")
	    ##cat("Ciclo: ",Cmin,"\n")
	    Cact <- Cmin
	    Lmin.old <- Lact <- Lmin
	    Lnei <- Inf
	    T <- T0
	    num1inT <- 0
	    ntry.zero <- 0
	    next
	}
	if (num1inT == max1inT & Lmin == Lmin.old) {
	    ntry.zero <- ntry.zero + 1
	    ##cat("Iteraciones sin mejora:", ntry.zero,"\n")
	    if (ntry.zero == try) { break }
	    Cact <- Cnei
	    Lact <- Lnei
	    Lnei <- Inf
	    T <- T0
	    num1inT <- 0
	}
    }

    list(tour = Cmin, distance = Lmin)
}

transpose_tour <- function(C,tr) {
    vtemp <- C[tr[1]]
    C[tr[1]] <- C[tr[2]]
    C[tr[2]] <- vtemp
    C
}


##' Distance gain when two cities in a TSP tour are interchanged, that
##'     is, the neighbors of the first become the neighbors of the
##'     second and vice versa.  It is used to detect favorable moves
##'     in a Lin-Kernighan-based routine for the TSP.
##'
##' It computes the gain in distance when interchanging two cities in
##'     a tour.  The transformation is akin to a 2-interchange; in
##'     fact, if the transposed vertices are neighbors in the tour or
##'     share a common neighbor, the transposition is a
##'     2-interchange.  If the transposed vertices in the tour do not
##'     share any neighbors, then the transposition is a pair of
##'     2-interchanges.
##'
##' This gain is used in [improve_tour_LinKer], where the
##'     transposition neighborhood is used instead of the variable
##'     k-opt neighborhood for simplicity.
##' 
##' @title Distance gain when transposing two cities in a tour
##' @param C Tour represented as a non-repeated vertex sequence.
##'     Equivalently, a permutation of the sequence from 1 to length(C). 
##' @param tr Transposition, represented as a pair of indices between
##'     1 and length(C).
##' @param d Distance matrix.
##' @return The gain in distance after performing transposition tr in
##'     tour C with distance matrix d.
##' @author Cesar Asensio
##' @seealso [improve_tour_LinKer], a where this function is used.
##' @examples 
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' compute_gain_transp(sample(n),c(4,23),d)  # -6.661
##' compute_gain_transp(sample(n),c(17,3),d)  #  4.698
##' 
##' @encoding UTF-8
##' @md 
##' @export 
compute_gain_transp <- function(C,tr,d) {
    n <- length(C)
    t1 <- min(tr)
    t2 <- max(tr)
    u <- C[t1]
    if (t1==1) {ua <- C[n] } else { ua <- C[t1 - 1] }
    if (t1==n) {up <- C[1] } else { up <- C[t1 + 1] }
    v <- C[t2]
    if (t2==1) {va <- C[n] } else { va <- C[t2 - 1] }
    if (t2==n) {vp <- C[1] } else { vp <- C[t2 + 1] }
    if (u==va | up==va) {
	g <- d[u,vp] + d[v,ua] - d[v,vp] - d[u,ua] 
    } else {
	g <- d[u,vp] + d[u,va] + d[v,up] + d[v,ua] - d[v,vp] - d[v,va] - d[u,up] - d[u,ua]
    }
    return(-g)
}

##' Random heuristic algorithm for TSP which performs chained 2-opt
##'     local search with multiple, random starting tours.
##'
##' Chained local search consists of starting with a random tour,
##'     improving it using 2-opt, and then perturb it using a random
##'     4-exchange.  The result is 2-optimized again, and then
##'     4-exchanged...  This sequence of chained
##'     2-optimizations/perturbations is repeated Nper times for each
##'     random starting tour.  The entire process is repeated Nit
##'     times, drawing a fresh random tour each iteration.
##'
##' The purpose of supplement the deterministic 2-opt algorithm with
##'     random additions (random starting point and random 4-exchange)
##'     is escaping from the 2-opt local minima.  Of course, more
##'     iterations and more perturbations might lower the result, but
##'     recall that no random algorithm can guarantee to find the
##'     optimum in a reasonable amount of time.
##'
##' This technique is most often applied in conjunction with the
##'     Lin-Kernighan local search heuristic.
##'
##' It should be warned that this algorithm calls Nper*Nit times the
##'     routine [improve_tour_2opt], and thus it is not especially efficient.
##' 
##' @title Chained 2-opt search with multiple, random starting tours
##' @param d Distance matrix of the TSP instance.
##' @param n Number of vertices of the TSP complete graph.
##' @param Nit Number of iterations of the algorithm, see details.
##' @param Nper Number of chained perturbations of 2-opt minima, see details.
##' @param log Boolean: Whether the algorithm should record the
##'     distances of the tours it finds during execution.  It
##'     defaults to FALSE.
##' @return A list with two components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, $distance contains the value of the distance
##'     covered by the tour.
##' @references Cook et al. \emph{Combinatorial Optimization} (1998)
##' @author Cesar Asensio
##' @seealso [perturb_tour_4exc] transforms a tour using a random
##'     4-exchange, [improve_tour_2opt] improves a tour using the 2-opt
##'     algorithm, [improve_tour_3opt] improves a tour using the 3-opt
##'     algorithm, [build_tour_nn_best] nearest neighbor heuristic,
##'     [build_tour_2tree] double-tree heuristic,
##'     [compute_tour_distance] computes tour distances,
##'     [compute_distance_matrix] computes a distance matrix,
##'     [plot_tour] plots a tour.
##' @examples 
##' ## Regular example with obvious solution (minimum distance 32)
##' m <- 6   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' bc <- search_tour_chain2opt(d, n, 5, 3)
##' bc     # Distance 48
##' plot_tour(z,bc)
##'
##' ## Random points
##' set.seed(1)
##' n <- 15
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' bc <- search_tour_chain2opt(d, n, 5, 3)
##' bc     # Distance 32.48669
##' plot_tour(z,bc)
##'
##' @encoding UTF-8
##' @md 
##' @export 
search_tour_chain2opt <- function(d, n, Nit, Nper, log=FALSE) {
    C0 <- 1:n
    Lmin <- Inf
    if (log) { Lran <- Lloc <- Lper <- c() }
    for (i in 1:Nit) {
	C <- sample(C0,n)
	L <- compute_tour_distance(C, d)
	if (L < Lmin) {
	    Cmin <- C
	    Lmin <- L
	}
	if (log) { Lran <- c(Lran,L) }
	for (j in 1:Nper) {
	    mej <- improve_tour_2opt(d, n, C)
	    C2 <- mej$tour
	    L2 <- mej$distance
	    if (L2 < Lmin) {
		Cmin <- C2
		Lmin <- L2
		cat("Iteration ", i, ", perturbation ", j, ", Lmin = ",
		    Lmin, "\n", sep="")
	    }
	    Cp <- perturb_tour_4exc(C2, C0, n)
	    Lp <- compute_tour_distance(Cp, d)
	    if (Lp < Lmin) {
		Cmin <- Cp
		Lmin <- Lp
	    }
	    if (log) {
		Lloc <- c(Lloc, L2) 
		Lper <- c(Lper, Lp)
	    }
	    C <- Cp
	}
    }
    if (log) {
	return(list(tour = Cmin, distance = Lmin,
		    Lran = Lran, Lloc = Lloc, Lper = Lper))
    } else {
	return(list(tour = Cmin, distance = Lmin))
    }
}

##' Previous, current, and next positions of a given index in a cycle.
##' 
##'
##' Given some position i in a n-lenght cycle, this function returns
##'     the triple c(i-1,i,i+1) taking into account that the next
##'     position of i=n is 1 and the previous position of i=1 is n.
##'     It is used to perform a 4-exchange in a cycle.
##' 
##' @title Previous, current, and next positions of a given index in a
##'     cycle.
##' @param i Position in a cycle
##' @param n Lenght of the cycle
##' @return A three component vector c(previous, current, next)
##' @author Cesar Asensio
##' @examples
##' neigh_index(6, 9)  # 5 6 7
##' neigh_index(9, 9)  # 8 9 1
##' neigh_index(1, 9)  # 9 1 2
##' 
##' @encoding UTF-8
##' @md 
##' @export 
neigh_index <- function(i, n) {
    if (i==1) { return(c(n,1,2)) }
    if (i==n) { return(c(n-1,n,1)) }
    return(c(i-1,i,i+1))
}

##' Next position to i in a cycle.
##'
##' In a cycle, the next slot to the i-th position is i+1 unless i=n.
##'     In this case, the next is 1.
##' 
##' @title Next position to i in a cycle
##' @param i Position in cycle 
##' @param n Lenght of cycle
##' @return The next position in cycle
##' @author Cesar Asensio
##' @examples
##' next_index(5, 7)   # 6
##' next_index(7, 7)   # 1
##' 
##' @encoding UTF-8
##' @md 
##' @export 
next_index <- function(i, n) {
    if (i==n) { return(1) } else { return(i+1) }
}

##' It performs a random 4-exchange transformation to a cycle.
##'
##' The transformation is carried out by randomly selecting four
##'     non-mutually incident edges from the cycle.  Upon eliminating
##'     these four edges, we obtain four pieces ci of the original
##'     cycle.  The 4-exchanged cycle is c1, c4, c3, c2.  This is a
##'     typical 4-exchange which cannot be constructed using
##'     2-exchanges and therefore it is used by local search routines
##'     as an escape from 2-opt local minima.
##' 
##' @title Random 4-exchange transformation
##' @param C Cycle to be 4-exchanged
##' @param V 1:n list, positions to draw from
##' @param n Number of vertices of the cycle
##' @return The 4-exchanged cycle.
##' @author Cesar Asensio
##' @examples
##' set.seed(1)
##' perturb_tour_4exc(1:9, 1:9, 9)   # 2 3 9 1 7 8 4 5 6
##' 
##' @encoding UTF-8
##' @md 
##' @export 
perturb_tour_4exc <- function(C, V, n) {
    p <- rep(1, n)
    i1 <- sample(V, 1)
    v1 <- neigh_index(i1, n)
    p[v1] <- 0
    i2 <- sample(V, 1, prob=p)
    v2 <- neigh_index(i2, n)
    p[v2] <- 0
    i3 <- sample(V, 1, prob=p)
    v3 <- neigh_index(i3, n)
    p[v3] <- 0
    i4 <- sample(V, 1, prob=p)
    is <- sort(c(i1,i2,i3,i4))
    v1 <- is[1];  u1 <- next_index(v1, n)
    v2 <- is[2];  u2 <- next_index(v2, n)
    v3 <- is[3];  u3 <- next_index(v3, n)
    v4 <- is[4];  u4 <- next_index(v4, n)
    c1 <- C[u1:v2]
    c2 <- C[u2:v3]
    c3 <- C[u3:v4]
    if (u4==1) {
	c4 <- C[1:v1]
    } else {
	c4 <- c(C[u4:n],C[1:v1])
    }
    return(c(c1,c4,c3,c2))
}

##' Genetic algorithm for TSP.  In addition to crossover and mutation,
##'     which are described below, the algorithm performs also 2-opt
##'     local search on offsprings and mutants.  In this way, this
##'     algorithm is at least as good as chained 2-opt search.
##'
##' The genetic algorithm consists of starting with a tour population,
##'     which can be provided by the user or can be random.  The
##'     initial population can be the output of a previous run of the
##'     genetic algorithm, thus allowing a chained execution.  Then
##'     the routine sequentially perform over the tours of the
##'     population the \strong{crossover}, \strong{mutation},
##'     \strong{local search} and \strong{selection} operations.
##'
##' The \strong{crossover} operation takes two tours and forms two
##'     offsprings trying to exploit the good structure of the
##'     parents; see [crossover_tours] for more information.
##'
##' The \strong{mutation} operation performs a "small" perturbation of
##'     each tour trying to escape from local optima.  It uses a
##'     random 4-exchange, see [perturb_tour_4exc] and
##'     [search_tour_chain2opt] for more information.
##'
##' The \strong{local search} operation takes the tours found by the
##'     crossover and mutation operations and improves them using the
##'     2-opt local search heuristic, see [improve_tour_2opt].  This
##'     makes this algorithm at least as good as chained local search,
##'     see [search_tour_chain2opt].
##'
##' The \strong{selection} operation is used when selecting pairs of
##'     parents for crossover and when selecting individuals to form
##'     the population for the next generation.  In both cases, it
##'     uses a probability exponential in the distance with rate
##'     parameter "beta", favouring the better fitted to be selected.
##'     Lower values of beta favours the inclusion of tours with worse
##'     fitting function values.  When selecting the next population,
##'     the selection uses \emph{elitism}, which is to save the best
##'     fitted individuals to the next generation; this is controlled
##'     with parameter "elite".
##'
##'
##' The usefulness of the crossover and mutation operations stems from
##'     its abitily to escape from the 2-opt local minima in a way
##'     akin to the perturbation used in chained local search
##'     [search_tour_chain2opt].  Of course, more iterations (Ngen)
##'     and larger populations (Npop) might lower the result, but
##'     recall that no random algorithm can guarantee to find the
##'     optimum of a given TSP instance.
##'
##' This algorithm calls many times the routines [crossover_tours],
##'     [improve_tour_2opt] and [perturb_tour_4exc]; therefore, it is
##'     not especially efficient when called on large problems or with
##'     high populations of many generations.  Input parameter "local"
##'     can be used to randomly select which tours will start local
##'     search, thus diminishing the run time of the algorithm.
##'     Please consider chaining the algorithm:  perform short runs,
##'     using the output of a run as the input of the next.
##' 
##' 
##' @title Genetic Algorithm for the TSP
##' @param d Distance matrix of the TSP instance.
##' @param n Number of vertices of the TSP complete graph.
##' @param Npop Population size.
##' @param Ngen Number of generations (iterations of the algorithm).
##' @param beta Control parameter of the crossing and selection
##'     probabilities.  It defaults to 1.
##' @param elite Number of better fitted individuals to pass on to the
##'     next generation.  It defaults to 2.
##' @param Pini Initial population.  If it is NA, a random initial
##'     population of Npop individuals is generated.  Otherwise, it
##'     should be a matrix; each row should be an individual (a
##'     permutation of the 1:n sequence) and then Npop is set to the
##'     number of rows of Pini.  This option allows to chain several
##'     runs of the genetic algorithm, which could be needed in the
##'     hardest cases.
##' @param local Average fraction of parents + offsprings + mutants
##'     that will be taken as starting tours by the local search
##'     algorithm [improve_tour_2opt].  It should be a number between
##'     0 and 1.  It defauls to 1.
##' @param verb Boolean to activate console echo.  It defaults to
##'     TRUE.
##' @param log Boolean to activate the recording of the distances of
##'     all tours found by the algorithm.  It defaults to FALSE.
##' @return A list with four components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, $distance contains the value of the distance
##'     covered by the tour, $generation contains the generation is
##'     which the minimum was found and $population contains the final
##'     tour population.  When log=TRUE, the output includes several
##'     lists of distances of tours found by the algorithm, separated
##'     by initial tours, offsprings, mutants, local minima and
##'     selected tours.
##' @references Hromkovic \emph{Algorithms for hard problems} (2004),
##'     Hartmann, Weigt, \emph{Phase transitions in combinatorial
##'     optimization problems} (2005).
##' @author Cesar Asensio
##' @importFrom stats runif
##' @seealso [crossover_tours] performs the crosover of two tours,
##'     [gauge_tour] transforms a tour into a canonical sequence for
##'     comparison, [search_tour_chain2opt] performs a chained 2-opt
##'     search, [perturb_tour_4exc] transforms a tour using a random
##'     4-exchange, [improve_tour_2opt] improves a tour using the
##'     2-opt algorithm, [improve_tour_3opt] improves a tour using the
##'     3-opt algorithm, [build_tour_nn_best] nearest neighbor
##'     heuristic, [build_tour_2tree] double-tree heuristic,
##'     [compute_tour_distance] computes tour distances,
##'     [compute_distance_matrix] computes a distance matrix,
##'     [plot_tour] plots a tour.
##' @examples 
##' ## Regular example with obvious solution (minimum distance 32)
##' m <- 6   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' bc <- search_tour_genetic(d, n, Npop = 5, Ngen = 3, local = 0.2)
##' bc     # Distance 32
##' plot_tour(z,bc)
##'
##' ## Random points
##' set.seed(1)
##' n <- 15
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' bg <- search_tour_genetic(d, n, 5, 3, local = 0.25)
##' bg     # Distance 32.48669
##' plot_tour(z,bg)
##'
##' @encoding UTF-8
##' @md 
##' @export 
search_tour_genetic <- function(d, n, Npop = 20, Ngen = 50,
				beta = 1, elite = 2, Pini = NA,
				local = 1, verb = TRUE, log = FALSE) {
    C0 <- 1:n

    ## Initial population
    if (is.na(Pini)) {
	welc <- "Generated"
	Pini <- matrix(0, nrow = Npop, ncol = n)
	for (i in 1:Npop) {
	    Cr <- sample(C0, n)
	    Cr <- gauge_tour(Cr, n)
	    Pini[i,] <- Cr
	}
    } else {
	welc <- "Received"
	Npop <- nrow(Pini)
    }

    ## Initialization
    Lmin <- Inf
    P <- Pini
    U  <- U0 <- matrix(0, ncol = n, nrow = Npop + Npop + 2*Npop + 4*Npop)
    U[1:Npop,] <- P
    Lu <- Lu0 <- rep(0,8*Npop)
    for (i in 1:Npop) {
	Lu[i] <- compute_tour_distance(P[i,], d)
    }
    K <- Npop + 1

    if (log) {
	Lini <- Lu[1:Npop]
	Loff <- Lmut <- Lvec <- Lsel <- c()
    }
    if (verb) cat(welc, " initial population (Npop = ", Npop, ")\n", sep="")

    ## Main loop
    for (i in 1:Ngen) {
	if (verb) cat("[Gen = ",i,"] : ", sep="")

	if (verb) cat("Crossover...")                 # Begin crossover 
	Noff <- 0
	Lpar <- Lu[1:Npop]
	pcross <- exp(-beta*(Lpar - min(Lpar))/mean(Lpar - min(Lpar)))
	for (m in 1:Npop) {
	    coup <- sample(1:Npop, 2, prob = pcross)
	    offs <- crossover_tours(P[coup[1],], P[coup[2],], d, n)
	    for (j in 1:2) {
		repe <- FALSE
		off <- offs[j,]
		off <- gauge_tour(off, n)
		for (ki in 1:(K-1)) {
		    if (sum(abs(off - U[ki,])) == 0) {
			repe <- TRUE
			break
		    }
		}
		if (!repe) {
		    Noff <- Noff + 1
		    U[K,] <- off
		    L <- compute_tour_distance(off, d)
		    Lu[K] <- L
		    if (log) { Loff <- c(Loff, L) }
		    K <- K + 1
		}
		if (Noff == Npop) { break }
	    }
	    if (Noff == Npop) { break }
	}
	if (verb) cat("(", Noff, " offsprings) ",sep="") # End crossover

	if (verb) cat("Mutation...")                    # Begin mutation 
	Ncan <- Npop + Noff
	Nmut <- 0
	for (m in 1:Ncan) {
	    mut <- perturb_tour_4exc(U[m, ], C0, n)
	    mut <- gauge_tour(mut, n)
	    repe <- FALSE
	    for (k in 1:(K-1)) {
		if (sum(abs(mut - U[k,])) == 0) {
		    repe <- TRUE
		    break
		}
	    }
	    if (!repe) {
		Nmut <- Nmut + 1
		U[K,] <- mut
		L <- compute_tour_distance(mut, d)
		Lu[K] <- L
		if (log) { Lmut <- c(Lmut, L) }
		K <- K + 1
	    }
	}
	if (verb) cat("(", Nmut, " mutants) ",sep="")     # End mutation

	if (verb) cat("Local search...")           # Begin local search 
	Ncan <- Npop + Noff + Nmut
	Nvec <- 0
	for (m in 1:Ncan) {
	    if (runif(1) > local) { next } # Selecting local search candidates
	    vecL <- improve_tour_2opt(d, n, U[m, ])
	    vec <- gauge_tour(vecL$tour, n)
	    repe <- FALSE
	    for (k in 1:(K-1)) {
		if (sum(abs(vec - U[k,])) == 0) {
		    repe <- TRUE
		    break
		}
	    }
	    if (!repe) {
		Nvec <- Nvec + 1
		U[K,] <- vec
		Lu[K] <- vecL$distance
		if (log) { Lvec <- c(Lvec, vecL$distance) }
		K <- K + 1
	    }
	}
	if (verb) cat("(", Nvec, " 2-minima) ",sep="") # End local search

	if (verb) cat("Selection...")                   # Begin selection 
	Ncan <- Npop + Noff + Nmut + Nvec
	Lcan <- Lu[1:Ncan]
	Lcs <- sort(Lcan, index.return = TRUE)
	imin <- Lcs$ix[1]
	if (Lcan[imin] < Lmin) {
	    Lmin <- Lcan[imin]
	    Cmin <- U[imin,]
	    gmin <- i
	}
	P[1:elite,] <- U[Lcs$ix[1:elite],]
	psel <- exp(-beta*(Lcan - Lmin)/mean(Lcan))
	psel[Lcs$ix[1:elite]] <- 0
	isel <- sample(1:Ncan, Npop-elite, prob=psel)
	P[(elite+1):Npop,] <- U[isel,]
	if (verb) cat("done [Lmin = ", Lmin, "]\n",sep="") # End selection

	## Preparing data for next generation
	U <- U0
	U[1:Npop,] <- P
	Lu <- Lu0
	for (i in 1:Npop) {
	    L <- compute_tour_distance(P[i,], d)
	    Lu[i] <- L
	    if (log) { Lsel <- c(Lsel, L)}
	}
	K <- Npop + 1
    }
    if (log) {
	return(list(tour = Cmin, distance = Lmin,
		    generation = gmin, population = P, Lini = Lini,
		    Loff = Loff, Lmut = Lmut, Lvec = Lvec, Lsel = Lsel))
    } else {
	return(list(tour = Cmin, distance = Lmin,
		    generation = gmin, population = P))
    }
}



##' Gauging a tour for easy comparison.
##'
##' A tour of \eqn{n} vertices is a permutation of the ordered
##'     sequence 1:\eqn{n}, and it is represented as a vector
##'     containing the integers from 1 to \eqn{n} in the permuted
##'     sequence.  As a subgraph of the complete graph with \eqn{n}
##'     vertices, it is assumed that each vertex is adjacent with the
##'     anterior and posterior ones, with the first and last being
##'     also adjacent.
##'
##' With respect to the TSP, a tour is invariant under cyclic
##'     permutation and inversion, so that there exists \eqn{(n-1)!/2}
##'     different tours in a complete graph of \eqn{n} vertices.  When
##'     searching for tours it is common to find the same tour under a
##'     different representation.  Therefore, we need to establish
##'     wheter two tours are equivalent or not.  To this end, we can
##'     "gauge" the tour by permuting cyclically its elements until
##'     the first vertex is at position 1, and fix the orientation so
##'     that the second vertex is less than the last.  Two equivalent
##'     tours will have the same "gauged" representation.
##'
##' This function is used in [search_tour_genetic] to discard repeated
##'     tours which can be found during the execution of the
##'     algorithm.
##' 
##' @title Gauging a tour
##' @param To Tour to be gauged, a vector containing a permutation of
##'     the 1:n sequence
##' @param n Number of elements of the tour T
##' @return The gauged tour.
##' @author Cesar Asensio
##' @seealso [search_tour_genetic] implements a version of the genetic
##'     algorithm for the TSP.
##' @examples
##' set.seed(2)
##' T0 <- sample(1:9,9)   #        T0 = 2 6 5 9 7 4 1 8 3
##' gauge_tour(T0, 9)     # gauged T0 = 1 4 7 9 5 6 2 3 8
##' 
##' @encoding UTF-8
##' @md 
##' @export 
gauge_tour <- function(To, n) {
    p <- which(To==1)             # Position of vertex 1
    if (p==1) {
	Tg <- To
    } else {
	Tg <- c(To[p:n],To[1:(p-1)])  # First vertex is 1
    }
    if (Tg[2] < Tg[n]) {               # Orientation OK
	return(Tg)
    } else {                      # Reverse orientation
	return(c(1,rev(Tg[-1])))
    }
}


##' Crossover operation used by the TSP genetic algorithm.  It takes
##'     two tours and it computes two "offsprings" trying to exploit
##'     the structure of the cycles, see below.
##'
##' In the genetic algorithm, the crossover operation is a
##'     generalization of local search in which two tours are combined
##'     somehow to produce two tours, hopefully different from their
##'     parents and with better fitting function values.  Crossover
##'     widens the search while trying to keep the good peculiarities
##'     of the parents.  However, in practice crossover almost never
##'     lowers the fitting function when parents are near the optimum,
##'     but it helps to explore new routes.  Therefore, it is always a
##'     good idea to complement crossover with some deterministic
##'     local search procedure which can find another local optima;
##'     crossover also helps in evading local minima.
##'
##' In this routine, crossover is performed as follows.  Firstly, the
##'     edges of the parents are combined in a single graph, and the
##'     repeated edges are eliminated.  Then, the odd degree vertices
##'     of the resulting graph are matched looking for a low-weight
##'     perfect matching using a greedy algorithm.  Adding the
##'     matching to the previous graph yields an Eulerian graph, as in
##'     Christofides algorithm, whose final step leads to the first
##'     offspring tour.  The second tour is constructed by recording
##'     the second visit of each vertex by the Eulerian walk, and
##'     completing the resulting partial tour with the nearest
##'     neighbor heuristic.
##' 
##' @title Crossover operation used by the TSP genetic algorithm
##' @param C1 Vertex numeric vector of the first parent tour.
##' @param C2 Vertex numeric vector of the second parent tour.
##' @param d Distance matrix of the TSP instance.  It is used in the
##'     computation of the low-weight perfect matching.
##' @param n The number of vertices of the TSP complete graph.
##' @return A two-row matrix containing the two offsprings as vertex
##'     numeric vectors. 
##' @author Cesar Asensio
##' @seealso [search_tour_genetic] implements a version of the genetic
##'     algorithm for the TSP.
##' @examples
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' c1 <- sample(1:n)
##' c2 <- sample(1:n)
##' c12 <- crossover_tours(c1, c2, d, n)
##' compute_tour_distance(c1, d)        # 114.5848
##' compute_tour_distance(c2, d)        # 112.8995
##' compute_tour_distance(c12[1,], d)   # 116.3589
##' compute_tour_distance(c12[2,], d)   # 111.5184
##' 
##' @encoding UTF-8
##' @md 
##' @export 
crossover_tours <- function(C1, C2, d, n) {
    eg <- cbind(rbind(C1, c(C1[-1], C1[1])), rbind(C2, c(C2[-1], C2[1])))
    g <- make_graph(eg, n=n, dir=FALSE)
    g <- simplify(g)   # Simple graph resulting from C1 U C2

    ## Odd degree vertices in g
    dg  <- degree(g)
    Vimp <- (1:n)[dg == 3]
    Ni <- length(Vimp)

    ## Greedy low weight perfect matching in K(Vimp)
    col <- rep(0, n)
    col[Vimp] <- 1
    M <- rep(0,Ni)
    i <- 1
    for (v in Vimp) {
	if (col[v] == 1) {
	    col[v] <- 2
	    b <- d[v,]
	    co1 <- (col != 1)
	    nc1 <- sum(co1)
	    b[co1] <- rep(Inf,nc1)
	    u <- which.min(b)
	    col[u] <- 2
	    M[c(i,i+1)] <- c(v,u)
	    i <- i+2
	} else { next }
    }

    ## Adding the matching produces a 4-regular graph
    g <- g %>% add_edges(M)

    ## To extract the offspring tours we construct an eulerian walk:
    eug <- find_euler(g, 1)$walk[,1]
    q <- length(eug)

    ## The first tour is constructed from the walk using the visiting
    ## order.  The second is started by the vertices which are
    ## repeated by the walk:
    col <- h <- ph <- rep(0, n)
    col[1] <- h[1] <- 1
    j <- k <- 1
    for (i in 2:q) {
	v <- eug[i]
	if (col[v] == 0) {
	    j <- j+1
	    h[j] <- v
	    col[v] <- 1
	} else {
	    if (col[v] == 1) {
		ph[k] <- v
		k <- k+1
		col[v] <- 2
	    } else { next }
	}
    }

    ## We complete the second tour using the nearest-neighbor
    ## heuristic:
    if (k<=n) {
	v <- ph[k-1]
	for (j in k:n) {
	    inc <- which(col == 2)
	    b <- d[v, ]
	    b[inc] <- Inf
	    v <- which.min(b)
	    col[v] <- 2
	    ph[j] <- v
	}
    }

    rbind(h,ph)
}

##' Ant colony optimization (ACO) heuristic algorithm to search for a
##'     low-distance tour of a TSP instance.  ACO is a random
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
##'     tour in a TSP instance, leaving a pheromone trail on the edges
##'     of the tour.  At each step, each ant decides the next step
##'     based on the pheromone level and on the distance of each
##'     neighboring edge.  In a single iteration, each ant completes a
##'     tour, and the best tour is recorded.  Then the pheromone level
##'     of the edges of the best tour are enhanced, and the remaining
##'     pheromones evaporate.
##' 
##' Default parameter values have been chosen in order to find the
##'     optimum in the examples considered below.  However, it cannot
##'     be guarateed that this is the best choice for all cases.  Keep
##'     in mind that no polynomial time exact algorithm can exist for
##'     the TSP, and thus harder instances will require to fine-tune
##'     the parameters.  In any case, no guarantee of optimality of
##'     tours found by this method can be given, so they might be
##'     improved further by other methods.
##'
##' 
##' @title Ant colony optimization algorithm for the TSP
##' @param d Distance matrix of the TSP instance.
##' @param n Number of vertices of the complete TSP graph.  It can be
##'     less than the number of rows of the distance matrix d.
##' @param K Number of tour-searching ants.  Defaults to 200.
##' @param N Number of iterations.  Defaults to 50.
##' @param beta Inverse temperature which determines the thermal
##'     probability in selecting the next vertex in tour.  High beta
##'     (low temperature) rewards lower distances (and thus it gets
##'     stuck sooner in local minima), while low beta (high
##'     temperature) rewards longer tours, thus escaping from local
##'     minima.  Defaults to 3.
##' @param alpha Exponent enhancing the pheromone trail.  High alpha
##'     means a clearer trail, low alpha means more options.  It
##'     defaults to 5.
##' @param dt Basic pheromone enhancement at each iteration.  It
##'     defaults to 1.
##' @param rho Parameter in the (0,1) interval controlling pheromone
##'     evaporation rate.  Pheromones of the chosen tour increase in
##'     dt*rho, while excluded pheromones diminish in 1-rho.  A rho
##'     value near 1 means select just one tour, while lower values of
##'     rho spread the probability and more tours can be explored.  It
##'     defaults to 0.05.
##' @param log Boolean.  When TRUE, it also outputs two vectos
##'     recording the performance of the algorithm.  It defaults to FALSE.
##' @return A list with two components: $tour contains a permutation
##'     of the 1:n sequence representing the tour constructed by the
##'     algorithm, $distance contains the value of the distance
##'     covered by the tour.  When log=TRUE, the output list contains
##'     also the component $Lant, best tour distance found in the current
##'     iteration, and component $Lopt, best tour distance found
##'     before and including the current iteration.
##' @seealso [compute_distance_matrix] computes matrix distances using
##'     2d points, [improve_tour_2opt] improves a tour using
##'     2-exchanges, [plot_tour] draws a tour
##' @author Cesar Asensio
##' @examples 
##' ## Regular example with obvious solution (minimum distance 32)
##' m <- 6   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' set.seed(2)
##' b <- search_tour_ants(d, n, K = 70, N = 20)
##' b$distance    # Distance 32 (optimum)
##' plot_tour(z,b)
##'
##' ## Random points
##' set.seed(1)
##' n <- 15
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' b <- search_tour_ants(d, n, K = 50, N = 20)
##' b$distance    # Distance 32.48669 
##' plot_tour(z,b)
##'
##' @encoding UTF-8
##' @md 
##' @export 
search_tour_ants <- function(d, n, K=200, N=50, beta=3,
			     alpha=5, dt=1, rho=0.05, log=FALSE) {
    tran <- exp(-beta*d[1:n,1:n])
    diag(tran) <- rep(0,n)
    pher <- matrix(1,ncol=n,nrow=n)
    Tmin <- c()
    Emin <- c()
    if (log) { Lant <- Lopt <- c() }
    for (i in 1:N) {       # Repetimos las exploraciones
	tours <- c()
	costs <- c()
	for (j in 1:K) {   # Enviamos K hormigas a explorar
	    tour <- c(1)   # Todos empiezan desde 1
	    left <- c(NA,2:n)
	    p <- tran[1,]*pher[1,]^alpha
	    E <- 0
	    for (k in 2:(n-1)) { # Se construye el tour paso a paso
		nex <- sample(left[!is.na(left)],size=1,prob=p[!is.na(left)])
		tour <- c(tour,nex)
		E <- E + d[tour[k-1],nex]
		left[nex] <- NA
		p <- tran[nex,]*pher[nex,]^alpha
	    }
	    nex <- left[!is.na(left)]
	    tour <- c(tour,nex)
	    E <- E + d[tour[n-1],nex] + d[nex,1]
	    tours <- rbind(tours,tour)  # Retenemos el tour...
	    costs <- c(costs,E)         # ...y su coste...
	}
	if (log) { Lant <- c(Lant, costs) }
	min <- which.min(costs)         # ...y seleccionamos el mejor
	Emin <- c(Emin,costs[min])
	Tmin <- rbind(Tmin,tours[min,])
	pher <- (1-rho)*pher            # Evaporar feromonas
	for (l in 2:n) {                # Refuerzo de feromonas 
	    pher[tours[min,c(l-1,l)],tours[min,c(l-1,l)]] <-
		pher[tours[min,c(l-1,l)],tours[min,c(l-1,l)]] + rho*dt
	}
	pher[tours[min,c(1,n)],tours[min,c(1,n)]] <-
	    pher[tours[min,c(1,n)],tours[min,c(1,n)]] + rho*dt

    }
    best <- which.min(Emin)
    if (log) {
	return(list(tour = Tmin[best,], distance = Emin[best],
		    Lant = Lant, Lopt = Emin))        
    } else {
	return(list(tour = Tmin[best,], distance = Emin[best]))
    }
}

##' It computes the 1-tree lower bound for an optimum tour for a TSP instance.
##'
##' It computes the 1-tree lower bound for an optimum tour for a TSP
##' instance from vertex 1.  Internally, it creates the graph Kn-{v1}
##' and invokes [mst] from package [igraph] to compute the minimum
##' weight spanning tree.  If optional argument "degree" is TRUE, it
##' returns the degree seguence of this internal minimum spanning
##' tree, which is very convenient when embedding this routine in the
##' Held-Karp lower bound estimation routine.
##' 
##' @title Computing the 1-tree lower bound for a TSP instance
##' @param d Distance matrix.
##' @param n Number of vertices of the TSP complete graph.
##' @param degree Boolean: Should the routine return the degree
##'     seguence of the internal minimum spanning tree? Defaults to
##'     FALSE.
##' @return The 1-tree lower bound -- A scalar if the optional
##'     argument "degree" is FALSE.  Otherwise, a list with the
##'     previous 1-tree lower bound in the $bound component and the
##'     degree sequence of the internal minimum spanning tree in the
##'     $degree component.
##' @author Cesar Asensio
##' @seealso [improve_tour_2opt] tour improving using 2-opt,
##'     [improve_tour_3opt] tour improving using 3-opt,
##'     [compute_lower_bound_HK] for Held-Karp lower bound estimates.
##' @examples
##' m <- 10   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' b <- build_tour_2tree(d, n)
##' b$distance           # Distance 57.868
##' bi <- improve_tour_2opt(d, n, b$tour)
##' bi$distance          # Distance 48 (optimum)
##' compute_lower_bound_1tree(d,n)       # 45
##'
##' ## Random points
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' compute_lower_bound_1tree(d,n)       # 31.4477
##' bn <- build_tour_nn_best(d,n)
##' b3 <- improve_tour_3opt(d,n,bn$tour)
##' b3$distance                          # 35.081
##'
##' @encoding UTF-8
##' @md 
##' @export 
compute_lower_bound_1tree <- function(d, n, degree=FALSE) {
    v1 <- 1
    F1 <- d[v1,-1];   i1 <- which.min(F1);  d1 <- F1[i1]
    F2 <- F1[-i1];    i2 <- which.min(F2);  d2 <- F2[i2]
    Kn1 <- make_full_graph(n) %>% delete_vertices(v1)
    eKn1 <- as_edgelist(Kn1)
    q <- gsize(Kn1)
    dKn1 <- d[-v1,-v1]
    dl <- rep(0,q)
    for (e in 1:q) { dl[e] <- dKn1[eKn1[e,1],eKn1[e,2]] }
    Kn1 <- set_edge_attr(Kn1,"weight",value=dl)
    T1 <- mst(Kn1)         # Minimum spanning tree of Kn - v1
    eT1 <- as_edgelist(T1)
    dT1 <- 0
    for (e in 1:nrow(eT1)) { dT1 <- dT1 + dKn1[eT1[e,1],eT1[e,2]] }
    if (degree) {
	return(list(bound = dT1 + d1 + d2, degree = degree(T1)))
    } else {
	return(dT1 + d1 + d2)
    }
}

##' Held-Karp lower bound estimate for the minimum distance of an
##' optimum tour in the TSP
##'
##' An implementation of the Held-Karp iterative algorithm towards the
##'     minimum distance tour of a TSP instance.  As the algorithm
##'     converges slowly, only an estimate will be achieved.  The
##'     accuracy of the estimate depends on the stopping requirements
##'     through the number of iteration blocks "block" and the number
##'     of iterations per block "it", as well as the smallest allowed
##'     multiplier size "tsmall": When it is reached the algorithm
##'     stops.  It is also crucial to a good estimate the quality of
##'     the upper bound "U" obtained by other methods.
##'
##' The Held-Karp bound has the following uses: (1) assessing the
##' quality of solutions not known to be optimal; (2) giving an
##' optimality proof of a given solution; and (3) providing the
##' "bound" part in a branch-and-bound technique.
##'
##' Please note that recommended computation of the Held-Karp bound
##' uses Lagrangean relaxation on an integer programming formulation
##' of the TSP, whereas this routine uses the Cook algorithm to be
##' found in the reference below.
##'
##' 
##' @title Held-Karp lower bound estimate
##' @param d Distance matrix of the TSP instance.
##' @param n Number of vertices of the TSP complete graph.
##' @param U Upper bound of the minimum distance.
##' @param tsmall Stop criterion: Is the step size decreases beyond
##'     this point the algorithm stops.  Defaults to 0.001.
##' @param it Iterates inside each block.  Some experimentation is
##'     required to adjust this parameter: If it is large, the run
##'     time will be larger; if it is small, the accuracy will
##'     decrease.
##' @param block Number of blocks of "it" iterations.  In each block
##'     the size of the multiplier is halved.
##' @return An estimate of the Held-Karp lower bound -- A scalar.
##' @author Cesar Asensio
##' @references Cook et al.  \emph{Combinatorial Optimization} (1998)
##'     sec. 7.3.
##' @seealso [compute_distance_matrix] computes distance matrix of a
##'     set of points, [build_tour_nn_best] builds a tour using the
##'     best-nearest-neighbor heuristic, [improve_tour_2opt] improves
##'     a tour using 2-opt, [improve_tour_3opt] improves a tour using
##'     3-opt, [compute_lower_bound_1tree] computes the 1-tree lower
##'     bound.
##' @examples 
##' m <- 10   # Generate some points in the plane
##' z <- cbind(c(rep(0,m), rep(2,m), rep(5,m), rep(7,m)), rep(seq(0,m-1),4))
##' n <- nrow(z)
##' d <- compute_distance_matrix(z)
##' b <- build_tour_nn_best(d, n)
##' b$distance           # Distance 48.6055
##' bi <- improve_tour_2opt(d, n, b$tour)
##' bi$distance          # Distance 48 (optimum)
##' compute_lower_bound_HK(d,n,U=48.61)                    # 45.927
##' compute_lower_bound_HK(d,n,U=48.61,it=20,tsmall=1e-6)  # 45.791
##'
##' ## Random points
##' set.seed(1)
##' n <- 25
##' z <- cbind(runif(n,min=1,max=10),runif(n,min=1,max=10))
##' d <- compute_distance_matrix(z)
##' bn <- build_tour_nn_best(d,n)
##' b3 <- improve_tour_3opt(d,n,bn$tour)
##' b3$distance                                                    # 35.08155
##' compute_lower_bound_HK(d, n, U=35.1)                           # 34.80512
##' compute_lower_bound_HK(d, n, U=35.0816, it=20)                 # 35.02892
##' compute_lower_bound_HK(d, n, U=35.0816, tsmall = 1e-5)         # 34.81119
##' compute_lower_bound_HK(d, n, U=35.0816, it=50, tsmall = 1e-9)  # 35.06789
##'
##' @encoding UTF-8
##' @md 
##' @export 
compute_lower_bound_HK <- function(d, n, U, tsmall = 0.001,
				   it = 0.2*n, block = 100) {
    K1t <- compute_lower_bound_1tree(d, n)
    h <- rep(0, n)
    dh <- d
    K0 <- Kmax <- 0
    alpha <- 2
    beta <- 0.5
    for (b in 1:block) {
	for (k in 1:it) {
	    for (i in 1:(n-1)) {
		for (j in (i+1):n) {
		    dh[i,j] <- dh[j,i] <- d[i,j] - h[i] - h[j]
		}
	    }
	    Kgt <- compute_lower_bound_1tree(dh, n, degree=TRUE)
	    Kt <- Kgt$bound + 2*sum(h)
	    if (Kt > K0) {
		Kmax <- Kt
	    }
	    gt <- c(2, Kgt$degree)
	    if (sum(gt == 2) == n) {    # Tour found
		##return(list(bound = Kmax, degree = gt, h=h))
		return(Kmax)
	    }
	    t <- alpha*(U - Kt)/sum((2 - gt)^2)
	    h <- h + t*(2 - gt)
	    if (t < tsmall) {  # Too small changes
		##return(list(bound = Kmax, degree = gt, h=h))
		return(Kmax)
	    }
	    K0 <- Kt
	    ##cat("b = ", b, ", it = ", k, ", K = ", Kt, ", t = ", t, ", g(T1) = ", gt, "\n",sep="")
	}
	alpha <- alpha*beta
    }
    ##return(list(bound = Kmax, degree = gt, h=h))
    Kmax
}

##' This routine performs a version of the Branch-and-Bound algorithm
##'     for the Traveling Salesperson Problem (TSP).  It is an exact
##'     algorithm with exponential worst-case running time; therefore,
##'     it can be run only with a very small number of cities.
##'
##' The algorithm starts at city 1 (to avoid the cyclic permutation
##'     tour equivalence) and the "branch" phase consists on the
##'     decision of which city follows next.  In order to avoid the
##'     equivalence between a tour and its reverse, it only considers
##'     those tours for which the second city has a smaller vertex id
##'     that the last.  With \eqn{n} cities, the total number of tours
##'     explored in this way is \eqn{(n-1)!/2}, which clearly is
##'     infeasible unless \eqn{n} is small.  Hence the "bound" phase
##'     estimates a lower bound on the distance covered by the tours
##'     which already are partially constructed.  When this lower
##'     bound grows larger than an upper bound on the optimum supplied
##'     by the user or computed on the fly, the search stops in this
##'     branch and the algorithm proceeds to the next.  This
##'     complexity reduction does not help in the worst case, though.
##'
##' This routine represents the tree search by iterating over the
##'     sucessors of the present tree vertex and calling itself when
##'     descending one level.  The leaves of the three are the actual
##'     tours, and the algorithm only reaches those tours whose cost
##'     is less than the upper bound provided.  By default, the
##'     algorithm will plot the tour found if the coordinates of the
##'     cities are supplied in the "z" input argument.
##'
##' When the routine takes too much time to complete, interrupting the
##'     run would result in losing the best tour found.  To prevent
##'     this, the routine can store the best tour found so far so that
##'     the user can stop the run afterwards.
##' 
##' @title Branch-and-Bound algorithm for the TSP
##' @param d Distance matrix of the TSP instance
##' @param n Number of vertices of the TSP complete graph
##' @param verb If detailed operation of the algorithm should be
##'     echoed to the console.  It defaults to FALSE
##' @param plot If tours found by the algorithm should be plotted
##'     using [plot_tour].  It defaults to TRUE
##' @param z Points to plot the tours found by the algorithm.  It
##'     defaults to NA; it should be set if plot is TRUE or else
##'     [plot_tour] will not plot the tours
##' @param tour Best tour found by the algorithm.  If the algorithm
##'     has ended its complete run, this is the optimum of the TSP
##'     instance.  This variable is used to store the internal state
##'     of the algorithm and it should not be set by the user
##' @param distance Distance covered by the best tour found.  This
##'     variable is used to store the internal state of the algorithm
##'     and it should not be set by the user
##' @param upper Upper bound on the distance covered by the optimum
##'     tour.  It can be provided by the user or the routine will use
##'     the result found by the heuristic [build_tour_nn_best]
##' @param col Vectors of "colors" of vertices.  This variable is used
##'     to store the internal state of the algorithm and it should not
##'     be set by the user
##' @param last Last vertex added to the tour being built by the
##'     algorithm.  This variable is used to store the internal state
##'     of the algorithm and it should not be set by the user
##' @param partial Partial tour built by the algorithm.  This variable
##'     is used to store the internal state of the algorithm and it
##'     should not be set by the user
##' @param covered Partial distance covered by the partial tour built
##'     by the algorithm.  This variable is used to store the internal
##'     state of the algorithm and it should not be set by the user
##' @param call Number of calls that the algorithm performs on itself.
##'     This variable is used to store the internal state of the
##'     algorithm and it should not be set by the user
##' @param save.best.result The time needed for a complete run of this
##'     algorithm may be exponentially large.  Since it only will
##'     return its results if it ends properly, we can save to a file
##'     the best result found by the routine at a given time when
##'     save.best.result = TRUE (default is FALSE).  Then, the user
##'     will be allowed to stop the run of the algorithm without
##'     losing the (possibly suboptimal) result.
##' @param filename.best.result The name of the file used to store the
##'     best result found so far when save.best.result = TRUE.  It
##'     defaults to "best_result_find_tour_BB.Rdata".  When loaded,
##'     this file will define the best tour in variable "Cbest".
##' @param order Numeric vector giving the order in which vertices
##'     will be search by the algorithm.  It defaults to NA, in which
##'     case the algorithm will take the order of the tour found by
##'     the heuristic [build_tour_nn_best].  If the user knows in
##'     advance some good tour and he/she wishes to use the order of
##'     its vertices, it should be taken into account that the third
##'     vertex used by the algorithm is the last vertex of the tour!
##' @return A list with nine components: $tour contains a permutation
##'     of the 1:n sequence representing the best tour constructed by
##'     the algorithm, $distance contains the value of the distance
##'     covered by the tour, which if the algorithm has ended properly
##'     will be the optimum distance.  Component $call is the number
##'     of calls the algorithm did on itself.  The remaining
##'     components are used to transfer the state of the algorithm in
##'     the search three from one call to the next; they are $upper
##'     for the current upper bound on the distance covered by the
##'     optimum tour, $col for the "vertex colors" used to mark the
##'     vertices added to the partially constructed tour, which is
##'     stored in $partial.  The distance covered by this partial tour
##'     is stored in $covered, the last vertex added to the partial
##'     tour is stored in $last, and the "save.best.result" and
##'     "filename.best.result" input arguments are stored in
##'     $save.best.result and $filename.best.result.
##' @author Cesar Asensio
##' @examples
##' ## Random points
##' set.seed(1)
##' n <- 10
##' z <- cbind(runif(n, min=1, max=10), runif(n, min=1, max=10))
##' d <- compute_distance_matrix(z)
##' bb <- find_tour_BB(d, n)
##' bb$distance                        # Optimum 26.05881
##' plot_tour(z,bb)
##' ## Saving tour to a file (useful when the run takes too long):
##' ## Can be stopped after tour is found
##' ## File is stored in tempdir(), variable is called "Cbest"
##' find_tour_BB(d, n, save.best.result = TRUE, z = z)
##' 
##' @encoding UTF-8
##' @md 
##' @export 
find_tour_BB <- function(d, n, verb = FALSE, plot = TRUE, z = NA,
			 tour     = rep(0,n),
			 distance = 0,
			 upper    = Inf,
			 col      = c(1, rep(0, n-1)),
			 last     = 1,
			 partial  = c(1, rep(NA, n-1)),
			 covered  = 0,
			 call     = 0,
			 save.best.result = FALSE, filename.best.result = "best_result_find_tour_BB.Rdata",
			 order = NA) {

    if (call == 0) {
	cat("Running branch-and-bound on a TSP with n =", n,
	    "vertices...\n")
    }

    C <- partial
    D <- covered
    U <- upper
    M <- col
    w <- last
    Cmin <- tour
    Dmin <- distance
    call <- call + 1 # Call count
    J <- call
    level <- sum(M)

    if (verb) {
	cat("///\n///   Entering call ", J, " at level ", level,
	    "\n///\n", sep="")
    } else {
	cat("Call ", J, "...\n", sep="")
    }

    if (is.infinite(U)) {
	bU <- build_tour_nn_best(d, n)
	Cmin <- bU$tour
	Cmin <- gauge_tour(Cmin, n)
	Dmin <- U <- bU$distance
    }

    if (is.na(order[1])) { order <- c(Cmin[1:2], Cmin[n], Cmin[3:(n-1)]) }

    dw <- d[w,]
    inc0 <- which(M == 1)
    ## sc <- sort(dw, index.return=TRUE, decreasing=TRUE) # Farthest first
    ## Vleft <- setdiff(sc$ix, inc0) # setdiff keeps order of sort(dw...)
    Vleft <- setdiff(order, inc0) # setdiff keeps order of "order"

    if (verb) {
	cat("(Call ", J, ") Vertices to explore from level ", level,
	    ": ", paste0(Vleft, collapse=" "),"\n",sep="")
    }

    for (v in Vleft) { ## BRANCH: We add a vertex and see what happens

	if (verb) {
	    cat("(Call ", J, ") Exploring from vertex ", w,
		" (at level ", level, ") to ", v, " (at level ", level + 1,
		")...", sep="")
	}

	## Compute the distance covered by adding v to the cycle
	if (level == 1) {
	    if (v == n) {
		if (verb) cat("ungauged\n")
		next
	    }
	    C[2] <- v
	    de <- d[1,v]
	}
	if (level == 2) {
	    if (v > C[2]) {
		C[n] <- v
		de <- d[1,v]
	    } else {
		if (verb) cat("ungauged\n")
		next
	    }
	}
	if (verb) cat("\n")
	if (level == 3) {
	    C[3] <- v
	    de <- d[C[2],C[3]]
	}
	if (level > 3) {
	    C[level] <- v
	    de <- d[w,v]
	}
	D <- D + de
	M[v] <- 1

	## Last step: Close the tour!
	if (level == n-1) {
	    D <- D + d[v, C[n]]
	    if (D < U) {
		U <- D      # Decrease the upper bound
		Cmin <- C   # Update the optimum candidate
		Dmin <- D
	    }
	    Lv <- m1 <- m2 <- 0    # Residual contributions vanish
	    if (verb) {
		cat("(Call ", J, ")",
		    "       Upper bound  U = ",  U, "\n",
		    rep(" ", 8), "          Distance  D = ",  D, "\n",
		    rep(" ", 8), "        Tour found  C = ",
		    paste0(replace(C, which(is.na(C)), "-"), collapse=" "),
		    "\n",sep="")
	    } else {
		cat("...tour found: Distance", Dmin, "\n")
	    }

	    if (save.best.result) {
		assign("Cbest", value = list(tour = Cmin, distance = Dmin))
		dirsep <- if(Sys.info()["sysname"] == "Windows") { "\\" } else { "/" }
		tfile <- paste(tempdir(), filename.best.result, sep = dirsep)
		save(Cbest, file = tfile)
		cat("\nSaved best result found in variable \"Cbest\" in file ", tfile, "\n", sep = "")
	    }

	    if (plot) {
		if (is.matrix(z)) {
		    plot_tour(z = z, h = list(tour = Cmin, distance = Dmin))
		}
	    }
	} else {

	    ## If the tour is not closed, we compute a lower BOUND:  
	    if (level == 1) { ex1 <-    1 } else { ex1 <- C[n] }
	    if (level <= 2) { ex2 <- C[2] } else { ex2 <- C[level] }

	    if (level == n-2) {
		u <- Vleft[Vleft != v]
		m1 <- d[ex1, u]         # Last edges
		m2 <- d[ex2, u]         #   "   "
		Lv <- 0                 # No residual edges
	    }

	    if (level == n-3) {
		pq <- Vleft[Vleft != v]
		Lv  <- d[pq[1], pq[2]]
		m1 <- min(d[ex1, pq[1]], d[ex1, pq[2]])
		m2 <- min(d[ex2, pq[1]], d[ex2, pq[2]])
	    }
	    if (level < n-3) {
		inc1 <- c(inc0, v)
		m1 <- min(d[ex1, -inc1])
		m2 <- min(d[ex2, -inc1])
		dHK <- d[-inc1, -inc1]
		Lv <- compute_lower_bound_HK(dHK, n - (level + 1),
					     U - m1 - m2 - D, it=20)
		Lv <- Lv - max(dHK)
	    }

	    Lv <- Lv + m1 + m2 + D

	    if (verb) {
		cat("(Call ", J, ")",
		    " Local lower bound Lv = ", Lv, "\n", rep(" ", 8),
		    "       Upper bound  U = ",  U, "\n", rep(" ", 8),
		    "  Partial distance  D = ",  D, "\n", rep(" ", 8),
		    "      Partial tour  C = ",
		    paste0(replace(C, which(is.na(C)), "-"), collapse=" "),
		    "\n", sep="")
	    }

	    ## We decide whether we descend the tree:
	    if (Lv > U) {
		if (verb) {
		    cat("(Call ", J, ") Still at level ", level,
			", next vertex...\n", sep="")
		}
	    } else {
		if (verb) {
		    cat("(Call ", J, ") Descending to level ",
			level + 1, "...\n",sep="")
		}

		S <- find_tour_BB(d, n, verb = verb, plot=plot, z=z,
				  tour     = Cmin,
				  distance = Dmin,
				  upper   = U, col     = M, last = v,
				  partial = C, covered = D, call = J,
				  save.best.result = save.best.result,
				  filename.best.result = filename.best.result,
				  order = order)
		Cmin <- S$tour
		Dmin <- S$distance
		U <- S$upper
		J <- S$call
		save.best.result <- S$save.best.result
		filename.best.result <- S$filename.best.result
	    }
	}

	## We restore the initial state for the next iteration:
	M[v] <- 0
	D <- D - de
	if (level == 1) { C[2] <- NA }
	if (level == 2) { C[n] <- NA }
	if (level > 2)  { C[level] <- NA }

    }

    ## At this point, we backtrack (if level > 1) or we are done!
    if (level == 1) {
	if (verb) {
	    cat("(Call ", J, ") End of run!\n", sep="")
	} else {
	    cat("Call ", J, ": End of run\n", sep="")
	}
    } else {
	if (verb) {
	    cat("(Call ", J, ") Returning from level ", level,
		"!\n", sep="")
	}
    }

    S <- list(tour = Cmin, distance = Dmin, upper   = U, col  = M,
	      last =    w, partial  =    C, covered = D, call = J,
	      save.best.result = save.best.result,
	      filename.best.result = filename.best.result)
    return(S)
}
