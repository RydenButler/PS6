sg.int <- function(g,...,lower,upper){
  require("SparseGrid")
  # Set lower bound over which to integrate (must be length of dimensions)
  lower <- floor(lower)
  # Set upper bound over which to integrate (must be length of dimensions)
  upper <- ceiling(upper)
  # Throw error if the lower bound is greater than the upper bound
  if (any(lower > upper)) stop("lower must be smaller than upper")
  # Create matrix of all combinations of supplied vectors
  gridss <- as.matrix(expand.grid(seq(lower[1], upper[1]-1, by=1),
                                  seq(lower[2], upper[2]-1, by=1)))
  # Creates nodes and weights for integration; 
  # KPU method is nested rule for unweighted integral over [0,1]
  # k is the accuracy level of integration
  sp.grid <- createIntegrationGrid( 'KPU', dimension = 2, k = 5)
  # Create node values for lower bounds of integration
  nodes <- gridss[1, ] + sp.grid$nodes
  # Store node weight values for lower bounds of interation
  weights <- sp.grid$weights
  # Create node values for all other points of integration
  for (i in 2:nrow(gridss)) {
    # Create node values for combination of points i
    nodes <- rbind(nodes, gridss[i, ] + sp.grid$nodes)  
    # Store node weight values for combination i
    weights <- c(weights, sp.grid$weights)
  }
  # Evaluate function over each set of nodes
  gx.sp <- apply(nodes, 1, g,...)
  # Weight values and sum rows
  val.sp <- gx.sp %*%weights
  # Return values
  val.sp
}

sg.int(sin, lower= c(-1,-1), upper = c(1,1))


