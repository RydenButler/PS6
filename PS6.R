library(testthat)
library(microbenchmark)
library(cubature)


# The original, two-dimensional function with comments
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

# A k-dimensional version of the integration function; includes minor additions
sg.mult <- function(g,...,lower,upper, dim){
  require("SparseGrid")
  # Check valid dim
  if(dim%%1 != 0) stop('dim may not be a decimal value')
  # Check valid starting values
  if(!is.vector(lower)) stop('lower must be a vector of length = dim')
  if(!is.vector(upper)) stop('upper must be a vector of length = dim')
  # Check valid starting value lengths
  if (length(lower) != dim) stop('lower must contain values for each dimension')
  if (length(upper) != dim) stop('upper must contain values for each dimension')
  # Set lower bound over which to integrate (must be length of dimensions)
  lower <- floor(lower)
  # Set upper bound over which to integrate (must be length of dimensions)
  upper <- ceiling(upper)
  # Throw error if the lower bound is greater than the upper bound
  if (any(lower > upper)) stop("lower must be smaller than upper")
  # Create list of all integration bounds
  all.dims <- lapply(1:dim, function(x) {seq(lower[x], upper[x]-1, by = 1)})
  # Create matrix of all combinations of supplied vectors
  gridss <- as.matrix(expand.grid(all.dims))
  # Creates nodes and weights for integration; 
  # KPU method is nested rule for unweighted integral over [0,1]
  # k is the accuracy level of integration
  sp.grid <- createIntegrationGrid( 'KPU', dimension = dim, k = 5)
  # Create node values for lower bounds of integration
  nodes <- gridss[1, ] + sp.grid$nodes
  # Store node weight values for lower bounds of interation
  weights <- rep(sp.grid$weights, nrow(gridss))
  # Create node values for all other points of integration
  nodes <- lapply(1:nrow(gridss), function(x) gridss[x, ] + sp.grid$nodes)
  nodes <- do.call(rbind, nodes)
  # Evaluate function over each set of nodes
  gx.sp <- apply(nodes, 1, g,...)
  # Weight values and sum rows
  val.sp <- gx.sp %*%weights
  # Return values
  return(as.numeric(val.sp))
}

# A parallelized version of the k-dimensional integration functions
sg.parallel <- function(g,...,lower,upper, dim){
  require("SparseGrid")
  require('parallel')
  no_cores <- detectCores()
  cl <- makeCluster(no_cores, type = 'FORK')
  # Check valid dim
  if(dim%%1 != 0) stop('dim may not be a decimal value')
  # Check valid starting values
  if(!is.vector(lower)) stop('lower must be a vector of length = dim')
  if(!is.vector(upper)) stop('upper must be a vector of length = dim')
  # Check valid starting value lengths
  if (length(lower) != dim) stop('lower must contain values for each dimension')
  if (length(upper) != dim) stop('upper must contain values for each dimension')
  # Set lower bound over which to integrate (must be length of dimensions)
  lower <- floor(lower)
  # Set upper bound over which to integrate (must be length of dimensions)
  upper <- ceiling(upper)
  # Throw error if the lower bound is greater than the upper bound
  if (any(lower > upper)) stop("lower must be smaller than upper")
  # Create list of all integration bounds
  all.dims <- parLapply(cl, 1:dim, 
                        function(x) {seq(lower[x], upper[x]-1, by = 1)})
  # Create matrix of all combinations of supplied vectors
  gridss <- as.matrix(expand.grid(all.dims))
  # Creates nodes and weights for integration; 
  # KPU method is nested rule for unweighted integral over [0,1]
  # k is the accuracy level of integration
  sp.grid <- createIntegrationGrid( 'KPU', dimension = dim, k = 5)
  # Store node weight values for lower bounds of interation
  weights <- rep(sp.grid$weights, nrow(gridss))
  # Create node values for all other points of integration
  nodes <- parLapply(cl, 1:nrow(gridss), 
                     function(x) gridss[x, ] + sp.grid$nodes)
  nodes <- do.call(rbind, nodes)
  # Evaluate function over each set of nodes
  gx.sp <- parApply(cl = cl, X = nodes, MARGIN = 1, FUN = g, ...)
  # Weight values and sum rows
  val.sp <- gx.sp %*%weights
  # Stop cluster to prevent errors
  stopCluster(cl)
  # Return values
  return(as.numeric(val.sp))
}


### Create function for testing ###
# 2-D function
fn2 <- function(x) x[1] + x[2]^2
# 3-D function
fn3 <- function(x) x[1] + x[2]^2 + x[3]^3
# Integral from -1 to 1
#fn3dx (x[1]^2)/2 + x[1]x[2]^2 + x[1]x[3]^3 
#-1 -> 1| .5 + x[2]^2 + x[3]^3 - .5 + x[2]^2 + x[3]^3
#fn3dy 2/3*(x[2]^3) + 2x[2]x[3]^3 
#-1 -> 1| 2/3 + 2x[3]^3 + 2/3 + 2x[3]^3
#fn3dz 4/3*x[3] + x[3]^4 | 4/3 + 1 + 4/3 - 1 = 8/3
# 5-D function
fn5 <- function(x) x[1] + x[2]^2 + x[3]^3 + x[4]^4 + x[5]^5
# Integral from -1 to 1
#fn5dx: (x[1]^2)/2 + x[1]x[2]^2 + x[1]x[3]^3  + x[1]x[4]^4 + x[1]x[5]^5
#-1 -> 1| 1/2 + x[2]^2 + x[3]^3 + x[4]^4 + x[5]^5 - 1/2 + x[2]^2 + x[3]^3 + x[4]^4 + x[5]^4
#fn5dy: 2/3*(x[2]^3) + 2x[2]x[3]^3 + 2x[2]x[4]^4 + 2x[2]x[5]^5
#-1 -> 1|| 2/3 + 2x[3]^3 2x[4]^4 + 2x[5]^5 + 2/3 + 2x[3]^3 2x[4]^4 + 2x[5]^5
#fn5dz: 4/3*x[3] + x[3]^4 + 4x[3]x[4]^4 + 4x[3]x[5]^5 
#-1 -> 1|| 4/3 + 1 + 4x[4]^4 + 4x[5]^5 + 4/3 - 1 + 4x[4]^4 + 4x[5]^5
#fn5dzz: 8/3x[4] + 8/5*x[4]^5 + 8x[4]x[5]^5 
#-1 -> 1|8/3 + 8/5 + 8x[5]^5 + 8/3 + 8/5 + 8x[5]^5
#fndzzz: 128/15x[5] + 16/6*x[5]^6
#-1 -> 1| 128/15 + 16/6 + 128/15 - 16/6 = 256/15

# Run functions for first pass of checks
sg.int(fn2, lower= c(-1,-1), upper = c(1,1))
sg.mult(fn3, lower = rep(-1, 3), upper = rep(1, 3), dim = 3)
sg.mult(fn5, lower = rep(-1, 5), upper = rep(1, 5), dim = 5)
sg.parallel(fn3, lower = rep(-1, 3), upper = rep(1, 3), dim = 3)
sg.parallel(fn5, lower = rep(-1, 5), upper = rep(1, 5), dim = 5)

#### Conduct formal tests of equality ###

# Check non-parallel in 3 dimensions
expect_equal(sg.mult(fn3, lower = rep(-1, 3), upper = rep(1, 3), dim = 3),
             8/3, tolerance = 0.00001)
# Check parallel in 5 dimension
expect_equal(sg.parallel(fn5, lower = rep(-1, 5), upper = rep(1, 5), dim = 5),
             256/15, 
             tolerance = 0.00001)
### Conduct formal tests for validation checks ###

# Check valid dim
expect_error(sg.mult(fn3, lower = rep(-1,3), 
                     upper = rep(1, 3), dim = 3.5),
             'dim may not be a decimal value')
# Check valid lower with non-parallel
expect_error(sg.mult(fn3, lower = matrix(rep(-1, 3), ncol = 3), 
                     upper = rep(1, 3), dim = 3),
             'lower must be a vector of length = dim')
# Check valid upper with parallel
expect_error(sg.parallel(fn3, lower = rep(-1, 3), 
                         upper = matrix(rep(1, 3), ncol = 3), dim = 3),
             'upper must be a vector of length = dim')
# Check upper and lower lengths
expect_error(sg.parallel(fn3, lower = rep(-1, 3) , upper = c(1,1), dim = 3),
             'upper must contain values for each dimension')
expect_error(sg.parallel(fn3, lower = c(-1,-1) , upper = rep(1, 3), dim = 3),
             'lower must contain values for each dimension')
# Check lower & upper ordering 
expect_error(sg.mult(fn3, lower = c(-1,-1,1), upper = c(1,1,-1), dim = 3),
             'lower must be smaller than upper')

# Microbenchmark times
microbenchmark(sg.mult(fn3, lower = rep(-1, 3), 
                       upper = rep(1, 3), dim = 3),
               sg.parallel(fn3, lower = rep(-1, 3), 
                           upper = rep(1, 3), dim = 3))
microbenchmark(sg.mult(fn5, lower = rep(-1, 5), 
                       upper = rep(1, 3), dim = 5),
               sg.parallel(fn5, lower = rep(-1, 5), 
                           upper = rep(1, 5), dim = 5))
# The parallel function is longer in both instances.
# This may be a result of the simplicity of the calculations, causing the 
#      allocation to multiple nodes to take longer than the actual integration
# Create 10-D function to test if the time improves over dimensionality

#X <- paste0('x[', 1:10, ']^', 1:10, collapse = ' + ') 
fn10 <- function(x) (x[1]^1 + x[2]^2 + x[3]^3 + x[4]^4 + x[5]^5 + x[6]^6 + 
                       x[7]^7 + x[8]^8 + x[9]^9 + x[10]^10)

microbenchmark(sg.mult(fn10, lower = rep(-1, 10), 
                       upper = rep(1, 10), dim = 10),
               sg.parallel(fn10, lower = rep(-1, 10), 
                           upper = rep(1, 10), dim = 10))



### Integrate functions with adaptIntegrate ###
adaptIntegrate(fn3, c(-1,-1,-1), c(1,1,1))$integral
adaptIntegrate(fn5, c(-1,-1,-1,-1,-1), c(1,1,1,1,1))$integral

### Find maximum of functions ###
optim()
optimize()

