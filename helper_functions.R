
rank_by_distance2 <- function(frame, metric='mode'){
  
  # Define
  parameter.draws <- frame
  
  # Determine metric to define 'center'
  if (metric == 'mean'){
    # Mean of parameters
    pt <- apply(parameter.draws, 2, mean)
    
  } else if (metric == 'median'){
    # Median of parameters
    pt <- apply(parameter.draws, 2, median)
    
  } else if (metric == 'mode'){
    # Mode of parameters
    pt <- apply(parameter.draws, 2, mlv, method = "lientz", bw = 0.2)
    
  } else {
    # Throw error
    stop("Metric must be 'mean', 'median', or 'mode'.")
    
  }
  
  # vector of weights 
  weights = 1/apply(parameter.draws, 2, var)
  # Vector of distances
  distances <- apply(parameter.draws, 1, function(x) weights%*%(x - pt) ^ 2)
  
  # Append to data and sort
  parameter.draws$distance <- distances
  
  # Rank and re-index
  parameter.draws.ranked <- parameter.draws[order(parameter.draws$distance), ]
  row.names(parameter.draws.ranked) <- NULL
  
  # Return sorted results and the point used ranking
  return(list(parameter.draws.ranked, pt))
}

################################################################################

### Helper functions ###

# General Euler ODE solvers
euler <- function (t, dt, fun, ic)
{ 
  p <- length(ic)
  n <- t/dt 
  xmat <- matrix(0,ncol=p,nrow=n)
  x <- ic 
  xmat[1,] <- x
  for (i in 2:n) { 
    x <- x + fun(x)*dt
    xmat[i,] <- x
  } 
  ts(xmat,start=0,deltat=dt) 
} 

# SIRT model
ode.sirt <- function(x, b=beta, g=gamma, d=delta)
{
  c(
    -b * x[1] * x[2],
    b * x[1] * x[2] - g * x[2],
    g * x[2] - d * x[3],
    d * x[3]
  )
}

ode.sirt2 <- function(x, b=beta, g=gamma, d=delta, r=rho, e=eps,p=psi)
{
  c(
    g*(1+r+e+p-exp(-b*(x[1]-e))-x[1]),
    d * x[1]
  )
}


################################################################################
