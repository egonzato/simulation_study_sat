rtruncexp=function(n, rate, max) {
  # generalise to uniform
  if (rate == 0) {
    return(runif(n, 0, max))
  }
  # define u from uniform
  u=runif(n)
  # calculate truncated exponential
  -log(1 - u * (1 - exp(-rate * max)))/rate
}