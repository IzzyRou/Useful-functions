# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

####COLLECTION OF USEFUL FUNCTIONS#########

##STANDARD ERROR###
se <- function(x) sd(x)/sqrt(length(x))

##95%CI, t distribution or norm

CI<-function(n,s,a, type="normal"){
  if(type=="normal"){
    error <- qt(0.975,df=n-1)*s/sqrt(n)

  }else{error <-qnorm(0.975)*s/sqrt(n)}

  left <- a-error
  right <- a+error
  return(c(left,right))
}

##95%CI bootstrapping mean
bootstrap<-function(x,n,nboot=100000){
  # Data for the example 6
  # sample mean
  xbar = mean(x)
  # Generate 20 bootstrap samples, i.e. an n x 20 array of
  # random resamples from x
  tmpdata = sample(x,n*nboot, replace=TRUE)
  bootstrapsample = matrix(tmpdata, nrow=n, ncol=nboot)
  # Compute the means x
  bsmeans = colMeans(bootstrapsample)
  # Compute ????? for each bootstrap sample
  deltastar = bsmeans - xbar
  # Find the 0.1 and 0.9 quantile for deltastar

  d = quantile(deltastar, c(0.025, 0.975))

  # Calculate the 80% confidence interval for the mean.

  ci = xbar - c(d[2], d[1])

  return(cat("Confidence interval: ",ci, "\n"))
}
