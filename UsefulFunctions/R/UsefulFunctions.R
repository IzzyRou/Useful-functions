
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

####COLLECTION OF USEFUL FUNCTIONS#########

##STANDARD ERROR
se <- function(x) sd(x)/sqrt(length(x))

##95%CI, t distribution or norm from http://www.cyclismo.org/tutorial/R/confidence.html
#' normCI
#'
#' @param n number in popularion
#' @param s standard deviation
#' @param a mean
#' @param type "normal" or "t" at moment
#'
#' @return
#' @export
#'
#' @examples
CI<-function(n,s,a, type="normal"){
  if(type=="normal"){
    error <- qnorm(0.975)*s/sqrt(n)

  }else{error <-qt(0.975,df=n-1)*s/sqrt(n)}

  left <- a-error
  right <- a+error
  return(c(left,right))
}

######
#' bootstrapMean
#' 95%CI bootstrapping mean
#'
#' @param x data
#' @param n size of sample
#' @param nboot number bootstrap samples
#'
#' @return
#' @export
#'
#' @examples
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

### incidence per group per day
#' Title
#'
#' @param obs.date date case observed
#' @param obs.group name of group observed e.g. age class, sex
#'
#' @return
#' @export
#'
#' @examples
daily.incidence.group<-function(obs.date,obs.group){
  group<-obs.group
  date.range<-c(min(obs.date),max(obs.date))
  tt<-table(obs.date,obs.group)
  Groups<-colnames(tt)
  Dates<-(rownames(tt))
  date<-seq(date.range[1],date.range[2],1)
  f<-which(as.character(date) %in% Dates)
  incidence<-data.frame(date)
  for (i in 1:ncol(tt)){
    x<-rep(0,length(date))
    x[f]<-tt[,i]
    incidence<-cbind(incidence,x)
    names(incidence)[i+1]<-Groups[i]
  }
  return (incidence)
}

######
#' dailyFOI
#'Force of Infection
#' @param incidence timeseries incidence
#' @param n number
#' @param t time
#' @param SI serial interval distribution
#'
#' @return
#' @export
#'
#' @examples
daily.foi.group<-function(incidence,n,t,SI){
  foi<-incidence
  foi[,2:n]<-0
  WS<-rev(SI)%*%matrix(1,1,n)
  if (n>1){
    for (i in 1:t){
      f=max(c(1,(i-SItrunc)))
      foi[i,2:n]=colSums(incidence[f:i,2:n]*WS[((SItrunc+1)-(i-f)):(SItrunc+1),])
    }
  }else{
    for (i in 1:t){
      f=max(c(1,(i-SItrunc)))
      foi[i,2:n]=sum(incidence[f:i,2:n]*WS[((SItrunc+1)-(i-f)):(SItrunc+1),])
    }
  }
  return (foi)
}

####
#' discreteSI
#'creates discrete serial interval
#' @param mu
#' @param CV
#' @param SItrunc
#'
#' @return
#' @export
#'
#' @examples
DiscretizeSI <- function(mu,CV,SItrunc){
  SIdistr <- sapply(0:SItrunc, function(k) EpiEstim::DiscrSI(k, mu, CV*mu))
  SIdistr <- SIdistr/sum(SIdistr)
  return(SIdistr)
}

######Raleigh likelihood, hazard and survival functions ####

Hral<-function(a,t) a*t
Sral<-function(a,t) exp(-a*0.5*t^2)
Lral<-function(t,a) a*t*exp(-a*0.5*t^2)




#############
#' legendCol
#' continuous colour scale for base plot
#'
#' @param col
#' @param lev
#'
#' @return
#' @export
#'
#' @examples

#'plot(1:100,1:100, type="p", col = colfunc(100), pch=19)
#'colfunc <- colorRampPalette(c("snow1","snow2","snow3","seagreen","orange","firebrick","red"), space = "rgb",bias=1)#colors
#'colfunc2 <- colorRampPalette(c('yellow','orange','red','brown'),bias=1)
#'legend.col(colfunc(1000),c(1:100))
#'legend.col(colfunc2(1000),c(1:100))

legend.col <- function(col, lev){

  opar <- par

  n <- length(col)

  bx <- par("usr")

  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n

  xx <- rep(box.cx, each = 2)

  par(xpd = TRUE)
  for(i in 1:n){

    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])

  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}



####bob functions####

# -----------------------------------
#' safeRead
#'
#' Reads in text file (tab-delimited by default) and replaces bad characters with replacements.
#' \cr\cr
#' WARNING - can overwrite original file with new values if \code{overwrite==TRUE}, although this is not the default.
#'
#' @export

safeRead <- function(fileName, delim="\t", useHeader=TRUE, report=TRUE, badCharacters=c("-1.#IND","1.#INF","-999"), replacements=c(0,0,0), overwrite=FALSE) {

  # read in raw data. All values are read in as characters at this stage, and headers are ignored
  data <- read.delim(fileName, sep=delim, header=F, colClasses="character")

  # replace bad characters
  badCount <- 0
  for (j in 1:length(badCharacters)) {
    badCount <- badCount + sum(data==badCharacters[j])
    data[data==badCharacters[j]] <- replacements[j]
  }

  # report to console
  if (report & badCount>0) {
    cat(paste(fileName,": ",badCount," bad characters replaced\n",sep=""))
  }

  # now convert variables to sensible types
  if (useHeader) {
    df <- data.frame(2:nrow(data))
    for (i in 1:ncol(data)) {
      df <- cbind(df, type.convert(data[-1,i]))
    }
    df <- df[,-1]
    names(df) <- as.character(data[1,])
  } else {
    df <- data.frame(1:nrow(data))
    for (i in 1:ncol(data)) {
      df <- cbind(df, type.convert(data[,i]))
    }
    df <- df[,-1]
    names(df) <- paste("X",1:ncol(data),sep="")
  }

  # return or write to file
  if (overwrite) {
    write.table(df, fileName, sep=delim, row.names=FALSE, col.names=useHeader, quote=FALSE)
  } else {
    return(df)
  }
}

# -----------------------------------
#' merge.SpatialPolygonsDataFrame
#'
#' Reads in a shape file and data frame to be merged with the data in this shapefile. Merges while preserving the order of objects (ordinary merge operation causes polygons to become disassociated with data).
#'
#' @export

merge.SpatialPolygonsDataFrame <- function(shp, df) {

  # load sp (I think this is the correct package?)
  if (!"sp"%in%rownames(installed.packages()))
    install.packages('sp')
  require(sp)

  # extract data from shapefile and add key
  shp_data <- shp@data
  shp_data$mergeKey <- 1:nrow(shp_data)

  # merge with df and sort based on key
  m <- merge(shp_data,df,all.x=T)
  m <- m[order(m$mergeKey),]
  m <- subset(m,select=-mergeKey)

  # fix row names
  row.names(m) <- row.names(shp)

  # make final SpatialPolygonsDataFrame object
  s <- SpatialPolygonsDataFrame(geometry(shp), m)
  return(s)
}

# -----------------------------------
#' getPolyArea
#'
#' Reads in a shape file and extracts area of every polygon.
#'
#' @export

getPolyArea <- function(shp) {

  polys <- slot(shp, "polygons")
  output <- rep(NA,length(polys))
  for (i in 1:length(polys)) {
    output[i] <- slot(polys[[i]], "area")
  }
  return(output)
}

# -----------------------------------
#' rateRatio
#'
#' Computes point estimate and upper and lower confidence intervals on a ratio of rates. Default method is to enter raw counts and time periods, but if entering rates simply set time1=1 and time2=1.
#'
#' @export

rateRatio <- function(count1, time1, count2, time2, alpha=0.05) {

  point <- time2/time1*count1/count2
  if (count1==0 | time1==0 | count2==0 | time2==0) {
    if (count1==0) {
      return(list(point=0,LL=NaN,UL=NaN))
    } else {
      return(list(point=NaN,LL=NaN,UL=NaN))
    }
  }
  LL <- time2/time1*count1/(count2+1)*qf(0.025,2*(count2+1),2*count1)
  UL <- time2/time1*(count1+1)/count2*qf(1-0.025,2*(count1+1),2*count2)

  return(list(point=point,LL=LL,UL=UL))
}

# -----------------------------------
#' simQuantiles
#'
#' Runs a given stochastic simulation function many times, computing the mean and quantiles over replicates. Note that this method will only work with simulations that have a fixed time step, i.e. synchronous or hybrid simulations, and not with asynchronous simulations. In the hybrid case the maxIterations limit cannot be reached in any simulation.
#'
#' @param FUN the stochastic simulation function to use.
#' @param args a list of arguments to the function.
#' @param reps number of times to repeat the stochastic simulation.
#' @param quantiles which quantiles to compute over replicates.
#'
#' @export

simQuantiles <- function(FUN="SIS_stochastic_hybrid", args=list(), reps=1e2, quantiles=c(0.05,0.5,0.95)) {

  # run function once to get dimensions and variable names
  testOutput <- do.call(FUN, args)
  varNames <- setdiff(names(testOutput),"time")

  # repeat simulation many times and store in array
  simArray <- array(0, dim=c(nrow(testOutput), length(varNames), reps))
  simArray[,,1] <- as.matrix(testOutput[,varNames])
  if (reps>1) {
    for (i in 2:reps) {
      simArray[,,i] <- as.matrix(do.call(FUN, args)[,varNames])
    }
  }

  # compute mean and quantiles over replicates and store in data frame
  df <- data.frame(time=testOutput$time)
  for (i in 1:length(varNames)) {
    m <- rowMeans(simArray[,i,,drop=FALSE])
    q <- apply(simArray[,i,,drop=FALSE], 1, quantile, prob=quantiles)
    df_new <- as.data.frame(t(rbind(m,q)))
    names(df_new) <- paste(varNames[i], c("mean", paste("Q", quantiles, sep="")), sep="_")
    df <- cbind(df, df_new)
  }

  # return summary data frame
  return(df)
}

# -----------------------------------
#' SIS_analytical
#'
#' Returns analytical solution to deterministic SIS model. At equilibrium the number of infectives is given by \eqn{I* = N(1 - r/beta)}.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_analytical <- function(beta=1, r=0.25, I_init=10, N=1e3, times=0:100) {

  I <- (beta-r)/(beta/N + (beta-r-beta*I_init/N)/I_init*exp(-(beta-r)*times))

  output <- data.frame(time=times, S=N-I, I=I)
  return(output)
}

# -----------------------------------
#' SIS_deterministic
#'
#' Returns solution to deterministic SIS model using the \code{odin} package.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_deterministic <- function(beta=1, r=0.25, I_init=10, N=1e3, times=0:100) {

  # solve ode
  mod <- SIS_deterministic_odin2(beta=beta, r=r, I_init=I_init, N=N)
  output <- as.data.frame(mod$run(times))
  names(output)[1] <- 'time'

  return(output)
}

# -----------------------------------
#' SIS_stochastic_async
#'
#' Draw from asynchronous stochastic SIS model. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIS_stochastic_async <- function(beta=1, r=0.25, I_init=100, N=1e3, maxIterations=1e4) {

  # run model
  args <- list(beta=beta, r=r, I_start=I_init, N=N, maxIterations=maxIterations)
  rawOutput <- SIS_stochastic_async_cpp(args)

  # format output object
  I <- rawOutput$I
  S <- N-I
  t_vec <- rawOutput$t
  output <- data.frame(time=t_vec,S=S,I=I)
  output <- subset(output, I>=0)

  return(output)
}

# -----------------------------------
#' SIS_stochastic_hybrid
#'
#' Draw from stochastic SIS model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIS_stochastic_hybrid <- function(beta=1, r=0.25, I_init=100, N=1e3, times=0:100, maxIterations=1e4) {

  # run model
  args <- list(beta=beta, r=r, I_start=I_init, N=N, t_vec=times, maxIterations=maxIterations)
  rawOutput <- SIS_stochastic_hybrid_cpp(args)

  # format output object
  I <- rawOutput$I
  S <- N-I
  S[I<0] <- NA
  I[I<0] <- NA
  output <- data.frame(time=times,S=S,I=I)

  return(output)
}

# -----------------------------------
#' SIS_stochastic_sync
#'
#' Draw from synchronous stochastic SIS model. Return infectives at known time points. Note that results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_stochastic_sync <- function(beta=1, r=0.25, I_init=100, N=1e3, times=0:100) {

  # run model
  args <- list(beta=beta, r=r, I_start=I_init, N=N, t_vec=times)
  rawOutput <- SIS_stochastic_sync_cpp(args)

  # format output object
  I <- rawOutput$I
  S <- N-I
  S[I<0] <- NA
  I[I<0] <- NA
  output <- data.frame(time=times,S=S,I=I)

  return(output)
}

# -----------------------------------
#' SIR_deterministic
#'
#' Returns solution to deterministic SIR model using the \code{odin} package.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIR_deterministic <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, times=0:100) {

  # solve ode
  mod <- SIR_deterministic_odin(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N)
  output <- as.data.frame(mod$run(times))
  names(output)[1] <- 'time'

  return(output)
}

# -----------------------------------
#' SIR_stochastic_async
#'
#' Draw from asynchronous stochastic SIR model. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIR_stochastic_async <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, maxIterations=1e4) {

  # run model
  args <- list(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N, maxIterations=maxIterations)
  rawOutput <- SIR_stochastic_async_cpp(args)

  # format output object
  t_vec <- rawOutput$t
  S <- rawOutput$S
  I <- rawOutput$I
  R <- rawOutput$R
  output <- data.frame(time=t_vec, S=S, I=I, R=R)
  output <- subset(output, I>=0)

  return(output)
}

# -----------------------------------
#' SIR_stochastic_hybrid
#'
#' Draw from stochastic SIR model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIR_stochastic_hybrid <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, times=0:100, maxIterations=1e4) {

  # run model
  args <- list(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N, t_vec=times, maxIterations=maxIterations)
  rawOutput <- SIR_stochastic_hybrid_cpp(args)

  # format output object
  S <- rawOutput$S
  I <- rawOutput$I
  R <- rawOutput$R
  S[S<0] <- NA
  I[I<0] <- NA
  R[R<0] <- NA
  output <- data.frame(time=times, S=S, I=I, R=R)

  return(output)
}

# -----------------------------------
#' SIR_stochastic_sync
#'
#' Draw from synchronous stochastic SIR model. Return state of the system at known time points. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N. Results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIR_stochastic_sync <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, times=0:100) {

  # run model
  args <- list(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N, t_vec=times)
  rawOutput <- SIR_stochastic_sync_cpp(args)

  # format output object
  S <- rawOutput$S
  I <- rawOutput$I
  R <- rawOutput$R
  S[S<0] <- NA
  I[I<0] <- NA
  R[R<0] <- NA
  output <- data.frame(time=times, S=S, I=I, R=R)

  return(output)
}

# -----------------------------------
#' SLIR_deterministic
#'
#' Returns solution to deterministic SLIR model, where L is an incubation (lag) stage of defined length. Solves delay differential equation using the \code{odin} package.
#'
#' @param beta contact rate.
#' @param dur_lag length of time in incubation state.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SLIR_deterministic <- function(beta=0.5, dur_lag=1, r=0.25, I_init=10, R_init=0, N=1e3, times=0:100) {

  # solve ode
  mod <- SLIR_deterministic_odin(beta=beta, dur_lag=dur_lag, r=0.25, I_init=I_init, R_init=R_init, N=N)
  output <- as.data.frame(mod$run(times))
  names(output)[1] <- 'time'

  return(output)
}

# -----------------------------------
#' SLIR_stochastic_async
#'
#' Draw from asynchronous stochastic SLIR model, where L is an incubation (lag) stage of defined length. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
#'
#' @param beta contact rate.
#' @param dur_lag length of time in incubation state.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SLIR_stochastic_async <- function(beta=1, dur_lag=1, r=0.25, I_init=100, R_init=0, N=1e3, maxIterations=1e4) {

  # run model
  args <- list(beta=beta, dur_lag=dur_lag, r=r, I_init=I_init, R_init=R_init, N=N, maxIterations=maxIterations)
  rawOutput <- SLIR_stochastic_async_cpp(args)

  # format output object
  t_vec <- rawOutput$t
  S <- rawOutput$S
  L <- rawOutput$L
  I <- rawOutput$I
  R <- rawOutput$R
  output <- data.frame(time=t_vec, S=S, L=L, I=I, R=R)
  output <- subset(output, I>=0)

  return(output)
}

# -----------------------------------
#' SLIR_stochastic_hybrid
#'
#' Draw from stochastic SLIR model, where L is an incubation (lag) stage of defined length, using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
#'
#' @param beta contact rate.
#' @param dur_lag length of time in incubation state.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SLIR_stochastic_hybrid <- function(beta=1, dur_lag=1, r=0.25, I_init=100, R_init=0, N=1e3, times=0:100, maxIterations=1e4) {

  # run model
  args <- list(beta=beta, dur_lag=dur_lag, r=r, I_init=I_init, R_init=R_init, N=N, t_vec=times, maxIterations=maxIterations)
  rawOutput <- SLIR_stochastic_hybrid_cpp(args)

  # format output object
  S <- rawOutput$S
  L <- rawOutput$L
  I <- rawOutput$I
  R <- rawOutput$R
  S[S<0] <- NA
  L[L<0] <- NA
  I[I<0] <- NA
  R[R<0] <- NA
  output <- data.frame(time=times, S=S, L=L, I=I, R=R)

  return(output)
}


# -----------------------------------
# The following commands are needed to ensure that the roxygen2 package, which deals with documenting the package, does not conflict with the Rcpp package. Do not alter!

#' @useDynLib bobFunctions
#' @importFrom Rcpp evalCpp
NULL



# -----------------------------------
#' header
#'
#' Write header text to console that can be copied into a new script.
#'
#' @export

header <- function(name="Name") {

  # note - do not indent this section or it will show up in header() output!
  # note - indents can be inserted automatically when copying code from one editor to another, so be wary of this.
  s <- paste("
             # ",name,".R
             # Author: Isobel Routledge
             # Date: ",Sys.Date(),"
             # Purpose:
             # (this is an example header)
             # ------------------------------------------------------------------
             # load usefulFunctions package. If not already installed, this can be obtained from github via the devtools command install_github('bobverity/bobFunctions')
             library(usefulFunctions)
             ", sep="")

  cat(s)
}
# -----------------------------------
#' loadPackage
#'
#' Single command for installing (if not already installed) and loading a package. For example, the \code{loadPackage('devtools')} command would run the following code block:
#' \preformatted{
#'	if (!"devtools"\%in\%rownames(installed.packages())) {
#'		install.packages('devtools')
#' }
#' require("devtools")
#' }
#'
#' @export

loadPackage = function(string) {

  # create command string
  s1 = paste('if (!\'',string,'\'%in%rownames(installed.packages())) {',sep='')
  s2 = paste('install.packages(\'',string,'\')',sep='')
  s3 = '}'
  s4 = paste('require(',string,')',sep='')

  # evaluate string
  fullString = paste(s1,s2,s3,s4,sep='\n')
  eval(parse(text=fullString))
}

# -----------------------------------
#' plotOn
#'
#' Start optional save-to-file function.
#'
#' @export

plotOn = function(active=FALSE, fileroot='', filestem1='', filestem2='', fileindex='', type='pdf', width=6, height=6, res=300) {

  if (active) {
    if (type=='pdf') {
      pdf(file=paste(fileroot,filestem1,filestem2,fileindex,".pdf",sep=""), bg="white", width=width, height=height)
    }
    if (type=='png') {
      png(file=paste(fileroot,filestem1,filestem2,fileindex,".png",sep=""), bg="white", width=width*500, height=height*500,res=res)
    }
  }
}

# -----------------------------------
#' plotOff
#'
#' End optional save-to-file function.
#'
#' @export

plotOff = function(active=FALSE) {

  if (active) {
    dev.off()
  }
}

# -----------------------------------
#' zeroPad
#'
#' Add leading zeros to number. Can handle negative numbers. For more advanced options see help for the formatC function.
#'
#' @export

zeroPad = function(number,padding) {

  preamble <- ''
  if (substr(number,1,1)=='-') {
    preamble <- '-'
    number <- substr(number,2,nchar(x))
  }
  output <- paste(paste(rep(0,padding-nchar(number)),collapse=""),number,sep="")
  output <- paste(preamble,output,sep='')
  return(output)
}

# -----------------------------------
#' rdirichlet
#'
#' Draw from Dirichlet distribution with parameters alpha_vec (does not have to be symmetric).
#'
#' @export

rdirichlet <- function(alpha_vec) {

  Y <- rgamma(length(alpha_vec), shape=alpha_vec, scale=1)
  output <- Y/sum(Y)
  return(output)
}

# -----------------------------------
#' rdirichlet
#'
#' Draw from Dirichlet distribution with parameters alpha_vec (does not have to be symmetric).
#'
#' @export

rdirichlet <- function(alpha_vec) {

  Y <- rgamma(length(alpha_vec), shape=alpha_vec, scale=1)
  output <- Y/sum(Y)
  return(output)
}

# -----------------------------------
#' rinvgamma
#'
#' Draw from inverse-gamma distribution
#'
#' @export

rinvgamma <- function(n, alpha, beta) {
  ret <- 1/rgamma(n,shape=alpha,rate=beta)
  return(ret)
}

# -----------------------------------
#' dinvgamma
#'
#' Density of inverse-gamma distribution
#'
#' @export

dinvgamma <- function(x, alpha, beta, log=TRUE) {
  ret <- alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x
  if (log) {
    return(ret)
  } else {
    return(exp(ret))
  }
}

# -----------------------------------
#' dt_scaled
#'
#' Density of scaled student's t distribution.
#'
#' @export

dt_scaled = function(x, df, ncp, scale, log=FALSE) {

  output <- lgamma((df+1)/2)-lgamma(df/2)-0.5*log(pi*df*scale^2)-((df+1)/2)*log(1+1/df*((x-ncp)/scale)^2)
  if (log==FALSE)
    output <- exp(output)
  return(output)
}

# -----------------------------------
#' normal_loglike
#'
#' Log-probability of data x integrated over normal likelihood with standard deviation sigma and normal prior with mean priorMean and standard deviation priorSD. Probability of data is raised to the power beta, which has default value of beta=1 giving ordinary likelihood.
#'
#' Model specification:
#' x_i ~ Normal(mu, sigma)
#' mu ~ Normal(priorMean, priorSD)
#'
#' @export

normal_loglike <- function(x, sigma, priorMean, priorSD, beta=1) {
  n <- length(x)
  mu_postVar <- 1/(beta*n/sigma^2+1/priorSD^2)
  mu_postMean <- mu_postVar * (beta*sum(x)/sigma^2+priorMean/priorSD^2)
  ret <- -0.5*beta*n*log(2*pi*sigma^2) - 0.5*log(priorSD^2) + 0.5*log(mu_postVar) - 0.5*(beta*sum(x^2)/sigma^2 + priorMean^2/priorSD^2 - mu_postMean^2/mu_postVar)
  return(ret)
}

# -----------------------------------
#' logDescriptive
#'
#' Calculates descriptive stats on extremely large or extremely small values. Input vector of values in log-space. Output sample mean and variance, also in log-space.
#'
#' @export

logDescriptive = function(Y) {

  n <- length(Y)
  log_ybar <- log(mean(exp(Y-min(Y))))+min(Y)
  log_samplevar <- 2*min(Y) + log(n/(n-1)*mean(exp(2*Y-2*min(Y))-2*exp(Y+log_ybar-2*min(Y))+exp(2*log_ybar-2*min(Y))))
  return(c(log_ybar,log_samplevar))
}

# -----------------------------------
#' dBDI
#'
#' Probability mass function for the birth-death-immigration model. Note that this is simply a particular version of the negative-binomial likelihood. When applied to gene families birthRate=rate of gene duplication, deathRate=rate of gene deletion, immigrationRate=rate of de novo gene creation.
#'
#' @param x number of individuals or genes (i.e. value of random variable).
#' @param birthRate rate of birth.
#' @param deathRate rate of death.
#' @param immigrationRate rate of immigration.
#'
#' @export

dBDI = function(x, birthRate, deathRate, immigrationRate, log=FALSE) {

  dnbinom(x=x, size=immigrationRate/birthRate, prob=1-birthRate/deathRate, log=log)

}

# -----------------------------------
#' dlgamma
#'
#' Probability density function for log(X), where X is gamma(shape=alpha,rate=beta). The mean of this distribution is equal to the integral of log(x)*dgamma(x) over the interval [0,infinity), which does not have a simple analytical solution that I am aware of.
#'
#' @export

dlgamma <- function(x, alpha, beta, log=FALSE) {

  output <- alpha*x+alpha*log(beta)-lgamma(alpha)-beta*exp(x)
  if (log==FALSE)
    output <- exp(output)
  return(output)
}

# -----------------------------------
#' logSum
#'
#' Add numbers together in log space while avoiding under/overflow problems. Accepts vector inputs.
#'
#' @export

logSum <- function(logA, logB) {
  output <- rep(0,length(logA))

  output[logA<logB] <- ( logB + log(1 + exp(logA-logB)) )[logA<logB]
  output[logB<=logA] <- (logA + log(1 + exp(logB-logA)) )[logB<=logA]

  return(output)
}


# -----------------------------------
#' allSamps
#'
#' Output all possible sequences that are possible by sampling without replacement from a given list.
#'
#' @export

allSamps <- function(targetlist) {

  lengths <- mapply(length,targetlist)
  outputmat <- matrix(0, nrow=prod(lengths), ncol=length(lengths))
  for (i in 1:length(lengths)) {
    val1 <- prod(lengths[1:i][-i])
    val2 <- prod(lengths[i:length(lengths)][-1])
    outputmat[,i] <- rep(rep(targetlist[[i]],each=val2),times=val1)
  }
  return(outputmat)
}

# -----------------------------------
#' convertRadix
#'
#' Convert decimal number system to any other radix (base).
#'
#' @param x number to convert.
#' @param base radix to convert to.
#' @param fixlength if used, always output at least this many digits.
#'
#' @export

convertRadix <- function(x, base=2, fixlength=NA) {

  floornumber <- x
  output <- NULL
  while (floornumber!=0) {
    output <- c(floornumber%%base,output)
    floornumber <- floor(floornumber/base)
  }
  if (!is.na(fixlength) & fixlength>length(output)) {
    output <- c(rep(0,fixlength-length(output)),output)
  }
  return(output)
}


# -----------------------------------
#' rCRP
#'
#' Draw group and group frequencies from a Chinese restaurant process with concentration parameter theta.
#'
#' @export

rCRP <- function(n, theta=1) {

  if (n==1)
    return(list(group=1,groupFreqs=1))

  group <- rep(1,n)
  maxGroup <- 1
  groupFreqs <- rep(0,n)
  groupFreqs[1] <- 1
  probVec <- groupFreqs
  probVec[maxGroup+1] <- theta
  for (i in 2:n) {
    newGroup <- sample(n,1,prob=probVec)
    group[i] <- newGroup
    if (newGroup>maxGroup) {
      maxGroup <- newGroup
      groupFreqs[maxGroup] <- 1
      probVec[maxGroup] <- 1
      probVec[maxGroup+1] <- theta
    } else {
      groupFreqs[newGroup] <- groupFreqs[newGroup]+1
      probVec[newGroup] <- probVec[newGroup]+1
    }
  }
  return(list(group=group,groupFreqs=groupFreqs))
}

# -----------------------------------
#' rCRP2
#'
#' Draw from Chinese restaurant process multiple times using efficient stick-breaking construction.
#'
#' @export

rCRP2 = function(n, reps, theta=1) {

  # initialise objects
  n_left <- rep(n,reps)
  freqs <- NULL

  # draw from stick-breaking process until no new groups
  while (any(n_left>0)) {
    p <- rbeta(reps,1,theta)
    newFreqs <- rbinom(reps,n_left,p)
    freqs <- cbind(freqs,newFreqs)
    n_left <- n_left - newFreqs
  }

  # drop all-zero columns
  freqs <- freqs[,colSums(freqs)>0,drop=FALSE]

  return(freqs)
}

# -----------------------------------
#' colVars
#'
#' Calculate variance (sample or population) of columns of a matrix or data frame. If sampleVar==FALSE then calcualte relative to mean mu.
#'
#' @export

colVars <- function(mat, sampleVar=TRUE, mu=0) {

  n <- colSums(!is.na(mat))
  X <- colSums(mat,na.rm=T)
  X2 <- colSums(mat^2,na.rm=T)
  if (sampleVar) {
    output <- 1/(n-1)*(X2-X^2/n)
  } else {
    output <- 1/n*(X2-2*mu*X+n*mu^2)
  }
  return(output)
}

# -----------------------------------
#' rowVars
#'
#' Calculate variance (sample or population) of rows of a matrix or data frame. If sampleVar==FALSE then calcualte relative to mean mu.
#'
#' @export

rowVars = function(mat, sampleVar=TRUE, mu=0) {

  n <- rowSums(!is.na(mat))
  X <- rowSums(mat,na.rm=T)
  X2 <- rowSums(mat^2,na.rm=T)
  if (sampleVar) {
    output <- 1/(n-1)*(X2-X^2/n)
  } else {
    output <- 1/n*(X2-2*mu*X+n*mu^2)
  }
  return(output)
}


# -----------------------------------
#' grad
#'
#' Numerically calculate gradient at each point (x,y). The gradient at a point is calculated as the mean of the two gradients in either direction from the point, apart from at the ends of the vector where there is only a single value to consider.
#'
#' @export

grad <- function(x=1:length(y), y) {

  # check that same length
  if (length(y)!=length(x))
    stop('x and y must be same length')

  # calculate gradient from differences in x and y direction
  delta_y <- y[-1]-y[-length(y)]
  delta_x <- x[-1]-x[-length(x)]
  grad <- delta_y/delta_x

  # gradient at each point is mean of gradient in each direction from point (apart from at ends).
  grad_left <- c(grad[1],grad)
  grad_right <- c(grad,grad[length(grad)])
  grad_mean <- (grad_left+grad_right)/2

  return(grad_mean)
}

# -----------------------------------
#' cubeSpline_segment
#'
#' Calculate cubic spline of the form ax^3+bx^2+cx+d between the points (x[1],y[1]) and (x[2],y[2]) with gradients g[1] and g[2]. Return constants a, b, c and d.
#'
#' @export

cubeSpline_segment <- function(x, y, g) {

  # calculate coefficients
  a <- (y[2]-y[1]-0.5*(g[1]+g[2])*(x[2]-x[1]))/(2*x[1]^3-3*x[1]^2*x[2]+x[2]^3-0.5*(3*x[2]^2-3*x[1]^2)*(x[2]-x[1]))
  b <- (g[2]-g[1]-a*(3*x[2]^2-3*x[1]^2))/(2*x[2]-2*x[1])
  c <- g[1]-3*a*x[1]^2-2*b*x[1]
  d <- y[2]-a*x[2]^3-b*x[2]^2-g[1]*x[2]+3*a*x[1]^2*x[2]+2*b*x[1]*x[2]

  # return results
  return(list('a'=a,'b'=b,'c'=c,'d'=d))
}

# -----------------------------------
#' cubeSpline
#'
#' Calculate cubic spline passing through the points (x,y) with known gradients g at these points. The parameter 'interpoints' describes the number of points in each segment.
#'
#' @export

cubicSpline <- function(x, y, g, interpoints=10) {

  #xvec <- seq(x[1],x[n],l=(n-1)*interpoints)
  #yvec <- rep(0,(n-1)*interpoints)

  # interpolate each segment using cubic spline
  n <- length(x)
  yout <- xout <- NULL
  for (i in 2:n) {
    #whichPoints <- ((i-2)*interpoints+1):((i-1)*interpoints)
    coeffs <- cubeSpline_segment(x[(i-1):i], y[(i-1):i], g[(i-1):i])
    xvec <- seq(x[i-1], x[i], l=interpoints+1)
    #yvec[whichPoints] <- z$a*xvec[whichPoints]^3+z$b*xvec[whichPoints]^2+z$c*xvec[whichPoints]+z$d
    yvec <- coeffs$a*xvec^3 + coeffs$b*xvec^2 + coeffs$c*xvec + coeffs$d

    if (i==n) {
      xout <- c(xout, xvec)
      yout <- c(yout, yvec)
    } else {
      xout <- c(xout, xvec[-length(xvec)])
      yout <- c(yout, yvec[-length(yvec)])
    }
  }

  return(list('x'=xout,'y'=yout))
}

# -----------------------------------
#' logHarmonicMean
#'
#' Calculates harmonic mean of values z, where values are in log space. Accounts for underflow and returns value of log harmonic mean.
#'
#' @export

logHarmonicMean <- function(z) {

  m <- min(z)
  output <- m - log(sum(exp(-(z-m)))) + log(length(z))
  return(output)
}

# -----------------------------------
#' MCMCPlot
#'
#' Produces simple plot useful for visualising MCMC chains. Enter x vector to plot single chain, or x and y vectors to visualise correlation between chains.
#'
#' @export

MCMCPlot <- function(x, y=NULL, col=grey(0.7), pch=4, cex=0.4, xmin=NULL, xmax=NULL, ymin=min(x,na.rm=T), ymax=max(x,na.rm=T), xlab=NULL, main='', ylab=main) {

  if (is.null(y)) {
    if (is.null(xmin))
      xmin <- 1
    if (is.null(xmax))
      xmax <- length(x)
    if (is.null(xlab))
      xlab <- 'iteration'
    plot(x, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=xlab, ylab=ylab, main=main, pch=pch, cex=cex, col=col)
  } else {
    if (is.null(xmin))
      xmin <- min(x,na.rm=T)
    if (is.null(xmax))
      xmax <- max(x,na.rm=T)
    if (is.null(xlab))
      xlab <- ''
    plot(x, y, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=xlab, ylab=ylab, main=main, pch=pch, cex=cex, col=col)
  }
}

# -----------------------------------
#' MCMCPlot2
#'
#' Similar to \code{MCMCPlot()}, but aimed at very large numbers of samples where ordinary plotting would not be feasible. Works by binning points before plotting. Can only handle a single chain (cannot plot x vs. y).
#'
#' @export

MCMCPlot2 <- function(x, h_cells=100, v_cells=100, xmin=1, xmax=length(x), ymin=min(x,na.rm=T), ymax=max(x,na.rm=T), xlab='iteration', main='', ylab=main) {

  # make x,y data frame
  df <- data.frame(x=1:length(x),y=x)

  # make vectors of breaks and matrix z for storing final results
  v_breaks <- seq(ymin, ymax, l=v_cells+1)
  h_breaks <- seq(xmin, xmax, l=h_cells+1)
  v_mids <- (v_breaks[-1]+v_breaks[-length(v_breaks)])/2
  h_mids <- (h_breaks[-1]+h_breaks[-length(h_breaks)])/2
  z <- matrix(0,v_cells,h_cells)

  # split data based on vertical breaks
  df_split <- split(df,f=cut(df$y,v_breaks))

  # for each element, count frequency classes based on horizontal breaks
  for (i in 1:length(df_split)) {
    df_split2 <- split(df_split[[i]],f=cut(df_split[[i]]$x,h_breaks))
    z[i,] <- mapply(nrow,df_split2)
  }


  # plot final matrix
  palette <- colorRampPalette(colors()[c(1,131,107)])
  image(h_mids, v_mids, t(z), col=palette(100), xlab=xlab, ylab=ylab, main=main)
}

# -----------------------------------
#' densityPlot
#'
#' Produce kernel density plot of data x, with mean, median and 95% limits shown.
#'
#' @export

densityPlot <- function(x, xmin=NULL, xmax=NULL, ymin=0, ymax=NULL, xlab='', ylab='density', main='', meanOn=TRUE, medianOn=TRUE, quantileOn=TRUE, boxOn=TRUE, boxSize=0.8, meanTextOn=TRUE, medianTextOn=TRUE, quantileTextOn=TRUE, boxSigFig=3) {

  # get x limits from data if not specified
  if (is.null(xmin))
    xmin <- 1.5*min(x,na.rm=T)-0.5*max(x,na.rm=T)
  if (is.null(xmax))
    xmax <- 1.5*max(x,na.rm=T)-0.5*min(x,na.rm=T)

  # produce kernel density
  fx <- density(x,from=xmin,to=xmax)
  if (is.null(ymax))
    ymax <- 1.1*max(fx$y,na.rm=T)
  plot(fx, ylim=c(ymin,ymax), xlab=xlab, ylab=ylab, main=main)
  polygon(c(fx$x,rev(fx$x)), c(fx$y,rep(0,length(fx$y))), col=grey(0.7))

  # add lines for median and mean
  if (meanOn) {
    best_mean <- which.min(abs(fx$x-mean(x)))
    lines(rep(fx$x[best_mean],2),c(0,fx$y[best_mean]),lty=2)
  }
  if (medianOn) {
    best_median <- which.min(abs(fx$x-median(x)))
    lines(rep(fx$x[best_median],2),c(0,fx$y[best_median]))
  }

  # add coloured quantiles to density plot
  if (quantileOn) {
    quantiles <- quantile(x,prob=c(0.025,0.975))
    fx_left_x <- fx$x[fx$x<quantiles[1]]
    fx_left_y <- fx$y[fx$x<quantiles[1]]
    polygon(c(fx_left_x,rev(fx_left_x)), c(fx_left_y,rep(0,length(fx_left_y))), col=colors()[373])
    fx_right_x <- fx$x[fx$x>quantiles[2]]
    fx_right_y <- fx$y[fx$x>quantiles[2]]
    polygon(c(fx_right_x,rev(fx_right_x)), c(fx_right_y,rep(0,length(fx_right_y))), col=colors()[373])
  }

  # add legend
  if (boxOn) {

    legendText <- NULL
    legend_lty <- NULL
    if (meanTextOn & meanOn) {
      legendText <- c(legendText, paste('   mean=', signif(mean(x), boxSigFig), sep=''))
      legend_lty <- c(legend_lty, 2)
    }
    if (medianTextOn & medianOn) {
      legendText <- c(legendText, paste('median=', signif(median(x), boxSigFig), sep=''))
      legend_lty <- c(legend_lty, 1)
    }
    if (quantileTextOn & quantileOn) {
      legendText <- c(legendText, paste(' q0.025=', signif(quantiles[1], boxSigFig), sep=''))
      legendText <- c(legendText, paste(' q0.975=', signif(quantiles[2], boxSigFig), sep=''))
      legend_lty <- c(legend_lty, NA, NA)
    }

    legend('topright', legend=legendText, lty=legend_lty, seg.len=1, bg='#FFFFFF80', box.col='#00000050', cex=boxSize)
  }
}




# -----------------------------------
#' colPlot
#'
#' Produce simple plot to visualise a color palette.
#'
#' @export

colPlot <- function(colVec, size=8) {

  n <- length(colVec)
  plot(1:n, rep(0,n), col=colVec, pch=20, cex=size, axes=FALSE, xlab='', ylab='')
  text(1:n,rep(0,n),labels=1:n)
}

# -----------------------------------
#' colPlot
#'
#' Take vector of colours c in hexadecimal form and build in transparancy alpha.
#'
#' @export

transHex <- function(c, alpha) {

  z <- t(col2rgb(c))/255
  output <- mapply(function(x) {rgb(x[1],x[2],x[3],alpha=alpha)}, split(z,f=row(z)))
  return(output)
}

# -----------------------------------
#' first
#'
#' Return first value in a vector, or alternatively trim off the first value in a vector.
#'
#' @export

first <- function(x, trim=FALSE) {
  if (trim) {
    return(x[-1])
  } else {
    return(x[1])
  }
}

# -----------------------------------
#' last
#'
#' Return last value in a vector, or alternatively trim off the last value in a vector.
#'
#' @export

last <- function(x, trim=FALSE) {
  if (trim) {
    return(x[-length(x)])
  } else {
    return(x[length(x)])
  }
}

# -----------------------------------
#' latexTable
#'
#' Convert data frame x to LaTeX format with set decimal places, and write output to file.
#'
#' @export

latexTable <- function(x, digits=3, file='~/Desktop/tester.txt') {

  outVec <- rep(NA,nrow(x))
  s <- paste("%.",digits,"f",sep='')
  for (i in 1:nrow(x)) {
    outVec[i] <- paste(paste(sprintf(s,x[i,]),collapse=' & '),'\\\\')
  }
  write.table(outVec, file=file, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# -----------------------------------
#' smoothCols
#'
#' Read in continuous values between xmin and xmax and return colours associated with these values, taken from a smoothly varying scale.
#'
#' @param x values from which to obtain colours.
#' @param xmin minimum range of x.
#' @param xmax maximum range of x.
#' @param res number of colours in palette.
#' @param rawCols colours that make up the palette. Leave as NULL to use default values, taken from tim.colors().
#'
#' @export

smoothCols <- function(x, xmin=min(x,na.rm=T), xmax=max(x,na.rm=T), res=1e3, rawCols=NULL) {

  # load RColorBrewer if needed
  loadPackage('RColorBrewer')

  # get x between 0 and 1
  x <- (x-xmin)/(xmax-xmin)

  # define colours if needed
  if (is.null(rawCols))
    rawCols <- c('#00008F', '#00009F', '#0000AF', '#0000BF', '#0000CF', '#0000DF', '#0000EF', '#0000FF', '#0010FF', '#0020FF', '#0030FF', '#0040FF', '#0050FF', '#0060FF', '#0070FF', '#0080FF', '#008FFF', '#009FFF', '#00AFFF', '#00BFFF', '#00CFFF', '#00DFFF', '#00EFFF', '#00FFFF', '#10FFEF', '#20FFDF', '#30FFCF', '#40FFBF', '#50FFAF', '#60FF9F', '#70FF8F', '#80FF80', '#8FFF70', '#9FFF60', '#AFFF50', '#BFFF40', '#CFFF30', '#DFFF20', '#EFFF10', '#FFFF00', '#FFEF00', '#FFDF00', '#FFCF00', '#FFBF00', '#FFAF00', '#FF9F00', '#FF8F00', '#FF8000', '#FF7000', '#FF6000', '#FF5000', '#FF4000', '#FF3000', '#FF2000', '#FF1000', '#FF0000', '#EF0000', '#DF0000', '#CF0000', '#BF0000', '#AF0000', '#9F0000', '#8F0000', '#800000')

  # make smooth colours
  myPal <- colorRampPalette(rawCols)
  cols <- myPal(res+1)[floor(x*res)+1]

  return(cols)
}

# -----------------------------------
#' changeNames
#'
#' Change names of selected variables in a data frame. Note that for large data frames the raw code of this function should be copied into script rather than using this function, as there will be a memory cost in copying the data frame within the function.
#'
#' @export

changeNames <- function(df, oldNames=names(df), newNames) {

  names(df)[names(df)%in%oldNames] <- newNames
  return(df)
}

# -----------------------------------
#' reportConsole
#'
#' Use within a loop to print current iteration to console.
#'
#' @export

reportConsole <- function(i,imax,istep=1) {

  if (i%%istep==0) {
    s <- paste("iteration ",i," of ",imax,", (",round(i/imax*100),"%)\n",sep="")
    cat(s)
    flush.console()
  }
}

# -----------------------------------
#' dot
#'
#' Prints one or more dots to screen, useful for tracking progress within nested loops. Follow by newLine() to add carriage return.
#'
#' @export

dot <- function(n=1) {

  cat(paste(rep('.',n),collapse=''))
  flush.console()
}

# -----------------------------------
#' newLine
#'
#' Prints one or more carriage returns.
#'
#' @export

newLine <- function(n=1) {

  cat(paste(rep('\n',n),collapse=''))
  flush.console()
}

# -----------------------------------
#' monthDays
#'
#' Return days in each month of chosen year.
#'
#' @export

monthDays <- function(year) {

  date1 <- paste(year,"-01-01",sep="")
  date2 <- paste(year+1,"-01-01",sep="")
  output <- as.numeric(diff(seq(as.Date(date1),as.Date(date2),by="month")))
  return(output)
}

# -----------------------------------
#' is.int
#'
#' Check that character string can be converted to integer without issue (for example, the value 3.4 would return FALSE).
#'
#' @export

is.int <- function(x) {

  output <- FALSE
  if (!is.na(suppressWarnings(as.numeric(x))))
    output <- (as.numeric(x)%%1==0)
  return(output)
}

# -----------------------------------
#' errorBars
#'
#' Produce error bars at given positions. When \code{se=FALSE}, input should be \code{y1}=lower limit and \code{y2}=upper limit. When \code{se=TRUE}, input should be \code{y1}=mean and \code{y2}=standard error, in which case error bars represent two standard errors either side of the mean. Can handle vector inputs as well as scalars.
#'
#' @param y1 under first method (when \code{se=FALSE}) this value is the lower limit of the error bar. Under the second method (when \code{se=TRUE}) this value is the mean of the distribution that will be used to produce error bars (mean +- 1.96 standard deviations).
#' @param y2 as above, either the upper limit of the error bar, or the standard deviation of the distribution that will be used to produce error bars.
#' @param x horizontal position of error bars.
#' @param se whether to use values \code{y1} and \code{y2} as mean and standard deviation.
#' @param width length of error bar handles.
#' @param col colour of error bars.
#' @param lty line type of error bars.
#' @param lwd line width of error bars.
#'
#' @export

errorBars <- function(y1, y2, x=1:length(y1), se=FALSE, width=1, col=1, lty=1, lwd=1) {

  # calculate limits based on method
  if (se) {
    LL <- y1-2*y2
    UL <- y1+2*y2
  } else {
    LL <- y1
    UL <- y2
  }

  segments(x,LL,x,UL,lty=lty,lwd=lwd,col=col)
  segments(x-width/2,LL,x+width/2,LL,lty=lty,lwd=lwd,col=col)
  segments(x-width/2,UL,x+width/2,UL,lty=lty,lwd=lwd,col=col)
}

# -----------------------------------
#' incrementBins
#'
#' Increment vector x under binMax constraint. Return NULL at final increment.
#'
#' @export

incrementBins <- function(x,binMax=c(1,4,3,2)) {

  if (length(x)!=length(binMax))
    stop("x and binMax of different lengths")
  if (any(x>binMax))
    stop("some x greater than binMax")
  if (all(x==binMax))
    return()
  for (i in length(x):1) {
    x[i] <- x[i]+1
    if (x[i]<=binMax[i]) {
      break
    } else {
      x[i] <- 1
    }
  }
  return(x)
}

# -----------------------------------
#' removeCols
#'
#' Remove named columns from data frame.
#'
#' @export

removeCols <- function(df, nameVec) {

  output <- subset(df,select=setdiff(names(df),nameVec))
  return(output)
}

# -----------------------------------
#' matrixSmooth
#'
#' Interpolate rows and columns of a matrix alternately using spline. reps is the number of times to double the matrix size.
#'
#' @export

matrixSmooth <- function(M, reps=1) {

  for (i in 1:reps) {
    M2 <- cbind(M,M)
    for (i in 1:nrow(M)) {
      M2[i,] <- spline(M[i,],n=2*length(M[i,]))$y
    }
    M <- rbind(M2,M2)
    for (i in 1:ncol(M2)) {
      M[,i] <- spline(M2[,i],n=2*length(M2[,i]))$y
    }
  }
  return(M)
}

# -----------------------------------
#' minSpanTree
#'
#' Use Prim's algorithm to calculate minimum spanning tree. The 'data' argument must be a matrix or data frame with observations in rows and dimensions in columns. Distances are calculated as Euclidean distance over all the dimensions of the data. Returns multiple objects; 1) a data frame of nodePairs giving all pairwise links that make up the tree in the order that they were added, 2) an edgeList giving the nodes attached to each node, 3) edgeNum giving the number of edges associated with each node.
#' \cr\cr
#' Although this function may not be as fast as some alternatives (not tested), it has the advantage that distances are only calculated as needed, meaning there is no need to read in a huge matrix of pairwise distances.
#'
#' @export

minSpanTree <- function(data) {

  # extract basic properties of data
  data <- as.matrix(data)
  n <- nrow(data)
  dims <- ncol(data)

  # initialise objects for storing results
  nodePairs <- data.frame(node1=rep(0,n-1),node2=0)
  edgeList <- replicate(n,NULL)
  edgeNum <- rep(0,n)

  # initialise algorithm by measuring all distances relative to first point
  point1 <- data[1,]
  data <- data[-1,]
  d_running <- 0
  for (j in 1:dims) {
    d_running <- d_running + (data[,j]-point1[j])^2
  }

  # given a point A and a cluster B, the distance between A and B is equal to the *minimum* distance between A and any point in B. The vector bestNode stores *which* node in B the point A links to
  bestNode <- rep(1,n-1)
  # during the algorithm the data object is trimmed, and so the row index no longer stores the unique ID of the data point. The vector pos fixes this by maintaining a record of unique IDs.
  pos <- 2:n

  # run Prims algorithm
  for (i in 1:(n-1)) {
    # node2 is the node with the minimum distance to the current tree. node1 is the node within the tree that node2 links to
    nodeIndex <- which.min(d_running)
    node1 <- bestNode[nodeIndex]
    node2 <- pos[nodeIndex]

    # store nodes and edges
    nodePairs$node1[i] <- node1
    nodePairs$node2[i] <- node2
    edgeList[[node1]] <- c(edgeList[[node1]], node2)
    edgeList[[node2]] <- c(edgeList[[node2]], node1)
    edgeNum[node1] <- edgeNum[node1]+1
    edgeNum[node2] <- edgeNum[node2]+1

    # calculate new distance of all points from node2
    d_new <- 0
    for (j in 1:dims) {
      d_new <- d_new + (data[,j]-data[nodeIndex,j])^2
    }
    # recalculate the *minimum* distance between the current tree and all nodes, and if the new distance is shorter then reflect that in bestNode
    bestNode[d_running>d_new] <- node2
    d_running <- (d_running<=d_new)*d_running + (d_running>d_new)*d_new

    # drop new node from data set and all other objects
    data <- data[-nodeIndex,,drop=FALSE]
    d_running <- d_running[-nodeIndex]
    pos <- pos[-nodeIndex]
    bestNode <- bestNode[-nodeIndex]
  }
  # return multiple objects
  return(list('nodePairs'=nodePairs, 'edgeList'=edgeList, 'edgeNum'=edgeNum))
}


# -----------------------------------
#' pieCharts
#'
#' Adds pie charts to an existing plot. Proportions do not need to sum to one.
#'
#' @export

pieCharts <- function(x, y, proportions, radius=0.2, lwd=1, border_col=1, seg_col=bobRainbow(), smoothness=50, x_stretch=1) {

  theta <- seq(0,2*pi,l=smoothness)
  for (i in 1:length(x)) {
    segs <- length(proportions[[i]])
    z <- c(0, cumsum(proportions[[i]]/sum(proportions[[i]])))
    if (segs==1) {
      polygon(x[i]+x_stretch*radius*cos(theta), y[i]+radius*sin(theta), lwd=lwd, border=border_col, col=seg_col[1])
    } else {
      for (j in 2:(segs+1)) {
        theta_sub <- c(2*pi*z[j-1], theta[theta>(2*pi*z[j-1]) & theta<(2*pi*z[j])], 2*pi*z[j])
        polygon(c(x[i],x[i]+x_stretch*radius*sin(theta_sub),x[i]), c(y[i],y[i]+radius*cos(theta_sub),y[i]), lwd=lwd, border=border_col, col=seg_col[j-1])
      }
    }
  }
}

# -----------------------------------
#' win
#'
#' Shortcut for running \code{par(mfrow=c(rows,cols))}, which takes slightly longer to type!
#'
#' @export

win = function(rows=1, cols=1) {
  par(mfrow=c(rows,cols))
}

# -----------------------------------
#' animateOn
#'
#' First of two functions needed to create mpg from static images. Run this function, then run the code needed to create a sequence of plots, then run \code{animateOff()}. Individual plots will be saved as jpg with a four-character number on the end (for example myFile0012.jpg), and can either be saved in the local directory or in a temporary directory. The \code{animateOff()} function should have the same shared arguments as \code{animateOn()} (for example the same fileName).
#' \cr\cr
#' Note that this function requires ImageMagick to be installed, and has only been tested on Windows. If running for the first time you will likely need to add ImageMagick to the path using something like the following:
#' \code{Sys.setenv(PATH=paste(Sys.getenv("PATH"),"C:\\Program Files\\ImageMagick-6.9.3-Q16",sep=";"))}
#'
#' @param active allows function to be selectively de-activated using a variable.
#' @param fileName the filename of input files (without extension or number). The same name is used for the final video file.
#' @param saveFrames if FALSE then images are saved to a temporary location before being deleted in \code{animateOff()}.
#' @param width width of video.
#' @param height height of video.
#' @param units units of jpg files that make up video (see ?jpeg).
#' @param pointsize pointsize of jpg files that make up video (see ?jpg).
#' @param quality quality of jpg files that make up video (see ?jpg).
#' @param bg background of jpg files that make up video (see ?jpg).
#' @param res resolution of jpg files that make up video (see ?jpg).
#'
#' @export

animateOn <- function(active=FALSE, fileName, saveFrames=FALSE, width=480, height=480, units="px", pointsize=12, quality=75, bg="white", res=NA) {

  # exit if not active
  if (!active)
    return()

  # set file path locally or in temp file
  filePath <- paste(fileName,"%04d.jpeg",sep="")
  if (!saveFrames)
    filePath <- paste(tempdir(),"/",filePath,sep="")

  # open jpg filestream
  jpeg(filePath, width=width, height=height, units=units, pointsize=pointsize, quality=quality, bg=bg, res=res)
}

# -----------------------------------
#' animateOff
#'
#' Second of two functions needed to create mpg from static images. This function should be run after running \code{animateOn()} and creating a series of plots. The \code{animateOff()} function should have the same shared arguments as \code{animateOn()}.
#'
#' @export

animateOff <- function(active=FALSE, fileName, saveFrames=FALSE, frameRate=25) {

  # exit if not active
  if (!active)
    return()

  # close file stream
  dev.off()

  # convert images to mpg
  filePath <- paste(fileName,"%04d.jpeg",sep="")
  if (!saveFrames)
    filePath <- paste(tempdir(),"/",filePath,sep="")
  myCommand <- paste("ffmpeg -y -r ",frameRate," -i ",filePath," ",fileName,".mp4",sep="")
  shell(myCommand,intern=TRUE)
}

# -----------------------------------
#' coordText
#'
#' Defines a square region from 0-1 in both x and y and inserts text at designated location. Useful for e.g. adding panel labels to figures.
#'
#' @export

coordText <- function(x=0.05, y=0.95, text='A)', cex=1.5) {

  pars <- par(new=T,xpd=T,mar=c(0,0,0,0))
  plot(0, type='n', axes=F, xlab=NA, ylab=NA, xlim=c(0,1), ylim=c(0,1))
  text(x,y,text,cex=cex)
  par(pars)
}

# -----------------------------------
#' fastRead
#'
#' Shortcut function for reading in data quickly using the data.table package.
#'
#' @export

fastRead <- function(fileName, header='auto') {

  loadPackage('data.table')
  data <- as.data.frame(fread(fileName, header=header))
  return(data)
}

# -----------------------------------
#' multiPanel
#'
#' This function is not really a proper function (it is not intended to return anything). Rather, it is a place to store this useful bit of code that can adapted as needed. Copy this code over, and then drop whatever individual plots are needed into the inner loop.
#'
#' @export

multiPanel <- function() {

  # setup multi-panel plotting parameters
  plot_rows <- 2
  plot_cols <- 2
  outer_margin <- c(4,4,2,1)
  tickLength <- -0.5
  x_range <- c(-5,5)
  y_range <- c(0,300)
  x_axis_at <- seq(-4,4,2)
  y_axis_at <- seq(0,250,50)
  sub_main <- paste('submain',1:(plot_rows*plot_cols))
  x_lab <- 'x axis'
  y_lab <- 'y axis'
  x_cex <- 0.8
  y_cex <- 0.8

  # create multi-panel plot
  par_store <- par(mfrow=c(plot_rows, plot_cols), mar=c(0,0,0,0), oma=outer_margin, tcl=tickLength)
  index <- 0
  for (i in 1:plot_rows) {
    for (j in 1:plot_cols) {
      index <- index+1

      # put individual plots here
      hist(rnorm(1e3), xlim=x_range, ylim=y_range, ann=FALSE, axes=FALSE)
      title(sub_main[index],line=-1)

      # add axes and box
      if (i==plot_rows)
        axis(1, at=x_axis_at)
      if (j==1)
        axis(2, at=y_axis_at)
      box()
    }
  }
  mtext(x_lab, side=1, outer=TRUE, cex=x_cex, line=2.2)
  mtext(y_lab, side=2, outer=TRUE, cex=y_cex, line=2.2)
  par(par_store)

}

# -----------------------------------
#' imageFix
#'
#' Produces an image plot, but has option for changing orientation. Values of \code{orientation} from 1 to 4 rotate the image clockwise.
#'
#' @export

imageFix <- function(z, x=1:nrow(z), y=1:ncol(z), orientation=1, ...) {

  if (orientation==1) {
    image(x, y, z, ...)
  } else if (orientation==2) {
    image(y, x, t(z[nrow(z):1,]), ...)
  } else if (orientation==3) {
    image(x, y, z[nrow(z):1,ncol(z):1], ...)
  } else if (orientation==4) {
    image(y, x, t(z[,ncol(z):1]), ...)
  }

}

# -----------------------------------
#' filledContour2
#'
#' Produces pretty alternative to ordinary filled contour.
#'
#' @export

filledContour2 <- function(z, x=NULL, y=NULL, l=11, col=bobRedBlue2(), orientation=1, zmin=min(z,na.rm=TRUE), zmax=max(z,na.rm=TRUE), main=NA, xlab="x", ylab="y", xlab_line=3, ylab_line=3) {

  # rotate z as needed
  if (orientation==2) {
    z <- t(z[nrow(z):1,])
  } else if (orientation==3) {
    z <- z[nrow(z):1,ncol(z):1]
  } else if (orientation==4) {
    z <- t(z[,ncol(z):1])
  }

  # set x and y based on matrix dimensions
  if (is.null(x))
    x <- 1:nrow(z)
  if (is.null(y))
    y <- 1:ncol(z)

  # produce plot
  myLevels <- seq(zmin, zmax, l=l+1)
  myCols <- smoothCols(1:l,rawCols=col)

  parStore <- par(mar=c(1+xlab_line, 1+ylab_line, 3, 1))
  image(x, y, z, zlim=c(zmin, zmax), col=myCols, xlab=NA, ylab=NA, main=main)
  title(xlab=xlab, line=xlab_line)
  title(ylab=ylab, line=ylab_line)
  contour(x, y, z, levels=myLevels, drawlabels=FALSE, add=TRUE)
  par(parStore)
}

# -----------------------------------
#' vec2mat
#'
#' Produce matrix from two vectors. dim=1 returns the x-matrix, dim=2 the y-matrix.
#'
#' @export

vec2mat <- function(x, y, dim) {

  if (dim==1) {
    output <- matrix(rep(x,each=length(y)),length(y))
  } else {
    output <- matrix(rep(y,length(x)),length(y))
  }
  return(output)
}

# -----------------------------------
#' bin2D
#'
#' Read in x and y data, along with vectors of breaks. Output 2D bin counts and vectors of midpoints.
#'
#' @export

bin2D <- function(x, y, x_breaks, y_breaks) {

  # bin data in both x and y
  freq <- as.data.frame(table(findInterval(x,x_breaks),findInterval(y,y_breaks)))
  freq[,1] <- as.numeric(as.character(freq[,1]))
  freq[,2] <- as.numeric(as.character(freq[,2]))

  # get rid of bins outside of range
  freq <- freq[freq[,1]>0 & freq[,1]<length(x_breaks) & freq[,2]>0 & freq[,2]<length(y_breaks),]

  # fill in 2D matrix
  freq2D <- matrix(0,length(y_breaks)-1,length(x_breaks)-1)
  freq2D[cbind(freq[,1],freq[,2])] <- freq[,3]

  # calculate midpoints
  x_mids <- (x_breaks[-1]+x_breaks[-length(x_breaks)])/2
  y_mids <- (y_breaks[-1]+y_breaks[-length(y_breaks)])/2

  # output all as list
  output <- list(x_mids= x_mids,y_mids= y_mids,z=freq2D)

  return(output)
}

# -----------------------------------
#' safeDivide
#'
#' Divide two numbers, but return 0 rather than NaN if both numerator and denominator are zero.
#'
#' @export

safeDivide <- function(a,b) {
  output <- a/b
  output[a==0 & b==0] <- 0
  return(output)
}

# -----------------------------------
#' exit
#'
#' Force-exits R.
#'
#' @export
exit <- function() {
  exit_cpp()
}

# -----------------------------------
#' ribbon
#'
#' Adds ribbon to plot. Set upper and lower limits of ribbon, or middle values and thickness.
#'
#' @param y1 lower limit of ribbon, or alternatively middle value if upperLower=FALSE
#' @param y2 upper limit of ribbon, or alternatively thickness if upperLower=FALSE
#' @param x x-values
#' @param upperLower whether y-values represent upper and lower values of the ribbon
#' @param density density of shading (leave NA for no shading)
#' @param border colour of border (leave NA for no border)
#' @param col colour of ribbon, or colour of shading lines if density>0
#'
#' @export
ribbon <- function(y1, y2, x=1:length(y1), upperLower=TRUE, density=NA, border=NA, col='#FF000020') {

  # choose limits based on method choice
  if (upperLower) {
    y_lower <- y1
    y_upper <- y2
  } else {
    y_lower <- y1-y2/2
    y_upper <- y1+y2/2
  }

  # build poly coordinates
  poly_x <- c(x, rev(x))
  poly_y <- c(y_lower, rev(y_upper))

  # add polygon to plot
  polygon(poly_x, poly_y, density=NA, border=NA, col=col)
}


# -----------------------------------
#' is_number
#'
#' Check that a given string can be interpreted as a number with as.numeric() without returning an NA value.
#'
#' @export
is_number <- function(x) {
  ret <- !is.na(suppressWarnings(as.numeric(x)))
  ret[is.na(x)] <- NA
  ret
}
