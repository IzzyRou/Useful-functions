
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

### incidence per group per day
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

###### Force of Infection
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

#### discrete SI
DiscretizeSI <- function(mu,CV,SItrunc){
  SIdistr <- sapply(0:SItrunc, function(k) DiscrSI(k, mu, CV*mu))
  SIdistr <- SIdistr/sum(SIdistr)
  return(SIdistr)
}

#####continuous colour scale for base plot
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

#par(mar=c(5.1,4.1,4.1,2.1)+1)
#par(mfrow=c(1,2))
#legend.col(colfunc(1000),likelihood[!wh])
#legend.col(colfunc(1000),A[!wh])

#colfunc <- colorRampPalette(c("snow1","snow2","snow3","seagreen","orange","firebrick","red"), space = "rgb",bias=1)#colors
#colfunc <- colorRampPalette(c('yellow','orange','red','brown'),bias=1)

