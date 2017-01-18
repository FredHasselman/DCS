#' @title Initialise It
#' @description Load and/or install R packages
#'
#' @param need    A vector of package names to be loaded. The wrapper functions have a predefinded \code{need} list and can be used as shortcuts (see details).
#' @param inT    Logical. If \code{TRUE} (default), packages in \code{need} wil be installed if they are not available on the system.
#'
#' @details \code{in.IT} will check if the Packages in the list argument \code{need} are installed on the system and load them. If \code{inT=TRUE} (default), it will first install the packages if they are not present and then proceed to load them.
#'
# The wrapper functions have a predefined need list:
#
# \itemize{
# \item \strong{in.IO}: Load I/O and data handling tools, \code{need = c("plyr","reshape2","RCurl","httr","dplyr","rio")}
# \item \strong{in.PLOT}: Load tools for plotting \code{need = c("lattice","latticeExtra","gplots","ggplot2","grid","gridExtra","scales","effects","RColorBrewer")}
# \item \strong{in.NLTS}: Load tools for Nonlinear Time Series Analysis, \code{need = c("fractaldim","fractalrock","RTisean","tsDyn","tseries","tseriesChaos")}
# \item \strong{in.SN}: Load tools for signal analysis (Matlab style) \code{need = c("pracma","signal","EMD","hht","matlab","oce")}
# }
#'
#' @export
#'
#' @author Fred Hasselman
#'
#' @family initialise packages
#'
#' @examples
#' in.IT(c("reshape2", "plyr", "dplyr"))
in.IT <- function(need=NULL,inT=TRUE){
  ip <- .packages(all.available=TRUE)
  if(any((need %in% ip)==FALSE)){
    if(inT==TRUE){
      install.packages(need[!(need %in% ip)])
    } else {
      cat('Package(s):\n',paste(need[(need %in% ip)==FALSE],sep='\n'),'\nnot installed.\nUse in.IT(c("packagename1","packagename2",...),inT=TRUE)')
      need <- need[(need %in% ip)==TRUE]
    }
  }
  ok <- sapply(1:length(need),function(p) require(need[[p]],character.only=TRUE))
}
# #' @title Wrapper for in.IT: Load I/O and data handling tools
# #'
# #' @rdname in.IT
# #'
# #' @export
# #'
# #' @examples
# #' in.IO()
# in.IO <- function(need = c("xlsx","plyr","doBy","reshape2","RCurl","XML","httr","dplyr")){
#     in.IT(need)
# }
# #' @title Wrapper for in.IT: Parallel computing tools
# #'
# #' @rdname in.IT
# #'
# #' @export
# #'
# #' @examples
# #' in.PAR()
# in.PAR <- function(need = c("parallel","doParallel","foreach")){
#     in.IT(need)
# }
# #' @title Wrapper for in.IT: Tools for plotting
# #'
# #' @rdname in.IT
# #'
# #' @export
# #'
# #' @examples
# #' in.PLOT()
# in.PLOT <- function(need = c("lattice","latticeExtra","gplots","ggplot2","grid","gridExtra","scales","effects","RColorBrewer")){
#     in.IT(need)
# }
# #' @title Wrapper for in.IT: Nonlinear Time Series packages
# #'
# #' @rdname in.IT
# #'
# #' @export
# #'
# #' @examples
# #' in.NLTS()
# in.NLTS <- function(need = c("fractaldim","fractalrock","RTisean","tsDyn","tseries","tseriesChaos")){
#     in.IT(need)
# }
# #' @title Wrapper for in.IT: Signal analysis packages
# #'
# #' @rdname in.IT
# #'
# #' @export
# #'
# #' @examples
# #' in.SN()
# in.SIGN <- function(need=c("pracma","signal","EMD","hht","matlab","oce")){
#     in.IT(need)
# }

#' @title Un-initialise It
#' @description Unload and/or uninstall R packages.
#' @param loose    A vector of package names to be unloaded.
#' @param unT    Logical. If \code{TRUE}, packages in \code{loose} wil be un-installed if they are available on the system.
#'
#' @details \code{un.IT}will check if the Packages in the list argument \code{loose} are installed on the system and unload them. If \code{unT=TRUE} it will first unload the packages if they are loaded, and then proceed to uninstall them.
#'
#' @export
#'
#' @author Fred Hasselman
#'
#' @family initialise packages
#'
#' @examples
#' \dontrun{un.IT(loose = c("reshape2", "plyr", "dplyr"), unT = FALSE)}
un.IT <- function(loose,unT=FALSE){
  dp <- .packages()
  if(any(loose %in% dp)){
    for(looseLib in loose[(loose %in% dp)]){detach(paste0("package:",looseLib), unload=TRUE,character.only=TRUE)}
  }
  rm(dp)
  if(unT==TRUE){
    dp <- .packages(all.available=TRUE)
    if(any(loose %in% dp)){remove.packages(loose[(loose %in% dp)])}
  }
}

#' Autocatlytic Growth: Iterating differential equations (maps)
#'
#' @title Autocatlytic Growth: Iterating differential equations (maps)
#'
#' @param Y0    Initial value.
#' @param r    Growth rate parameter.
#' @param k    Carrying capacity.
#' @param N    Length of the time series.
#' @param type    One of: "driving" (default), "damping", "logistic", "vanGeert1991".
#'
#' @return A timeseries object of length N.
#' @export
#'
#' @author Fred Hasselman
#'
#' @family autocatalytic growth functions
#'
#' @examples
#' # The logistic map in the chaotic regime
#' growth.ac(Y0 = 0.01, r = 4, type = "logistic")
growth.ac <- function(Y0 = 0.01, r = 1, k = 1, N = 100, type = c("driving", "damping", "logistic", "vanGeert")[1]){
  # Create a vector Y of length N, which has value Y0 at Y[1]
  if(N>1){
    Y <- as.numeric(c(Y0, rep(NA,N-2)))
    # Conditional on the value of type ...
    switch(type,
           # Iterate N steps of the difference function with values passed for Y0, k and r.
           driving  = sapply(seq_along(Y), function(t) Y[[t+1]] <<- r * Y[t] ),
           damping  = k + sapply(seq_along(Y), function(t) Y[[t+1]] <<- - r * Y[t]^2 / k),
           logistic = sapply(seq_along(Y), function(t) Y[[t+1]] <<- r * Y[t] * ((k - Y[t]) / k)),
           vanGeert = sapply(seq_along(Y), function(t) Y[[t+1]] <<- Y[t] * (1 + r - r * Y[t] / k))
    )}
  return(ts(Y))
}

#' Conditional Autocatlytic Growth: Iterating differential equations (maps)
#'
#' @param Y0 Initial value
#' @param r Growth rate parameter
#' @param k Carrying capacity
#' @param cond Conditional rules passed as a data.frame of the form: cbind.data.frame(Y = ..., par = ..., val = ...)
#' @param N Length of the time series
#'
#' @export
#'
#' @author Fred Hasselman
#'
#' @family autocatalytic growth functions
#'
#' @examples
#' # Plot with the default settings
#' xyplot(growth.ac.cond())
#'
#' # The function such that it can take a set of conditional rules and apply them sequentially during the iterations.
#' # The conditional rules are passed as a `data.frame`
#' (cond <- cbind.data.frame(Y = c(0.2, 0.6), par = c("r", "r"), val = c(0.5, 0.1)))
#' xyplot(growth.ac.cond(cond=cond))
#'
#' # Combine a change of `r` and a change of `k`
#' (cond <- cbind.data.frame(Y = c(0.2, 1.99), par = c("r", "k"), val = c(0.5, 3)))
#' xyplot(growth.ac.cond(cond=cond))
#'
#' # A fantasy growth process
#' (cond <- cbind.data.frame(Y = c(0.1, 1.99, 1.999, 2.5, 2.9), par = c("r", "k", "r", "r","k"), val = c(0.3, 3, 0.9, 0.1, 1.3)))
#' xyplot(growth.ac.cond(cond=cond))
growth.ac.cond <- function(Y0 = 0.01, r = 0.1, k = 2, cond = cbind.data.frame(Y = 0.2, par = "r", val = 2), N = 100){
  # Create a vector Y of length N, which has value Y0 at Y[1]
  Y <- c(Y0, rep(NA, N-1))
  # Iterate N steps of the difference equation with values passed for Y0, k and r.
  cnt <- 1
  for(t in seq_along(Y)){
    # Check if the current value of Y is greater than the threshold for the current conditional rule in cond
    if(Y[t] > cond$Y[cnt]){
      # If the threshold is surpassed, change the parameter settings by evaluating: cond$par = cond$val
      eval(parse(text = paste(cond$par[cnt], "=", cond$val[cnt])))
      # Update the counter if there is another conditional rule in cond
      if(cnt < nrow(cond)){cnt <- cnt + 1}
    }
    # Van Geert growth model
    Y[[t+1]] <- Y[t] * (1 + r - r * Y[t] / k)
  }
  return(ts(Y))
}

#' plotRP
#'
#' @param crqaOutput    List output from \code{\link[crqa]{crqa}}
#'
#' @return
#' @export
#' @author Fred Hasselman
#' @description Creates a recurrence plot from the sparse matrix output generated by \code{\link[crqa]{crqa}}.
#' @examples
plotRP.crqa <- function(crqaOutput){
  require(Matrix)
  require(lattice)
  require(grid)
  require(gridExtra)

  AUTO <- FALSE

  # Is is crqa output?
  if(all(names(crqaOutput)%in%c("RR","DET","NRLINE","maxL","L","ENTR","rENTR","LAM","TT","RP"))){

    # RP dimensions
    nr <- dim(crqaOutput$RP)[1]
    nc <- dim(crqaOutput$RP)[2]

    RP <- crqaOutput$RP

    # AUTO or CROSS?
    if(class(RP)=="dtCMatrix"){
      message("\nRP in crqa output is a Triangular Sparse Matrix, this implies auto-recurrence...\n")
      AUTO <- TRUE
    }
    if(isSymmetric(RP)){
      message("\nRP in crqa output is a Symmetric Sparse Matrix, this implies auto-recurrence: \n >> RP Measures will include the Line of Identity\n >> Use crqa() with side = 'upper' or 'lower' to get correct measures.\n")
      AUTO <- TRUE}

    if(nc*nr>10^6){
      message(paste("\nLarge RP with",nc*nr,"elements... \n >> This could take some time to plot!\n"))
    }

    ifelse(AUTO,
           Titles <- list(main="Auto Recurrence Plot", X = expression(Y[t]), Y = expression(Y[t])),
           Titles <- list(main="Auto Recurrence Plot", X = expression(X[t]), Y = expression(Y[t]))
    )

    # # Thresholded or Distance matrix?
    # ifelse(length(unique(Matrix:::as.vector(RP)))>2,
    #        distPallette <- colorRampPalette(c("red", "white", "blue"))( length(unique(as.vector(RP)))),
    #        distPallette <- colorRampPalette(c("white", "black"))(2)
    # )

    distPallette <- colorRampPalette(c("white", "black"))(2)

    lp <- levelplot(Matrix::as.matrix(RP),
                    colorkey = FALSE, region = TRUE, col.regions = distPallette, useRaster = TRUE, aspect = nc/nr,
                    main = Titles$main,
                    xlab = Titles$X,
                    ylab = Titles$Y)

    txt <- grid.text(paste0("\n   RR: ", round(crqaOutput$RR, digits = 1),
                            "\n  DET: ", round(crqaOutput$DET, digits = 1),
                            "\nmeanL: ", round(crqaOutput$maxL, digits = 1),
                            "\n maxL: ", round(crqaOutput$NRLINE, digits = 1),
                            "\n  LAM: ", round(crqaOutput$LAM, digits = 1),
                            "\n   TT: ", round(crqaOutput$TT, digits = 1),
                            "\n  ENT: ", round(crqaOutput$ENT, digits = 1),
                            "\n rENT: ", round(crqaOutput$rENTR, digits = 1)),
                     .02,.7, just = "left", draw = FALSE, gp = gpar(fontfamily = "mono", cex = .7))
    # ,)

    grid.arrange(lp, txt, ncol=2, nrow=1, widths=c(4, 1), heights=c(4))

  } else {
    warning("\nInput is not list output from function crqa().\n")
  }

}

plotRP.fnn <- function(FNNoutput){
  plot(FNNoutput["combined",],type="b",pch=16, cex=2, col="grey80", ylim=c(0,100), xaxt="n",
       xlab = "Embedding Dimension", ylab = "False Nearest Neighbours")
  lines(FNNoutput["atol",],type="b",pch="a",col="grey30", lty=2)
  lines(FNNoutput["rtol",],type="b",pch="r", col="grey30",lty=2)
  Axis(side=1,at=seq_along(FNNoutput[1,]),labels = dimnames(FNNoutput)[[2]])
  legend("topright",c("Combined","atol","rtol"), pch = c(16,97,114), lty = c(1,2,2), col = c("grey80","grey30","grey30"), pt.cex=c(2,1,1))
}

# GRAPH PLOTTING ---------------------------------------------------------------
graph2svg <- function(TDM,pname){

  # Create weighted Term-Term matrix
  tTM <- as.matrix(TDM)
  TTM <- tTM %*% t(tTM)
  TTM <- log1p(TTM)

  g <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
  g <- simplify(g)

  # Remove vertices that were used in the search query
  Vrem <- which(V(g)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
  g <- (g - V(g)$name[Vrem])

  # Set colors and sizes for vertices
  V(g)$degree <- degree(g)
  rev         <- scaleRange(log1p(V(g)$degree))
  rev[rev<=0.3]<-0.3

  V(g)$color       <- rgb(scaleRange(V(g)$degree), 1-scaleRange(V(g)$degree),  0, rev)
  V(g)$size        <- 10*scaleRange(V(g)$degree)
  V(g)$frame.color <- NA

  # set vertex labels and their colors and sizes
  V(g)$label       <- V(g)$name
  V(g)$label.color <- rgb(0, 0, 0, rev)
  V(g)$label.cex   <- scaleRange(V(g)$degree)+.1

  # set edge width and color
  rew <- scaleRange(E(g)$weight)
  rew[rew<=0.3]<-0.3

  E(g)$width <- 2*scaleRange(E(g)$weight)
  E(g)$color <- rgb(.5, .5, 0, rew)
  set.seed(958)

  svg(paste(pname,sep=""),width=8,height=8)
  plot(g, layout=layout.fruchterman.reingold(g))
  dev.off()

  return(g)
}

# Plot vertex neighbourhood
hoodGraph2svg <- function(TDM,Vname,pname){

  # Create weighted Term-Term matrix
  tTM <- as.matrix(TDM)
  TTM <- tTM %*% t(tTM)
  TTM <- log1p(TTM)

  ig <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
  ig <- simplify(ig)

  # Remove vertices that were used in the search query
  Vrem <- which(V(ig)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
  ig <- (ig - V(ig)$name[Vrem])

  # This is a deletion specific for the Neighbourhood graphs
  Vrem <- which(V(ig)$name %in% c("~rdsp~","~imp~","~som~","~bod~","~mlt~"))
  ig   <- ig - V(ig)$name[Vrem]

  idx <- which(V(ig)$name==Vname)
  sg  <- graph.neighborhood(ig, order = 1, nodes=V(ig)[idx], mode = 'all')[[1]]

  # set colors and sizes for vertices
  V(sg)$degree <- degree(sg)

  rev<-scaleRange(log1p(V(sg)$degree))
  rev[rev<=0.3]<-0.3

  V(sg)$color <- rgb(scaleRange(V(sg)$degree), 1-scaleRange(log1p(V(sg)$degree*V(sg)$degree)),  0, rev)

  V(sg)$size        <- 35*scaleRange(V(sg)$degree)
  V(sg)$frame.color <- NA

  # set vertex labels and their colors and sizes
  V(sg)$label       <- V(sg)$name
  V(sg)$label.color <- rgb(0, 0, 0, rev)
  V(sg)$label.cex   <- scaleRange(V(sg)$degree)

  # set edge width and color
  rew<-scaleRange(E(sg)$weight)
  rew[rew<=0.3]<-0.3

  E(sg)$width <- 6*scaleRange(E(sg)$weight)
  E(sg)$color <- rgb(.5, .5, 0, rew)

  idV <- which(V(sg)$name==Vname)
  idE <- incident(sg,V(sg)[[idV]])
  E(sg)$color[idE] <- rgb(0, 0, 1 ,0.8)

  set.seed(958)

  idx <- which(V(sg)$name==Vname)
  svg(paste(pname,sep=""),width=8,height=8)
  plot(sg,layout=layout.star(sg,center=V(sg)[idx]))
  dev.off()

  return(sg)
}

sliceTS<-function(TSmat,epochSz=1) {
  # Slice columns of TSmat in epochs of size = epochSz
  require(plyr)

  N<-dim(TSmat)
  return(llply(seq(1,N[1],epochSz),function(i) TSmat[i:min(i+epochSz-1,N[1]),1:N[2]]))
}

fltrIT <- function(TS,f){
  # Apply filtfilt to TS using f (filter settings)
  require("signal")

  return(filtfilt(f=f,x=TS))

}

SWtest0 <- function(g){
  Nreps <- 10;
  histr  <- vector("integer",Nreps)
  target<- round(mean(degree(g)))
  now   <- target/2
  for(i in 1:Nreps){
    gt      <- watts.strogatz.game(dim=1, size=length(degree(g)), nei=now, 0)
    histr[i] <- round(mean(degree(gt)))
    ifelse(histr[i] %in% histr,break,{
      ifelse(histr[i]>target,{now<-now-1},{
        ifelse(histr[i]<target,{now<-now+1},{
          break})
      })
    })
  }
  return(gt)
}


# SWtestV <- function(g,N){
#  return(list(cp=transitivity(g,type="global"),cpR=transitivity(rewire(g,mode=c("simple"),niter=N),type="global"),lp=average.path.length(g), lpR=average.path.length(rewire(g,mode=c("simple"),niter=N))))
# }

SWtestE <- function(g,p=1,N=20){
  values <- matrix(nrow=N,ncol=6,dimnames=list(c(1:N),c("cp","cpR","cp0","lp","lpR","lp0")))

  for(n in 1:N) {
    gt<-SWtest0(g)
    values[n,] <- c(transitivity(g,type="localaverage"),transitivity(rewire(g,each_edge(p=p)),type="localaverage"),transitivity(gt,type="localaverage"),average.path.length(g),average.path.length(rewire(g,each_edge(p=p))),average.path.length(gt))}
  values[n,values[n,]==0] <- NA #values[n,values[n,]==0]+1e-8}

  values   <- cbind(values,(values[,1]/values[,2])/(values[,4]/values[,5]),(values[,1]/values[,3]),(values[,4]/values[,6]),((values[,1]/values[,3])/values[,2])/((values[,4]/values[,6])/values[,5]))
  valuesSD <- data.frame(matrix(apply(values[,1:10],2,sd,na.rm=T),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  valuesAV <- data.frame(matrix(colMeans(values[,1:10],na.rm=T),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  return(list(valuesAV=valuesAV,valuesSD=valuesSD,valuesSE=valuesSD/sqrt(N)))
}

PLFsmall <- function(g){
  reload <- FALSE
  if("signal" %in% .packages()){
    warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
    reload <- TRUE
    detach("package:signal", unload=TRUE)}

  if(length(V(g))>100){warning("Vertices > 100, no need to use PLFsmall, use a binning procedure");break}
  d <- degree(g,mode="all")

  y <- hist(d,breaks=0.5:(max(d)+0.5),plot=FALSE)$counts
  #y<-y[y>0]
  if(length(y)<2){
    warning("Less than 2 points in Log-Log regression... alpha=0")
    alpha <- 0
  } else {
    if(length(y)==2){
      warning("Caution... Log-Log slope is a bridge (2 points)")
      chop <- 0
    } else {
      chop <- 1
    }
    alpha <- coef(lm(rev(log1p(y)[1:(length(y)-chop)]) ~ log1p(1:(length(y)-chop))))[2]
  }
  if(reload==TRUE){library(signal,verbose=FALSE,quietly=TRUE)}

  return(alpha)
}

plotSW <- function(n,k,p){

  g <- watts.strogatz.game(1, n, k, p)

  V(g)$degree <- degree(g)

  # set colors and sizes for vertices
  rev<-scale.R(log1p(V(g)$degree))
  rev[rev<=0.2]<-0.2
  rev[rev>=0.9]<-0.9
  V(g)$rev <- rev$x

  V(g)$color       <- rgb(V(g)$rev, 1-V(g)$rev,  0, 1)
  V(g)$size        <- 25*V(g)$rev

  # set vertex labels and their colors and sizes
  V(g)$label       <- ""

  E(g)$width <- 1
  E(g)$color <- rgb(0.5, 0.5, 0.5, 1)

  return(g)
}

plotBA <- function(n,pwr,out.dist){
  #require("Cairo")

  g <- barabasi.game(n,pwr,out.dist=out.dist,directed=F)
  V(g)$degree <- degree(g)

  # set colors and sizes for vertices
  rev<-scale.R(log1p(V(g)$degree))
  rev[rev<=0.2] <- 0.2
  rev[rev>=0.9] <- 0.9
  V(g)$rev <- rev$x

  V(g)$color    <- rgb(V(g)$rev, 1-V(g)$rev,  0, 1)
  V(g)$size     <- 25*V(g)$rev
  # V(g)$frame.color <- rgb(.5, .5,  0, .4)

  # set vertex labels and their colors and sizes
  V(g)$label <- ""

  E(g)$width <- 1
  E(g)$color <- rgb(0.5, 0.5, 0.5, 1)

  return(g)
}


FDrel <- function(g){
  d<-degree(g,mode="all")
  nbreaks <- round(length(V(g))/2)-1
  y<-hist(d,breaks=nbreaks,plot=F)$density
  y<-y[y>0]
  return(FD <- -sum(y*log2(y))/-(log2(1/length(y))))
}

sa2fd <- function(sa, ...) UseMethod("sa2fd")

sa2fd.default <- function(sa, ...){
  cat("No type specified.")
}


#' Informed Dimension estimate from Spectral Slope (aplha)
#'
#' @description Conversion formula: From periodogram based self-affinity parameter estimate (\code{sa}) to an informed estimate of the (fractal) dimension (FD).
#' @param sa Self-Affinity parameter estimate based on PSD slope (e.g., \code{\link{fd.psd}})).
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#' @export
#'
#' @details The spectral slope will be converted to a dimension estimate using:
#'
#' \deqn{D_{PSD}\approx\frac{3}{2}+\frac{14}{33}*\tanh\left(Slope * \ln(1+\sqrt{2})\right)}
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @family SA to FD converters
#'
#' @examples
#' # Informed FD of white noise
#' sa2fd.psd(0)
#'
#' # Informed FD of Brownian noise
#' sa2fd.psd(-2)
#'
#' # Informed FD of blue noise
#' sa2fd.psd(2)
sa2fd.psd <- function(sa){return(round(3/2 + ((14/33)*tanh(sa*log(1+sqrt(2)))), digits = 2))}


#' Informed Dimension estimate from DFA slope (H)
#'
#' @description Conversion formula: Detrended Fluctuation Analysis (DFA) estimate of the Hurst exponent (a self-affinity parameter \code{sa}) to an informed estimate of the (fractal) dimension (FD).
#'
#' @param sa Self-Afinity parameter estimate based on DFA slope (e.g., \code{\link{fd.sda}})).
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#'
#' @export
#'
#' @details The DFA slope (H) will be converted to a dimension estimate using:
#'
#' \deqn{D_{DFA}\approx 2-(\tanh(\log(3)*sa)) }{D_{DFA} ≈ 2-(tanh(log(3)*sa)) }
#'
#' @family SA to FD converters
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @examples
#' # Informed FD of white noise
#' sa2fd.dfa(0.5)
#'
#' # Informed FD of Pink noise
#' sa2fd.dfa(1)
#'
#' # Informed FD of blue noise
#' sa2fd.dfa(0.1)
sa2fd.dfa <- function(sa){return(round(2-(tanh(log(3)*sa)), digits = 2))}


#' Informed Dimension estimate from SDA slope.
#'
#' @description Conversion formula: Standardised Dispersion Analysis (SDA) estimate of self-affinity parameter (\code{SA}) to an informed estimate of the fractal dimension (FD).
#'
#' @param sa Self-afinity parameter estimate based on SDA slope (e.g., \code{\link{fd.sda}})).
#'
#' @details
#'
#'
#' Note that for some signals different PSD slope values project to a single SDA slope. That is, SDA cannot distinguish between all variaties of power-law scaling in the frequency domain.
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#' @export
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @family SA to FD converters
#'
#' @examples
#' # Informed FD of white noise
#' sa2fd.sda(-0.5)
#'
#' # Informed FD of Brownian noise
#' sa2fd.sda(-1)
#'
#' # Informed FD of blue noise
#' sa2fd.sda(-0.9)
sa2fd.sda <- function(sa){return(1-sa)}



# fd estimators ---------------------------------------------------------------------------------------------------

fd <- function(y, ...) UseMethod("fd")

fd.default <- function(y, ...){

  cat("No type specified.\nReturning exponential growth power law.")

  r = 1.01
  y <- growth.ac(Y0=0.001, r=r, N=2048, type = "driving")
  tsp(y) <-c(1/500,2048/500,500)
  bulk <- log1p(hist(y,plot = F, breaks = seq(0,max(y),length.out = 129))$counts)
  size <- log1p(seq(0,2047,length.out = 128))
  id<-bulk==0

  lmfit <- lm(bulk[!id] ~ size[!id])

  old <- ifultools::splitplot(2,1,1)
  plot(y, ylab = "Y", main = paste0('Exponential growth  sap: ', round(coef(lmfit)[2],digits=2), ' | r:', r))
  ifultools::splitplot(2,1,2)
  plot(size[!id],bulk[!id], xlab="Size = log(bin(Time))", ylab = "Bulk = logbin(Y)", pch=21, bg="grey60", pty="s")
  lines(size[!id], predict(lmfit),lwd=4,col="darkred")
  #legend("bottomleft",c(paste0("Range (n = ",sum(powspec$size<=0.25),")"), paste0("Hurvic-Deo estimate (n = ",nr,")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
  par(old)
}


# PSD -------------------------------------------------------------------------------------------------------------

#' @title Power Spectral Density Slope (PSD).

#' @description Estimate Alpha, Hurst Exponent and Fractal Dimension through log-log slope.
#'
#' @param y    A numeric vector or time series object.
#' @param normalize    Normalize the series (default).
#' @param detrend    Subtract linear trend from the series (default).
#' @param plot    Return the log-log spectrum with linear fit (default).
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @return A list object containing:
#' \itemize{
#' \item A data matrix \code{PLAW} with columns \code{freq.norm}, \code{size} and \code{bulk}.
#' \item Estimate of scaling exponent \code{alpha} based on a fit over the lowest 25\% frequencies (\code{low25}), or using the HD estimate \code{HD}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @family FD estimators
#'
#' @export
#'
#' @details Calls function \code{\link[sapa]{SDF}} to estimate the scaling exponent of a timeseries based on the periodogram frequency spectrum. After detrending and normalizing the signal (if requested), \code{SDF} is called using a Tukey window (\code{raised cosine \link[sapa]{taper}}).
#'
#' A line is fitted on the periodogram in log-log coordinates. Two fit-ranges are used: The 25\% lowest frequencies and the Hurvich-Deo estimate (\code{\link[fractal]{HDEst}}).
#'
#' @examples
#' fd.psd(rnorm(2048), plot = TRUE)
fd.psd <- function(y, fs = NULL, normalize = TRUE, dtrend = TRUE, plot = FALSE){
  require(pracma)
  require(fractal)
  require(sapa)
  require(ifultools)

  if(!is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- ts(y, frequency = fs)
    cat("\n\nfd.psd:\tSample rate was set to 1.\n\n")
  }

  N             <- length(y)
  # Simple linear detrending.
  if(dtrend)    y <- ts(pracma::detrend(as.vector(y), tt = 'linear'), frequency = fs)
  # Normalize using N instead of N-1.
  if(normalize) y <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE)*sqrt((N-1)/N))

  # Number of frequencies estimated cannot be set! (defaults to Nyquist)
  # Use Tukey window: cosine taper with r = 0.5

  # fast = TRUE ensures padding with zeros to optimize FFT to highly composite number.
  # However, we just pad to nextPow2, except if length already is a power of 2.
  npad <- 1+(stats::nextn(N,factors=2)-N)/N
  npad <- stats::nextn(N)

  # if(N==npad) npad = 0
  # psd  <- stats::spec.pgram(y, fast = FALSE, demean=FALSE, detrend=FALSE, plot=FALSE, pad=npad, taper=0.5)

  Tukey <- sapa::taper(type="raised cosine", flatness = 0.5, n.sample = npad)
  psd   <- sapa::SDF(y, taper. = Tukey, npad = npad)

  powspec <- cbind.data.frame(freq.norm = attr(psd, "frequency")[-1], size = attr(psd, "frequency")[-1]*frequency(y), bulk = as.matrix(psd)[-1])

  # First check the global slope for anti-persistent noise (GT +0.20)
  # If so, fit the line starting from the highest frequency
  nr     <- length(powspec[,1])
  lsfit  <- lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
  glob   <- coef(lsfit)[2]

  # General guideline: fit over 25% frequencies
  # If signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure
  nr <- fractal::HDEst(NFT = length(powspec[,1]), sdf = psd)

  exp1 <- fractal::hurstSpec(y, sdf.method = "direct", freq.max = 0.25, taper. = Tukey )
  exp2 <- fractal::hurstSpec(y, sdf.method = "direct", freq.max = powspec$freq.norm[nr], taper. = Tukey)

  ifelse((glob > 0.2), {
    lmfit1 <- lm(log(rev(powspec$bulk[powspec$size<=0.25])) ~ log(rev(powspec$size[powspec$size<=0.25])))
    lmfit2 <- lm(log(rev(powspec$bulk[1:nr])) ~ log(rev(powspec$size[1:nr])))
  },{
    lmfit1 <- lm(log(powspec$bulk[powspec$size<=0.25]) ~ log(powspec$size[powspec$size<=0.25]))
    lmfit2 <- lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
  })

  if(plot){
    old<- ifultools::splitplot(2,1,1)
    plot(y,ylab = "Y", main = paste0('Lowest 25%    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:', round(exp1,digits=2), ' | FD:',round(sa2fd.psd(coef(lmfit1)[2]),digits=2),'\nHurvic-Deo    sap: ', round(coef(lmfit2)[2],digits=2), ' | H:', round(exp2,digits=2), ' | FD:',round(sa2fd.psd(coef(lmfit2)[2]),digits=2)))
    ifultools::splitplot(2,1,2)
    plot(log(powspec$bulk) ~ log(powspec$size), xlab="log(Frequency)", ylab = "log(Power)")
    lines(log(powspec$size[powspec$size<=0.25]), predict(lmfit1),lwd=3,col="darkred")
    lines(log(powspec$size[1:nr]), predict(lmfit2),lwd=3,col="darkblue")
    legend("bottomleft",c(paste0("lowest 25% (n = ",sum(powspec$size<=0.25),")"), paste0("Hurvic-Deo estimate (n = ",nr,")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
    par(old)
  }

  return(list(
    PLAW  = powspec,
    low25 = list(sap = coef(lmfit1)[2], H = exp1, FD = sa2fd.psd(coef(lmfit1)[2]), fitlm1 = lmfit1),
    HD    = list(sap = coef(lmfit2)[2], H = exp2, FD = sa2fd.psd(coef(lmfit2)[2]), fitlm2 = lmfit2),
    info  = psd)
  )
}


# SDA -------------------------------------------------

#' fd.sda
#'
#' @title Standardised Dispersion Analysis (SDA).
#'
#' @param y    A numeric vector or time series object.
#' @param normalize    Normalize the series (default).
#' @param plot    Return the log-log spectrum with linear fit (default).
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @return A list object containing:
#' \itemize{
#' \item A data matrix \code{PLAW} with columns \code{freq.norm}, \code{size} and \code{bulk}.
#' \item Estimate of scaling exponent \code{sap} based on a fit over the standard range (\code{fullRange}), or on a user defined range \code{fitRange}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @export
#'
#' @family FD estimators
#' @examples
fd.sda <- function(y, fs = NULL, normalize = TRUE, dtrend = FALSE, scales = dispersion(y)$scale, fitRange = c(scales[1], scales[length(scales)-2]), plot = FALSE){
  require(pracma)
  require(fractal)

  if(!is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- ts(y, frequency = fs)
    cat("\n\nfd.sda:\tSample rate was set to 1.\n\n")
  }

  N             <- length(y)
  # Simple linear detrending.
  if(dtrend)    y <- ts(pracma::detrend(as.vector(y), tt = 'linear'), frequency = fs)
  # Normalize using N instead of N-1.
  if(normalize) y <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE)*sqrt((N-1)/N))

  bins          <- which(fitRange[1]==scales):which(fitRange[2]==scales)
  out           <- dispersion(y, front = FALSE)
  lmfit1        <- lm(log(out$sd) ~ log(out$scale))
  lmfit2        <- lm(log(out$sd[bins]) ~ log(out$scale[bins]))

  if(plot){
    old<- ifultools::splitplot(2,1,1)
    plot(y,ylab = "Y", main = paste0('Full    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:', round(1+coef(lmfit1)[2],digits=2), ' | FD:',round(sa2fd.sda(coef(lmfit1)[2]),digits=2),'\nRange    sap: ', round(coef(lmfit2)[2],digits=2), ' | H:', round(1+coef(lmfit1)[2],digits=2), ' | FD:',round(sa2fd.sda(coef(lmfit2)[2]),digits=2)))
    ifultools::splitplot(2,1,2)
    plot(log(out$sd) ~ log(out$scale), xlab="log(Bin Size)", ylab = "log(SD)")
    lines(log(out$scale), predict(lmfit1),lwd=3,col="darkred")
    lines(log(out$scale[bins]), predict(lmfit2),lwd=3,col="darkblue")
    legend("bottomleft",c(paste0("Full (n = ",length(out$scale),")"), paste0("Range (n = ",length(bins),")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
    par(old)
  }

  return(list(
    PLAW  =  cbind.data.frame(freq.norm = frequency(y)/scales, size = out$scale, bulk = out$sd),
    fullRange = list(sap = coef(lmfit1)[2], H = 1+coef(lmfit1)[2], FD = sa2fd.sda(coef(lmfit1)[2]), fitlm1 = lmfit1),
    fitRange  = list(sap = coef(lmfit2)[2], H = 1+coef(lmfit2)[2], FD = sa2fd.sda(coef(lmfit2)[2]), fitlm2 = lmfit2),
    info = out)
  )
}


# DFA ---------------------------------------------

#' fd.dfa
#'
#' @title Detrended Fluctuation Analysis (DFA)
#'
#' @param y    A numeric vector or time series object.
#' @param normalize    Normalize the series (default).
#' @param detrend    Subtract linear trend from the series (default).
#' @param dmethod     Method to use for detrending, see \code{\link[fractal]{DFA}}.
#' @param plot    Return the log-log spectrum with linear fit (default).
#'
#'
#' @return Estimate of Hurst exponent (slope of \code{log(bin)} vs. \code{log(RMSE))} and an FD estimate based on Hasselman(2013)
#' A list object containing:
#' \itemize{
#' \item A data matrix \code{PLAW} with columns \code{freq.norm}, \code{size} and \code{bulk}.
#' \item Estimate of scaling exponent \code{sap} based on a fit over the standard range (\code{fullRange}), or on a user defined range \code{fitRange}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @export
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @family FD estimators
#' @examples
fd.dfa <- function(y, fs = NULL, dtrend = "poly1", normalize = FALSE, sum.order = 1, scale.max=trunc(length(y)/4), scale.min=4, scale.ratio=2^(1/4), overlap = 0, plot = FALSE){
  require(pracma)
  require(fractal)

  reload <- FALSE
  if("signal" %in% .packages()){
    warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
    reload <- TRUE
    detach("package:signal", unload=TRUE)
  }



  if(!is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- ts(y, frequency = fs)
    cat("\n\nfd.dfa:\tSample rate was set to 1.\n\n")
  }

  N             <- length(y)
  # Normalize using N instead of N-1.
  if(normalize) y <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE)*sqrt((N-1)/N))

  out1 <- fractal::DFA(y, detrend=dtrend, sum.order=sum.order, scale.max=trunc(length(y)/2), scale.min=2, scale.ratio=2, overlap = 0, verbose=FALSE)
  out2 <- fractal::DFA(y, detrend=dtrend, sum.order=sum.order, scale.max=scale.max, scale.min=scale.min, scale.ratio=scale.ratio, overlap = overlap, verbose=FALSE)

  lmfit1        <- lm(log(attributes(out1)$stat) ~ log(attributes(out1)$scale))
  lmfit2        <- lm(log(attributes(out2)$stat) ~ log(attributes(out2)$scale))

  if(plot){
    plot.new()
    old <- ifultools::splitplot(2,1,1)
    plot(y,ylab = "Y", main = paste0('Full    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:',
                                     round(attributes(out1)$logfit[]$coefficients['x'] ,digits=2), ' | FD:',
                                     round(sa2fd.dfa(coef(lmfit1)[2]),digits=2),'\nRange    sap: ',
                                     round(coef(lmfit2)[2],digits=2), ' | H:',
                                     round(attributes(out2)$logfit[]$coefficients['x'] ,digits=2), ' | FD:',
                                     round(sa2fd.dfa(coef(lmfit2)[2]),digits=2)
    )
    )
    ifultools::splitplot(2,1,2)
    plot(log(attributes(out1)$stat) ~ log(attributes(out1)$scale), xlab="log(Bin Size)", ylab = "log(RMSE)")
    lines(log(attributes(out1)$scale), predict(lmfit1),lwd=3,col="darkred")
    lines(log(attributes(out2)$scale), predict(lmfit2),lwd=3,col="darkblue")
    legend("topleft",c(paste0("Full (n = ",length(attributes(out1)$scale),")"), paste0("Range (n = ",length(attributes(out2)$scale),")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
    par(old)
  }

  if(reload==TRUE){library(signal,verbose=FALSE,quietly=TRUE)}

  return(list(
    PLAW  =  cbind.data.frame(freq.norm = scale.R(attributes(out1)$scale*frequency(y)), size = attributes(out1)$scale, bulk = attributes(out1)$stat),
    fullRange = list(sap = coef(lmfit1)[2], H = attributes(out1)$logfit[]$coefficients['x'] , FD = sa2fd.dfa(coef(lmfit1)[2]), fitlm1 = lmfit1),
    fitRange  = list(sap = coef(lmfit2)[2], H = coef(lmfit2)[2], FD = sa2fd.dfa(coef(lmfit2)[2]), fitlm2 = lmfit2),
    info = list(out1,out2))
  )
}



#' Detrended Fluctuation Analysis
#'
#' @param signal    An input signal.
#' @param qq    A vector containing a range of values for the order of fluctuation \code{q}.
#' @param mins    Minimum scale to consider.
#' @param maxs    Maximum scale to consider.
#' @param ressc
#' @param m
#'
#' @return
#' @export
#'
#' @examples
DFA1 <- function(signal,mins=4,maxs=12,ressc=30,m=1){
  require(pracma)
  require(plyr)
  require(dplyr)

  #   reload <- FALSE
  #   if("signal" %in% .packages()){
  #     warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
  #     reload <- TRUE
  #     detach("package:signal", unload=TRUE)
  #   }
  scale     <- round(2^(seq(mins,maxs,by=((maxs-mins)/ressc))))
  segv      <- numeric(length(scale))
  RMS_scale <- vector("list",length(scale))
  # qRMS      <- vector("list",length(qq))
  # Fq        <- vector("list",length(qq))
  # qRegLine  <- vector("list",length(qq))
  # Hq        <- numeric(length(qq))

  Y        <- cumsum(signal-mean(signal))
  TSm      <- as.matrix(cbind(t=1:length(Y),y=Y))
  Hglobal  <- monoH(TSm,scale)
}


# Multi-Fractal DFA -----------------------------------------------------------------------------------------------

#' Multi-fractal Detrended Fluctuation Analysis
#'
#' @param signal    An input signal.
#' @param qq    A vector containing a range of values for the order of fluctuation \code{q}.
#' @param mins    Minimum scale to consider.
#' @param maxs    Maximum scale to consider.
#' @param ressc
#' @param m
#'
#' @return
#' @export
#'
#' @examples
MFDFA <- function(signal,qq=c(-10,-5:5,10),mins=6,maxs=12,ressc=30,m=1){
  require(pracma)
  require(plyr)
  require(dplyr)

  #   reload <- FALSE
  #   if("signal" %in% .packages()){
  #     warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
  #     reload <- TRUE
  #     detach("package:signal", unload=TRUE)
  #   }
  scale     <- round(2^(seq(mins,maxs,by=((maxs-mins)/ressc))))
  segv      <- numeric(length(scale))
  RMS_scale <- vector("list",length(scale))
  qRMS      <- vector("list",length(qq))
  Fq        <- vector("list",length(qq))
  qRegLine  <- vector("list",length(qq))
  Hq        <- numeric(length(qq))

  Y        <- cumsum(signal-mean(signal))
  TSm      <- as.matrix(cbind(t=1:length(Y),y=Y))
  Hglobal  <- monoH(TSm,scale)

  Hadj <- 0
  if((Hglobal>1.2)&(Hglobal<1.8)){
    Y <- diff(signal)
    Hadj=1}
  if(Hglobal>1.8){
    Y <- diff(diff(signal))
    Hadj <- 2}
  if(Hglobal<0.2){
    Y <- cumsum(signal-mean(signal))
    Hadj <- -1}
  if(Hadj!=0){TSm  <- as.matrix(cbind(t=1:length(Y),y=cumsum(Y-mean(Y))))}

  for(ns in seq_along(scale)){
    RMS_scale[[ns]] <- ldply(sliceTS(TSm,scale[ns]),function(sv){return(sqrt(mean(detRend(sv[,2]))^2))})
    for(nq in seq_along(qq)){
      qRMS[[nq]][1:length(RMS_scale[[ns]]$V1)] <- RMS_scale[[ns]]$V1^qq[nq]
      Fq[[nq]][ns] <- mean(qRMS[[nq]][1:length(RMS_scale[[ns]]$V1)])^(1/qq[nq])
      if(is.inf(log2(Fq[[nq]][ns]))){Fq[[nq]][ns]<-NA}
    }
    Fq[[which(qq==0)]][ns] <- exp(0.5*mean(log(RMS_scale[[ns]]^2)))
    if(is.inf(log2(Fq[[which(qq==0)]][ns]))){Fq[[which(qq==0)]][ns]<-NA}
  }

  fmin<-1
  fmax<-which(scale==max(scale))
  #for(nq in seq_along(qq)){Hq[nq] <- lm(log2(Fq[[nq]])~log2(scale))$coefficients[2]}
  Hq <- ldply(Fq,function(Fqs){lm(log2(Fqs[fmin:fmax])~log2(scale[fmin:fmax]),na.action=na.omit)$coefficients[2]})

  tq <- (Hq[,1]*qq)-1
  hq <- diff(tq)/diff(qq)
  Dq <- (qq[1:(length(qq)-1)]*hq) - (tq[1:(length(qq)-1)])

  if(reload==TRUE){library(signal,verbose=FALSE,quietly=TRUE)}

  return(list(q=qq,Hq=Hq,tq=tq,hq=hq,Dq=Dq,Hglobal=Hglobal,Hadj=Hadj))
}


monoH <- function(TSm,scale){
  dfaRMS_scale <- vector("list",length(scale))
  F2 <- numeric(length(scale))
  for(ns in seq_along(scale)){
    dfaRMS_scale[[ns]] <- ldply(sliceTS(TSm,scale[ns]),function(sv){return(sqrt(mean(detRend(sv[,2]))^2))})
    F2[ns] <- mean(dfaRMS_scale[[ns]]$V1^2)^(1/2)
    if(is.inf(log2(F2[ns]))){F2[ns] <- NA}
  }
  return(lm(log2(F2)~log2(scale),na.action=na.omit)$coefficients[2])
}

detRend <- function(TS, Order=1){
  detR <- lm(TS~stats::poly(1:length(TS), degree=Order))$residuals
  return(detR)
}


#
# set.seed(100)
# z <- dispersion(rnorm(1024))
# plot(log(z$scale),log(z$sd))
# #

# trace(detRend,edit=T)
# seq(1,length(X),by=4096)
#
# z<-sliceTS(TSm,scale[1])
# z[[1]][,2]
#
# Hglobal <-
#
# segments <- laply(scale,function(s) floor(length(X)/s))
# IDv <- llply(segments,slice.index,)
# segv <- function(X,segments){
#
#
#
#   seq((((v-1)*scale[ns])+1),(v*scale[ns]),length=scale[ns])
#
# }
# for(ns in seq_along(scale)){
#   segv[ns] <- floor(length(X)/scale[ns])
#   for(v in 1:segv[ns]){
#     Index <- seq((((v-1)*scale[ns])+1),(v*scale[ns]),length=scale[ns])
#     Cslope<- polyfit(Index,X[Index],m)
#     fit   <- polyval(Cslope,Index)
#     RMS_scale[[ns]][v] <- sqrt(mean((X[Index]-fit)^2))
#     rm(Cslope, fit, Index)
#   }
#   for(nq in seq_along(qq)){
#     qRMS[[nq]][1:segv[ns]] <- RMS_scale[[ns]]^qq[nq]
#     Fq[[nq]][ns] <- mean(qRMS[[nq]][1:segv[ns]])^(1/qq[nq])
#   }
#   Fq[[which(qq==0)]][ns] <- exp(0.5*mean(log(RMS_scale[[ns]]^2)))
# }
#
# for(nq in seq_along(qq)){
#   Cslope <- polyfit(log2(scale),log2(Fq[[nq]]),1)
#   Hq[nq] <- Cslope[1]
#   qRegLine[[nq]] <- polyval(Cslope,log2(scale))
#   rm(Cslope)
# }
#
# tq <- (Hq*qq)-1
# hq <- diff(tq)/diff(qq)
# Dq <- (qq[1:(length(qq)-1)]*hq) - (tq[1:(length(qq)-1)])
#
# plot(hq,Dq,type="l")
#
#
# qq<-c(-10,-5,seq(-2,2,.1),5,10)


# PLOTS -------------------------------------------------------------------
#
#
#' gg.theme
#'
#' @param type      One of \code{"clean"}, or \code{"noax"}
#' @param useArial    Use the Arial font (requires \code{.afm} font files in the \code{afmPath})
#' @param afmPATH    Path to Arial \code{.afm} font files.
#'
#' @details Will generate a \code{"clean"} ggplot theme, or a theme without any axes (\code{"noax"}).
#'
#' Some scientific journals explicitly request the Arial font should be used in figures. This can be achieved by using \code{.afm} font format (see, e.g. http://www.pure-mac.com/font.html).
#'
#' @return A theme for \code{ggplot2}.
#' @export
#'
#' @examples
#' library(ggplot2)
#' g <- ggplot(data.frame(x = rnorm(n = 100), y = rnorm(n = 100)), aes(x = x, y = y)) + geom_point()
#' g + gg.theme()
#' g + gg.theme("noax")
gg.theme <- function(type=c("clean","noax"),useArial = F, afmPATH="~/Dropbox"){
  require(ggplot2)

  if(length(type)>1){type <- type[1]}

  if(useArial){
    set.Arial(afmPATH)
    bf_font="Arial"
  } else {bf_font="Helvetica"}

  switch(type,
         clean = theme_bw(base_size = 16, base_family=bf_font) +
           theme(axis.text.x     = element_text(size = 14),
                 axis.title.y    = element_text(vjust = +1.5),
                 panel.grid.major  = element_blank(),
                 panel.grid.minor  = element_blank(),
                 legend.background = element_blank(),
                 legend.key = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 axis.line  = element_line(colour = "black")),
         noax = theme(line = element_blank(),
                      text  = element_blank(),
                      title = element_blank(),
                      plot.background = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
  )
}

#' gg.plotHolder
#'
#' @param useArial    Use the Arial font (requires \code{.afm} font files in the \code{afmPath})
#' @param afmPATH    Path to Arial \code{.afm} font files.
#'
#' @return A blank \code{ggplot2} object that can be used in concordance with \code{grid.arrange}.
#' @export
#'
#' @examples
#' # Create a plot with marginal distributions.
#' library(ggplot2)
#' library(scales)
#'
#' df <- data.frame(x = rnorm(n = 100), y = rnorm(n = 100), group = factor(sample(x=c(0,1), size = 100, replace = TRUE)))
#'
#' scatterP <- ggplot(df, aes(x = x, y =y, colour = group)) + geom_point() + gg.theme()
#' xDense <- ggplot(df, aes(x = x, fill = group)) + geom_density(aes(y= ..count..),trim=FALSE, alpha=.5) + gg.theme("noax") + theme(legend.position = "none")
#' yDense <- ggplot(df, aes(x = y, fill = group)) + geom_density(aes(y= ..count..),trim=FALSE, alpha=.5) + coord_flip() + gg.theme("noax") + theme(legend.position = "none")
#'
#' library(gridExtra)
#' grid.arrange(xDense, gg.plotHolder(), scatterP, yDense, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
gg.plotHolder <- function(useArial = F,afmPATH="~/Dropbox"){
  require(ggplot2)
  ggplot() +
    geom_blank(aes(1,1)) +
    theme(line = element_blank(),
          text  = element_blank(),
          title = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    )
}

set.Arial <- function(afmPATH="~/Dropbox"){
  # Set up PDF device on MAC OSX to use Arial as a font in Graphs
  if(nchar(afmPATH>0)){
    if(file.exists(paste0(afmPATH,"/Arial.afm"))){
      Arial <- Type1Font("Arial",
                         c(paste(afmPATH,"/Arial.afm",sep=""),
                           paste(afmPATH,"/Arial Bold.afm",sep=""),
                           paste(afmPATH,"/Arial Italic.afm",sep=""),
                           paste(afmPATH,"/Arial Bold Italic.afm",sep="")))
      if(!"Arial" %in% names(pdfFonts())){pdfFonts(Arial=Arial)}
      if(!"Arial" %in% names(postscriptFonts())){postscriptFonts(Arial=Arial)}
      return()
    } else {disp(header='useArial=TRUE',message='The directory did not contain the *.afm version of the Arial font family')}
  } else {disp(header='useArial=TRUE',message='Please provide the path to the *.afm version of the Arial font family')}
}


plot.loglog <- function(fd.OUT){
  require(ggplot2)
  require(scales)
  g <- ggplot2::ggplot(fd.OUT$PLAW, aes(x=size,y=bulk), na.rm=T) +
    scale_x_log10(breaks = log_breaks(n=abs(diff(range(round(log10(fd.OUT$PLAW$size)))+c(-1,1))),base=10),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = range(round(log10(fd.OUT$PLAW$size)))+c(-1,1)) +
    scale_y_log10(breaks = log_breaks(n=abs(diff(range(round(log10(fd.OUT$PLAW$bulk)))+c(-1,1))),base=10),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = range(round(log10(fd.OUT$PLAW$bulk)))+c(-1,1)) +
    geom_point() +
    geom_abline(intercept = fd.OUT[[2]]$fitlm1$coefficients[[1]], slope = fd.OUT[[2]]$fitlm1$coefficients[[2]], colour = "red", size = 2) +
    ggtitle(paste("Regression over ",length(fd.OUT[[2]]$fitlm1$fitted.values)," frequencies/bins",sep=""))+
    xlab("Frequency (log10)")+ylab("Power (log10)") +
    annotation_logticks() +
    annotate("text",x=10^-2,y=10^5,label=paste("Slope = ",round(fd.OUT[[2]]$alpha,digits=2),sep="")) +
    gg.theme("clean")
  return(g)
}

# Variability Analyses --------------------------------------------------------------------------------------------------------------------------

#
# #' PSDslope
# #'
# #' @param y    A time series object, or a vector that can be converted to a time series object.
# #' @param fs    Sample frequency (defults to 1).
# #' @param nfft    Number of frequencies to estimate (defaults to next power of 2)
# #' @param fitRange    Vector of length 2 with range of frequencies to perform log-log fit.
# #' @param plot    Plot the log-log spectrum and slope.
# #'
# #' @return
# #' @export
# #'
# #' @examples
# #'
# PSDslope <- function(y  = ts(rnorm(n = 1024), frequency = 1),
#                      fs = frequency(y),
#                      nfft = 2^(nextpow2(length(y)/2)),
#                      fitRange = c(1,round(.1*nfft)),
#                      plot = FALSE){
#   require(oce)
#   require(signal)
#   if(!is.ts(y)){ts(y, frequency = fs)}
#
#   win <- signal::hamming(n=nfft)
#
#   perioGram <- oce::pwelch(x = y, window = win, fs = frequency(y), nfft = nfft, plot = FALSE)
#   spec <- data.frame(Frequency = perioGram$freq, Power = perioGram$spec)
#   spec[1,1:2] <- NA
#   fit <- lm(log10(spec$Power[fitRange[1]:fitRange[2]])~log10(spec$Power[fitRange[1]:fitRange[2]]))
#   return(list(spec = spec,
#               slope = fit)
#   )
# }

#' Scale.R
#'
#' @description Rescale a vector to a user defined range defined by user.
#'
#' @param x     Input vector or data frame.
#' @param mn     Minimum value of original, defaults to \code{min(x, na.rm = TRUE)}.
#' @param mx     Maximum value of original, defaults to \code{max(x, na.rm = TRUE)}.
#' @param hi     Minimum value to rescale to, defaults to \code{0}.
#' @param lo     Maximum value to rescale to, defaults to \code{1}.
#'
#'
#' @details Three uses:
#' \enumerate{
#' \item scale.R(x)             - Scale x to data range: min(x.out)==0;      max(x.out)==1
#' \item scale.R(x,mn,mx)       - Scale x to arg. range: min(x.out)==mn==0;  max(x.out)==mx==1
#' \item scale.R(x,mn,mx,lo,hi) - Scale x to arg. range: min(x.out)==mn==lo; max(x.out)==mx==hi
#' }
#'
#' @return
#' @export
#'
#' @examples
#' # Works on numeric objects
#' somenumbers <- cbind(c(-5,100,sqrt(2)),c(exp(1),0,-pi))
#'
#' scale.R(somenumbers)
#' scale.R(somenumbers,mn=-100)
#
#' # Values < mn will return < lo (default=0)
#' # Values > mx will return > hi (default=1)
#' scale.R(somenumbers,mn=-1,mx=99)
#'
#' scale.R(somenumbers,lo=-1,hi=1)
#' scale.R(somenumbers,mn=-10,mx=101,lo=-1,hi=4)
scale.R <- function(x,mn=min(x,na.rm=T),mx=max(x,na.rm=T),lo=0,hi=1){
  x <- as.data.frame(x)
  u <- x
  for(i in 1:dim(x)[2]){
    mn=min(x[,i],na.rm=T)
    mx=max(x[,i],na.rm=T)
    if(mn>mx){warning("Minimum (mn) >= maximum (mx).")}
    if(lo>hi){warning("Lowest scale value (lo) >= highest scale value (hi).")}
    ifelse(mn==mx,{u[,i]<-rep(mx,length(x[,i]))},{
      u[,i]<-(((x[i]-mn)*(hi-lo))/(mx-mn))+lo
      id<-complete.cases(u[,i])
      u[!id,i]<-0
    })
  }
  return(u)
}

# Rmd2htmlWP <- function(infile, outfile, sup = T) {
#   require(markdown)
#   require(knitr)
#   mdOpt <- markdownHTMLOptions(default = T)
#   mdOpt <- mdOpt[mdOpt != "mathjax"]
#   mdExt <- markdownExtensions()
#   mdExt <- mdExt[mdExt != "latex_math"]
#   if (sup == T) {
#     mdExt <- mdExt[mdExt != "superscript"]
#   }
#   knit2html(input = infile, output = outfile, options = c(mdOpt), extensions = c(mdExt))
# }

# MULTIPLOT FUNCTION ------------------------------------------------------------------------------------------------------------------
#
# [copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ ]
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multi.PLOT <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' TRY ... CATCH
#'
#' @details
#'  In longer simulations, aka computer experiments,
#'  you may want to
#'  1) catch all errors and warnings (and continue)
#'  2) store the error or warning messages
#'
#'  Here's a solution  (see \R-help mailing list, Dec 9, 2010):
#'
#' Catch *and* save both errors and warnings, and in the case of
#' a warning, also keep the computed result.
#'
#' @title tryCatch both warnings (with value) and errors
#' @param expr an \R expression to evaluate
#' @return a list with 'value' and 'warning', where value' may be an error caught.
#' @author Martin Maechler; Copyright (C) 2010-2012  The R Core Team
try.CATCH <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}
