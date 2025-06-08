## Method II
# @h the maximum number of meaningful effects except the overall mean
# @X X is the model matrix without the all-one column
library(stringr)

GenerateTWOint <- function(M = list("main"=NULL,"intr"=NULL),string, h, M.store = list()) {
  if (is.null(string) || length(M$main)==h)
    return(M.store)
  
  store <- M.store
  l <- length(store)
  sl <- length(string)
  
  for (i in 1:sl) {
    l<- l+1
    if (is.null(M$main)){
      newM <- list("main"=c(string[i]),"intr"=NULL)
    }else {
      newM <- list("main"=c(M$main,string[i]),"intr"=NULL)
      if(length(M$main)==1) {
        newM$intr <- paste0(M$main,":",string[i])
      } else {newM$intr <- c(M$intr,paste0(M$main,":",string[i]),paste0(M$intr,":",string[i]))}
    }
    store[[l]] <- newM
    
    if (i==sl) 
      return(store)
    
    newstring <- string[(i+1):sl]
    store <- GenerateTWOint(newM,newstring,h,store)
    l <- length(store)
  }
}

hardPruning <- function(M.store=list(),h) {
  if(length(M.store)==0)
    return(NULL)
  
  Prun.store <- list()
  count.PM <- 0
  l <- length(M.store)
  for(i in 1:l){
    main <- M.store[[i]]$main
    intr <- M.store[[i]]$intr
    
    if(length(main)==h) {
      count.PM <- count.PM + 1
      Prun.store[[count.PM]]$main <- M.store[[i]]$main
      Prun.store[[count.PM]]$intr <- NULL
    } else if (length(main) + length(intr) <= h) {
      count.PM <- count.PM + 1
      Prun.store[[count.PM]] <- M.store[[i]]
    } else {
      k <- h - length(main)
      intr.cand <- combn(intr,k)
      no.intr.cand <- ncol(intr.cand)
      for(j in 1:no.intr.cand) {
        count.PM <- count.PM + 1
        Prun.store[[count.PM]]$main <- M.store[[i]]$main
        Prun.store[[count.PM]]$intr <- as.vector(intr.cand[,j])
      }
    } 
  }
  return(Prun.store)
}

getFullModel <- function(X) {
  # generate the full model
  d <- ncol(X)
  namc <- colnames(X)
  if(is.null(namc)) {
    colnames(X) <- paste0("X",1:d)
    namc <- colnames(X)
  }

  main <- namc[1:2]
  intr <- paste0(namc[1],":",namc[2])
  
  for(i in 3:d) {
    intr <- c(intr, paste0(main,":",namc[i]),paste0(intr,":",namc[i]))
    main <- c(main, namc[i])
  }
  
  return(list("main"=main,"intr"=intr))
}

getFormula <- function(M.list=list("main"=NULL,"intr"=NULL)) {
  total <- c(M.list$main,M.list$intr)
  if(length(total)==0)
    return( as.formula("y~1") )
  
  return(as.formula(paste0("y~", paste0(total,collapse = "+"))))
}

get.All.Formula <- function(M.store) {
  return(lapply(M.store, getFormula))
}

MallowCp <- function(partial_res, n, p,hat_v) {
  # `partial_res` is the the residual sum of squares on a training set of data
  #         under the partial model with `p` predictors
  # `n` is the sample size; `p` is the number of predictors
  # `hat_v` is the estimate of the variance under the full model
  (partial_res + 2*p*hat_v^2)/n
}

BIC <- function(partial_res, n, p) {
  n*log(partial_res/n) + (p+2)*log(n)
}

r_sq <- function(res,tot,df_res,df_tot) {
  1-(res/df_res)/(tot/df_tot)
}

GlobalSearching <- function(X,y,h) { 
  #-# `h` is the number of predictors in the model except the all-one column
  # Outputs
  p <- list("no.main"=NULL, "formula"=NULL,"Cp"=NULL,"BIC"=NULL,"r_sq"=NULL)
  
  # get column names of `X`
  np <- dim(X)
  n <- np[[1]]
  p <- np[[2]]
  namc <- colnames(X)
  if(is.null(namc)) {
    colnames(X) <- paste0("X",1:p)
    namc <- colnames(X)
  }
  # mean y
  my <- mean(y)
  tot <- sum((y-my)^2)
  
  # smaller or equal to `h`-size models
  raw.store <- GenerateTWOint(string = namc, h=h)
  hardPruning <- hardPruning(raw.store,h=h)
  ALLSET <- get.All.Formula(raw.store)
  numMdl <- length(ALLSET)
  for(i in 1:numMdl) {
    p$no.main <- c(p$no.main, length(hardPruning[[i]]$main))
  }
  p$formula <- list()
  
  # the estimate of the variance under the full model
  DX <- as.data.frame(X)
  DX$y <- y
  full.f <- getFormula( getFullModel(X))
  ff <- lm(full.f,data=DX)
  hat_v <- sum(ff$residuals^2)/ff$df.residual
  
  # calculate the Cp, BIC and r_sq
  for(mI in 1:numMdl) {
    p$formula[[mI]] <- ALLSET[[mI]] 
    f <- lm(ALLSET[[mI]], data =DX)
    p$Cp <- c(p$Cp, MallowCp( sum(f$residuals^2),n,p,hat_v ))
    p$BIC <- c(p$BIC, BIC(sum(f$residuals^2),n,p ))
    p$r_sq <- c(p$r_sq,r_sq(sum(f$residuals^2)),tot,f$df.residual,n-1)
  }
  
  # return results
  nCp <- which.max(p$Cp)
  nBIC <- which.max(p$BIC)
  nr_sq <- which.max(p$r_sq)
  p$maxCp <- p$formula[[nCp]]
  p$maxBIC <- p$formula[[nBIC]]
  p$maxr_sq <- p$formula[[nr_sq]]
  return(p)
}
