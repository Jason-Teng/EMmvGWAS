library(doParallel) # parallel computing
library(foreach) # parallel computing 
library(RSpectra) # spectral decomposition

###################  function for GWAS  ##################
eigen_trans <- function(K, npc = nrow(K)){
  message("Starting Eigen Decomposition of GRM...")
  start_time <- Sys.time()
  ######## eigen decomposition
  r <- suppressWarnings(eigs_sym(K, k = npc))
  # U <- r$vectors
  r$values[r$values<1e-5] =1e-5
  # D = r$values
  # H = sqrt(solve(diag(D)))%*%t(U)
  end_time <- Sys.time()
  message(paste("Eigen decomposition completed in", format(end_time-start_time)))
  return(r)
}


vc_estimation <- function(P, r, name="gwas", X = rep(1,nrow(P))){
  # “-o [prefix]” specifies output file prefix.
  ######### transformation
  message("Starting variance component estimation...")
  start_time <- Sys.time()
  
  if (is.vector(P)) {
    Y_sub <- matrix(P, ncol = 1)
  }
  if (is.null(colnames(P))) {
    colnames(P) <- paste0("trait", seq_len(ncol(P)))
  }
  
  n <- nrow(P)
  m <- ncol(P)
  U <- r$vectors
  D <- r$values
  H <- sqrt(solve(diag(D)))%*%t(U)
  
  Y <- H%*%P
  X <- H%*%X
  W <- D
  # Z0 <- H%*%t(Genotypes)
  
  # write.csv(data_wide, "data_wide_hybrid.csv",row.names = F)

  # EM algorithm
  aa <- diag(m)
  ee <- diag(m)
  maxiter <- 1000
  minerr <- 1e-8
  iter <- 0
  err <- 1e8  ## should update err
  www <- matrix(numeric(0), ncol=2+m)  # To store iteration results if needed
  
  while(iter < maxiter && err > minerr) {
    xx <- matrix(0, m, m)
    xy <- matrix(0, m, 1)
    yy <- matrix(0, m, m)
    bb <- matrix(0, m, m)
    
    for (j in 1:n) {
      Xj = diag(X[j],m)
      tt <- solve(aa * W[j] + ee+ diag(m) * 1e-5) #
      xx <- xx + W[j] * Xj %*% tt %*% Xj
      xy <- xy + W[j] * Xj %*% (tt %*% Y[j,])
    }
    
    mu <- solve(xx) %*% xy
    
    for (j in 1:n) {
      Xj = diag(X[j],m)
      tt <- solve(aa * W[j] + ee+ diag(m) * 1e-5 ) #+ diag(m) * 1e-5
      ex <- W[j] * aa %*% tt %*% (Y[j,] - Xj %*% mu)
      yy <- yy + W[j] * (Y[j,] - Xj %*% mu - ex)%*% Y[j,]
      va <- aa - W[j] * aa %*% tt %*% aa
      bb <- bb + (ex %*% t(ex) + va)
    }
    
    ee1 <- yy / (n - 1)
    aa1 <- bb / n
    
    ### info for each iteration
    www <- rbind(www, c(iter, err, t(mu)))#,ll
    err <- (sum((aa1 - aa)^2) + sum((ee1 - ee)^2)) / (2 * m^2)
    ee <- ee1
    aa <- aa1
    iter <- iter + 1
  }
  lambda = (aa)%*%solve(ee)
  write.csv(aa,paste(name, "G.csv",sep = ""),row.names = F)
  write.csv(ee,paste(name, "R.csv",sep = ""),row.names = F)
  write.csv(mu,paste(name, "beta.csv",sep = ""), row.names = F)
  end_time <- Sys.time()
  message(paste("Variance component estimation completed in", format(end_time-start_time)))
  return(lambda)
}


semi_gwas <- function(P, r, Z, lambda, name="gwas", X=rep(1,nrow(P)),num_cores=1){
  # Z: genexsample
  message("Starting GWAS...")
  start_time <- Sys.time()
  
  if (is.vector(P)) {
    P <- matrix(P, ncol = 1)
  }
  if (is.null(colnames(P))) {
    colnames(P) <- "trait"
  }
  p <- nrow(Z) 
  n <- nrow(P)
  m <- ncol(P)
  
  U <- r$vectors
  D = r$values
  H = sqrt(solve(diag(D)))%*%t(U)
  # Z<- H%*%t(Z)
  Y <- H%*%P
  X<- H%*%X
  W = D
  
  ## scanning the genome
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  result_scan <- foreach(k = 1:p, .combine = rbind, .packages = "Matrix") %dopar% {
    Z_k <- H%*%Z[k, ]
    # Z_k <- Z[,k] # this one for transform genotype first
    Axx <- matrix(0,m,m)
    Axz=matrix(0,m,m)
    Axy=matrix(0,m,1)
    Azy=matrix(0,m,1)
    Azz=matrix(0,m,m)
    Ayy=matrix(0,m,m)
    for (j in 1:n){
      Xj = diag(X[j,],m)
      Zj = diag(Z_k[j],m)
      tt=solve(lambda*W[j]+diag(m)*(1+1e-5))
      Axx=Axx+W[j]*Xj%*%tt%*%Xj
      Axy=Axy+W[j]*Xj%*%(tt%*%(Y[j,]))
      Axz=Axz+W[j]*Xj%*%tt%*%Zj
      Azy=Azy+W[j]*Zj%*%(tt%*%(Y[j,]))
      Azz=Azz+W[j]*Zj%*%tt%*%Zj
      # Axy=Axy+W[j]*Xj%*%(tt%*%(Y[j,]))
      Ayy=Ayy+W[j]*(tt%*%(Y[j,]))%*%t(Y[j,])
    }
    Azx <- Axz
    Ayx <- Axy
    Czz <- Azz - Azx %*% solve(Axx) %*% Axz
    Czy <- Azy - Azx %*% solve(Axx) %*% Axy
    Cyy <- Ayy - (solve(Axx) %*% Axy) %*% t(Ayx)
    Cyz <- Czy
    
    gg <- solve(Czz) %*% Czy
    Crr <- Cyy - solve(Czz) %*% Czy %*% t(Cyz)
    ee <- Crr / (n - 2)
    vv <- solve(Czz) %*% ee
    vvi <- solve(vv)
    wald <- t(gg) %*% vvi %*% gg
    pvalue <- exp(pchisq(wald,df=m, lower.tail = FALSE, log.p = TRUE))
    log10p <- -log10(pvalue)
    
    # Return the results as a vector
    c(k, t(gg), wald, pvalue, log10p)
  }
  # Stop the parallel cluster
  stopCluster(cl)
  
  colnames(result_scan)<- c("snp", colnames(Y), "wald", "p","log10p")
  result_scan = as.data.frame(result_scan)
  end_time <- Sys.time()
  message(paste("GWAS completed in", format(end_time-start_time)))
  write.csv(result_scan,paste(name, "_association.csv",sep = ""),row.names = FALSE)
  return(result_scan)
}

gwas_all <- function(K, P, Z, name="gwas", X=rep(1,nrow(P)),num_cores=1){
  r <- eigen_trans(K)
  lambda <- vc_estimation(P, r, name, X) 
  res <- semi_gwas(P, r, Z, lambda,name, X,num_cores)
  return(res)
}