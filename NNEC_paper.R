if(system.file(package = "FNN") == "") install.packages("FNN")
if(system.file(package = "knn.covertree") == "") install.packages("knn.covertree")
if(system.file(package = "Matrix") == "") install.packages("Matrix")
library(FNN)
library(knn.covertree)
library(Matrix)


### Nearest Neighbour Equilibrium Clustering (NNEC)
### An implementation of the algorithm described in
### the work "Nearest Neighbour Equilibrium Clustering",
### D.P. Hofmeyr, Arxiv...
## Arguments:
# X = matrix (n x d) with observations row-wise. 
# nn_opts = integer vector of options for number of nearest neighbours
# lambda_opts = (strictly positive) numeric vector of options for balancing threshold
# itmax = maximum number of iterations for finding a single equilibrium cluster
# cycleLength = number of previous iterations over which to check for cycling

NNEC <- function(X, nn_opts = 5*2:5, lambda_opts = 1+0:10/5, itmax = 100, cycleLength = 5, distance = 'Euclidean'){
  n <- nrow(X)
  d <- ncol(X)
  
  ### find neighbours up to maximum nn
  if(distance=='Euclidean') nns <- get.knn(X, max(nn_opts))
  else if(distance=='cosine') nns <- list(nn.index = find_knn(X-matrix(colMeans(X), n, d, byrow = TRUE), max(nn_opts), distance = 'cosine')$index)
  
  best <- -Inf ## tracks the quality of the best solution found so far
  
  ### iterate over different values for nn and lambda
  for(ni in 1:length(nn_opts)){
    nn <- nn_opts[ni]
    
    ### tL is a matrix with points' neighbours stored column-wise
    tL <- sparseMatrix(j = rep(1:n, nn), i = nns$nn.index[,1:nn], x = 1, dims = c(n, n))
    
    for(di in 1:length(lambda_opts)){
      lam <- lambda_opts[di]
      
      ### C stores equilibrium cluster membership strengths
      C <- matrix(0, n, 0)
      
      ### add new equilibrium clusters until every point is assigned to at least one
      while(min(Cl <- rowSums(C))==0){
        clus <- numeric(n) ## start with an empty cluster
        ix0 <- which(Cl==0)[which.max(rowSums(tL)[Cl==0])]
        clus[ix0] <- 1 ## add seed point with the greatest number of "reverse neighbours"
        clusOld <- sapply(1:cycleLength, function(r) clus+1) ## check for convergence/cycling
        it <- 1
        while(min(apply(clusOld, 2, function(z) mean((z-clus)^2)))>1e-7 & it < itmax){
          clusOld[,-1] <- clusOld[,-cycleLength]
          clusOld[,1] <- clus
          ### update cluster
          sc <- sum(clus)
          if(sc==0) break
          else if(sc==1) clus <- (tL[clus==1,] > sc/n*lam*nn)
          else clus <- (colSums(tL[clus==1,]) > sc/n*lam*nn)
          
          it <- it + 1
        }
        
        ### add equilibrium cluster and membership strengths to the collection
        cc <- colSums(tL*clus)/nn - mean(clus)*lam
        C <- cbind(C, cc*(cc>0))
        
        if(C[ix0,ncol(C)] == 0){
          ### if seed point is no longer in the cluster add a singleton cluster to stop re-selection as a seed
          ### Allocate membership strength to small eps > 0 so that if the seed
          ### belongs naturally to another equilibrium cluster it will rather be allocated
          ### there, but not set to zero otherwise it will be re-selected as seed
          C <- cbind(C, 0)
          C[ix0,ncol(C)] <- 1e-10
        }
        
      }
      
      ### normalise membership strengths
      C <- C/rowSums(C)
      
      ### if the quality of the solution is better than the best found so far, then store the solution
      prf <- mean(apply(C, 1, max), na.rm = T)
      if(prf > best){
        nnbest <- nn
        lambda <- lam
        best <- prf
        ddiff <- C
      }
    }
  }
  ### return final clustering assignment, normalised membership strength and best settings
  list(C = ddiff, cluster = as.numeric(as.factor(apply(ddiff, 1, which.max))), nn = nnbest, lambda = lambda)
}
