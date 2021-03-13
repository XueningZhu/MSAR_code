
library(MASS)
library(Matrix)
# library(poweRlaw)                                                                                     ### for power-law distribution

library(methods)
require(VGAM)

rpldis = function (n, xmin, alpha, discrete_max = 10000) 
{
    if (length(n) > 1L) 
        n = length(n)
    xmin = floor(xmin)
    u = runif(n)
    if (discrete_max > 0.5) {
        constant = zeta(alpha)
        if (xmin > 1) 
            constant = constant - sum((1:(xmin - 1))^(-alpha))
        cdf = c(0, 1 - (constant - cumsum((xmin:discrete_max)^(-alpha)))/constant)
        dups = duplicated(cdf, fromLast = TRUE)
        if (any(dups)) 
            cdf = cdf[1:which.min(!dups)]
        rngs = as.numeric(cut(u, cdf)) + xmin - 1
        is_na = is.na(rngs)
        if (any(is_na)) 
            rngs[is_na] = floor((xmin - 0.5) * (1 - u[is_na])^(-1/(alpha - 
                1)) + 0.5)
    }
    else {
        rngs = floor((xmin - 0.5) * (1 - u)^(-1/(alpha - 1)) + 
            0.5)
    }
    rngs
}

getYmat<-function(D, W, N, X, B, Sige, is_nor)
{
  p = nrow(D)
  In = Diagonal(n = N, x = 1)
  if (is_nor)
    eps = rnorm(N*p)
  else
    eps = rt(N*p, df = 5)
  ee = eigen(Sige)
  eps = kronecker(t(sqrt(ee$values)*t(ee$vectors)), In)%*%eps
  Xbeta = as.vector(X%*%B)
  eps1 = eps+Xbeta
  DkW = kronecker(t(D),(W))
  #InvDkW = solve(Diagonal(N*p)-DkW)
  #Yvec = InvDkW%*%eps
  DkW2 = DkW%*%DkW; DkW3 = DkW2%*%DkW; DkW4 = DkW3%*%DkW
  
  #D2 = D%*%D; D3 = D2%*%D; D4 = D3%*%D
  #W2 = W%*%W; W3 = W2%*%W; W4 = W3%*%W
  #DkW2 = kronecker(t(D2),(W2)); DkW3 = kronecker(t(D3),(W3)); DkW4 = kronecker(t(D4),(W4))
  Yvec = eps1 + DkW%*%eps1 + DkW2%*%eps1+ DkW3%*%eps1+ DkW4%*%eps1
  Ymat = matrix(Yvec, ncol = p) 
  return(Ymat)
}



getYmat.gamma<-function(D, W, N, X, B, Sige, is_nor)
{
  p = nrow(D)
  In = Diagonal(n = N, x = 1)
  eps = rgamma(N*p, shape = 1.5, scale = 2)
  eps = eps - mean(eps)
  ee = eigen(Sige)
  eps = kronecker(t(sqrt(ee$values)*t(ee$vectors)), In)%*%eps
  Xbeta = as.vector(X%*%B)
  eps1 = eps+Xbeta
  DkW = kronecker(t(D),(W))
  #InvDkW = solve(Diagonal(N*p)-DkW)
  #Yvec = InvDkW%*%eps
  DkW2 = DkW%*%DkW; DkW3 = DkW2%*%DkW; DkW4 = DkW3%*%DkW
  
  #D2 = D%*%D; D3 = D2%*%D; D4 = D3%*%D
  #W2 = W%*%W; W3 = W2%*%W; W4 = W3%*%W
  #DkW2 = kronecker(t(D2),(W2)); DkW3 = kronecker(t(D3),(W3)); DkW4 = kronecker(t(D4),(W4))
  Yvec = eps1 + DkW%*%eps1 + DkW2%*%eps1+ DkW3%*%eps1+ DkW4%*%eps1
  Ymat = matrix(Yvec, ncol = p) 
  return(Ymat)
}


getYmat.heter<-function(D, W, N, X, B, Sige, is_nor)
{
  p = nrow(D)
  In = Diagonal(n = N, x = 1)
  if (is_nor)
    eps = rnorm(N*p)
  else
    eps = rt(N*p, df = 5)
  d = (rowSums(W>0))
  eps = eps*rep(d, p) #rnorm(N*p, sd = 3)
  ee = eigen(Sige)
  eps = kronecker(t(sqrt(ee$values)*t(ee$vectors)), In)%*%eps
  Xbeta = as.vector(X%*%B)
  eps1 = eps + Xbeta
  DkW = kronecker(t(D),(W))
  #InvDkW = solve(Diagonal(N*p)-DkW)
  #Yvec = InvDkW%*%eps
  DkW2 = DkW%*%DkW; DkW3 = DkW2%*%DkW; DkW4 = DkW3%*%DkW
  
  #D2 = D%*%D; D3 = D2%*%D; D4 = D3%*%D
  #W2 = W%*%W; W3 = W2%*%W; W4 = W3%*%W
  #DkW2 = kronecker(t(D2),(W2)); DkW3 = kronecker(t(D3),(W3)); DkW4 = kronecker(t(D4),(W4))
  Yvec = eps1 + DkW%*%eps1 + DkW2%*%eps1+ DkW3%*%eps1+ DkW4%*%eps1
  Ymat = matrix(Yvec, ncol = p) 
  return(Ymat)
}



getDyadW<-function(N, N1 = 10, delta = 1.2, normalize = T)                                             ### simulate Dyad network W; N1: number of mutual pairs 2N1*N, delta: P((0,1)) = P((1,0)) = 0.5*N^{delta}
{
  A = matrix(0, nrow = N, ncol = N)                                                                    ### use A to store network structure
  
  ######################################### mutual follower ###########################################
  ind = which(upper.tri(A), arr.ind = T)                                                               ### record all the index of upper triangular matrix
  indM = ind[sample(1:nrow(ind), N*N1),]                                                               ### sample N*N1 as mutual pairs in the upper triangular matrix
  A[indM] = 1                                                                                          ### the corresponding links are set to be 1
  A[indM[,2:1]] = 1                                                                                    ### the following matrix is set to be symmetric, as a result, mutual pairs are 2N1*N
  
  ######################################### single relationship #######################################
  ind1 = which(A==0&upper.tri(A), arr.ind = T)                                                         ### record all the zero pairs in the upper triangular matrix
  indS = ind1[sample(1:nrow(ind1), N^delta),]                                                          ### choose N^delta index as single relations
  tmp = sample(1:nrow(indS), floor(N^delta/2))                                                         ### randomly choose 0.5*N^delta to be in the lower triangular matrix
  indS[tmp,] = indS[tmp, 2:1]                                                                          ### change the corresponding index to be inverse
  A[indS] = 1                                                                                          ### the single following relation is set to be 
  diag(A) = 0                                                                                          ### aii = 0
  if (!normalize)
    return(A)
  W = A/rowSums(A)                                                                                     ### W is row-normalized
  W = as(W, "dgCMatrix")
  return(W)
}

getPowerLawW<-function(N, alpha, normalize = T)                                                        ### get power-law network W
{
  Nfollowers = rpldis(N, 1, alpha)                                                                     ### generate N random numbers following power-law(1, alpha): k1-kN
  A = sapply(Nfollowers, function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  diag(A) = 0
  ind = which(rowSums(A)==0)                                                                           ### in case some row sums are zero
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1                                                                ### for those node, randomly select 3 followees
  }
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  W = as(W, "dgCMatrix")
  return(W)
}


getBlockW<-function(N, Nblock, normalize = T)                                                          ### get block network
{
  if (N%%Nblock==0){                                                                                   ### if N mod Nblock is integer
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)                        ### obtain the diagnal block list
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),                       ### generate following relations within the blocks
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)        ### if N mod Nblock is not integer
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),                  ### generate following relations within the blocks
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),                 ### generate following relations within the blocks
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)                                                                           ### combine the blocks in matrix
  offDiag = which(isDiag == 0, arr.ind = T)                                                            ### to calculate the index of the off digonal indexes
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)                                          ### people between blocks have 0.3 prob to follow
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  ################ transform bA to be a symmetric matrix ##############################################
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  
  ind = which(rowSums(bA)==0)                                                                          ### in case some row sums are zero
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                                                               ### for those node, randomly select 3 followees
  }
  
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                                                                   ### row normalize bA
  W = as(W, "dgCMatrix")
  return(W)
}