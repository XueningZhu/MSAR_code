
#### This file includes functions which calculate the newton-raphson algorithm of LSE


### calculate the trace of a matrix
tr<-function(M)
{
  return(sum(diag(M)))
}



MSAR.Lse<-function(Ymat, X, W, D, B, Sige, verbose = F, infer = T, old = F, tol = 10^{-3})  ### LSE algorithm
{
  tim = system.time({
    N = nrow(W); p = ncol(Ymat); vecY = as.vector(Ymat); ww = crossprod(W); q = nrow(B)
    Omee = solve(Sige)
    Del1 = 1
    n.iter = 1
    flag = FALSE
    ind_d = 1:p^2; ind_b = (p^2+1):(p^2+p*q); ind_e0 = which(lower.tri(Sige,diag = T))
    obj0 = -1
    while (mean(abs(Del1))>tol & n.iter<=100)
    {
      if (verbose)
         cat(mean(abs(Del1)), " ")
      if (any(abs(diag(D))>1))
      {
        D = Diagonal(n = p, x = runif(p, -0.5, 0.5))
        B = Matrix(runif(p*q), nrow = q, ncol = p)
      }
      vecTheta = c(as.vector(D), as.vector(B))
      
      
      system.time({grad_hessian = lse.grad.hessian(D, Ymat, vecY, W, ww, N, p, X, B, Sige)})
      
      hessian = grad_hessian$hessian
      Del = solve(hessian)%*%grad_hessian$grad
      if (any(abs(diag(D))>1))
        vecTheta = vecTheta - Del*0.1
      else
        vecTheta = vecTheta - Del
      D = Matrix(matrix(vecTheta[ind_d], nrow = p))
      B = Matrix(matrix(vecTheta[ind_b], nrow = q))
      E = Ymat - W%*%Ymat%*%D - X%*%B
      Sige = t(E)%*%E/N
      Del1 = abs(grad_hessian$obj-obj0)
      obj0 = grad_hessian$obj
      n.iter = n.iter+1
    }
  })
  if (infer)
  {
    grad_hessian_cov = lse.grad.hessian(D, Ymat, vecY, W, ww, N, p, X, B, Sige, infer = infer, old = old)
  }else{
    grad_hessian_cov = NULL
  }
    
  return(list(D = D, B = B, Sige = Sige, theta = vecTheta, cov = grad_hessian_cov$cov,
              iter = n.iter, time = tim[3]))
}
