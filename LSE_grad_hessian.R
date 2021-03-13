
ifelse<-function(x, val1, val0)
{
  if (x)
  {
    return(val1)
  }else
  {
    return(val0)
  }
}

### calculate the gradient and hessian matrix of Q(theta)
lse.grad.hessian<-function(D, Ymat, vecY, W, ww, N, p, X, B, Sige, infer = F, old = F)
{
  p = ncol(Ymat)
  q = ncol(X)
  ############################################################
  #### matrices
  
  tW = t(W) # ww = tW%*%W
  
  ### identity matrix
  I = Diagonal(N*p, x = 1)
  Ip = Diagonal(p, x = 1)
  In = Diagonal(N, x = 1)
  Inp = Diagonal(N*p, x = 1)
  
  
  IX = kronecker(Ip, X)
  IXbeta = as.vector(X%*%B)
  
  Omee = solve(Sige)
  IOmee = kronecker(Omee, In)
  
  S = I - kronecker(t(D), W); tS = t(S)
  SYX = S%*%vecY - IXbeta # This \mE
  OmeeS = IOmee%*%S
  tSe = t(OmeeS)
  OmeSYX = tSe%*%SYX
  
  Ome = tSe%*%S
  m = 1/diag(Ome)
  
  MY = m*OmeSYX # this is F
  
  dww = diag(ww)
  vMY = as.vector(MY)
  mtSe = m*tSe
  
  # left1 = matrix(-m^2*vMY, ncol = p)
  # left2 = matrix(m*vMY, ncol = p)
  # left3 = matrix(OmeeS%*%(m*vMY), ncol = p)
  
  ### calculate gradient
  
  right1 = matrix(tSe%*%SYX, ncol = p)
  right2 = matrix(SYX, ncol = p)
  right3 = Ymat
  
  ###
  left1b = matrix(m^2*vMY, ncol = p)
  left2b = matrix(-m*vMY, ncol = p)
  
  
  #grad_D = rep(0, p^2)
  F_d = matrix(0, nrow = N*p, ncol = p^2)
  Q_db2 = Matrix(0, nrow = p^2, ncol = p*q)
  
  for (j1 in 1:p)
  {
    for (j2 in 1:p)
    {
      #cat(j1, j2, "\n")
      Ij2j1 = Matrix(0, nrow = p, ncol = p)
      Ij2j1[j2, j1] = 1
      Ij1j2 = t(Ij2j1)
      
      Vj1j21 = Ij1j2%*%Omee%*%t(D)+D%*%Omee%*%Ij2j1
      # right1j = (dww*right1)%*%diag(diag(Vj1j21))
    
      F_d1j = ifelse(p==1, -m^2*as.vector((dww*right1)%*%(Vj1j21)),
                     -m^2*as.vector((dww*right1)%*%diag(diag(Vj1j21))))
      
      # right2j = -t(W)%*%right2%*%t(Ij1j2%*%Omee)
      F_d2j = -m*as.vector(t(W)%*%right2%*%t(Ij1j2%*%Omee))
      
      # right3j = -W%*%right3%*%t(Ij2j1)
      F_d3j = -mtSe%*%as.vector(W%*%right3%*%t(Ij2j1))
      
      
      ii = (j2-1)*p+j1 # j2 is the column
      #grad_D[ii] = sum(left1*right1j) + sum(left2*right2j)+ sum(left3*right3j)
      F_d[,ii] = F_d1j + F_d2j + F_d3j[,1]
      
      ### for the calculation of F_db
      ### 1st term: 
      ### 1.1: diag(Vj1j21)*Omee \otimes (dww*X)
      ### 1.2: -diag(Vj1j21)*(Omee%*%t(D))\otimes (dww*W%*%X)
      Q_db1j = t(dww*X)%*%left1b%*%(diag(Vj1j21)*Omee)-
        t(dww*t(W)%*%X)%*%left1b%*%(diag(Vj1j21)*(D%*%Omee))
      ### 2nd term: -t(Omee%*%Ij2j1) \otimes t(W)%%X
      Q_db2j = -t(tW%*%X)%*%left2b%*%t(Omee%*%Ij2j1)
      Q_db2[ii,] = as.vector(Q_db1j + Q_db2j)*2
    }
  }
  
  F_b = -m*tSe%*%IX 
  
  ### calculate Q_{j1j2}^d and Q_{\beta} by (2.1)
  #grad_D = 2*grad_D
  grad_D = 2*colSums(vMY*F_d)
  grad_beta = 2*colSums(vMY*F_b)
  grad = c(grad_D, grad_beta)
  
  
  #############################################################
  #### 2nd order: Hessian
  
  ### modify the hessian matrix computation
  
  ### first term
  left11 = matrix(2*m^3*vMY, ncol = p)
  left12 = matrix(-m^2*vMY, ncol = p)
  right1 = matrix(OmeSYX, ncol = p)
  
  ### second term
  #### for Omega_k1k2
  left2 = matrix(-m^2*vMY, ncol = p)
  right21 = Ymat
  right22 = X%*%B
  
  ### 3rd term: just as the 2nd one
  ### 4th term
  left4 = matrix(m*vMY, ncol = p)
  right4 = Ymat
  
  Q_dd2 = matrix(0, p^2, p^2)
  
  system.time({
    for (j1 in 1:p)
    {
      for (j2 in 1:p)
      {
        for (k1 in 1:p)
        {
          for (k2 in 1:p)
          {
            jj = (j2-1)*p+j1 # jj is the row
            kk = (k2-1)*p+k1 # kk is the column
            if (jj <= kk)
            {
              ii = (kk-1)*p^2 + jj
              
              
              #j1 = 1;j2 = 1; k1= 1;k2 = 1
              ### get Ij1j2 Ij2j1 and Ik1k2 Ik2k1
              Ij2j1 = Matrix(0, nrow = p, ncol = p)
              Ij2j1[j2, j1] = 1
              Ij1j2 = t(Ij2j1)
              Ik2k1 = Matrix(0, nrow = p, ncol = p)
              Ik2k1[k2, k1] = 1
              Ik1k2 = t(Ik2k1)
              IIk2k1 = kronecker(Ik2k1, In)
              
              Vj1j21 = Ij1j2%*%Omee%*%t(D)+D%*%Omee%*%Ij2j1
              Vk1k21 = Ik1k2%*%Omee%*%t(D)+D%*%Omee%*%Ik2k1
              Vjk = Ij1j2%*%Omee%*%Ik2k1+Ik1k2%*%Omee%*%Ij2j1
              
              ### 1st term
              right11jk = ifelse(p==1, (dww^2*right1)%*%(Vj1j21)*diag(Vk1k21),
                                 (dww^2*right1)%*%diag(diag(Vj1j21)*diag(Vk1k21)))
              right12jk = ifelse(p==1, (dww*right1)%*%(Vjk), (dww*right1)%*%diag(diag(Vjk)))
              
              ### 2nd term
              right211jk = -(dww*tW)%*%right21%*%t(diag(Vj1j21)*(Ik1k2%*%Omee))
              right212jk = -(dww*W)%*%right21%*%t(diag(Vj1j21)*(Omee%*%Ik2k1))
              right213jk = (dww*ww)%*%right21%*%t(diag(Vj1j21)*Vk1k21)
              right22jk = (dww*tW)%*%right22%*%t(diag(Vj1j21)*t(Omee%*%Ik2k1))
              
              ### 3rd term
              right211kj = -(dww*tW)%*%right21%*%t(diag(Vk1k21)*(Ij1j2%*%Omee))
              right212kj = -(dww*W)%*%right21%*%t(diag(Vk1k21)*(Omee%*%Ij2j1))
              right213kj = (dww*ww)%*%right21%*%t(diag(Vk1k21)*Vj1j21)
              right22kj = (dww*tW)%*%right22%*%t(diag(Vk1k21)*t(Omee%*%Ij2j1))
              
              ### 4th term
              right4jk = ww%*%right4%*%t(Vjk)
              
              Q_dd2[jj,kk] = sum(left11*right11jk)+sum(left12*right12jk)+
                sum(left2*right211jk) + sum(left2*right212jk) + sum(left2*right213jk)+ sum(left2*right22jk)+
                sum(left2*right211kj) + sum(left2*right212kj) + sum(left2*right213kj)+ sum(left2*right22kj)+
                sum(left4*right4jk)
              
            }
          }
        }
        
      }
    }
  })
  
  Q_dd2[lower.tri(Q_dd2)] = t(Q_dd2)[lower.tri(Q_dd2)]
  
  ### calculate (2.8)--(2.10)
  Q_dd = 2*t(F_d)%*%F_d + Q_dd2*2
  Q_db = 2*t(F_d)%*%F_b+Q_db2
  Q_bb = 2*t(F_b)%*%F_b
  
  
  ### Hessian Matrix
  Q_h = Matrix(0, nrow = p^2+p*q+p^2, ncol = p^2+p*q+p^2)
  ind_d = 1:p^2; ind_b = (p^2+1):(p^2+p*q)
  
  Q_h[ind_d, ind_d] = Q_dd
  Q_h[ind_b, ind_b] = Q_bb
  
  Q_h[ind_d, ind_b] = Q_db
  Q_h[ind_b, ind_d] = t(Q_db)
  
  ind0 = c(ind_d, ind_b)
  Q_h = Q_h[ind0, ind0]
  
  if (infer) # if inference is required
  {
    tS = t(S)
    ### fast calculation of (I-t(D)\otimes W)^{-1}\approx I + dw+dw^2+dw^3...
    dw = kronecker(t(D), W)
    W2 = (W%*%W)
    W3 = W2%*%(W)
    dw2 = kronecker(t(D%*%D), W2)
    dw3 = kronecker(t(D%*%D%*%D),W3)
    S1 = Inp + dw + dw2 + dw3
    
    if (old)
    {
      Sig1 = MSAR.Lse.Sig1.old(vecY, W, ww, N, p, q, D, Sige, S, m, OmeeS, tSe, Ome, Omee, IX, IXbeta, SYX)
    }else{
      Sig1 = MSAR.Lse.Sig1(vecY, W, ww, N, p, q, D, Sige, Omee,
                           S, tS, S1, m, OmeeS, tSe, SYX, IX, IXbeta)
    }
    
    inv_Qh = solve(Q_h)
    covl = inv_Qh%*%Sig1%*%inv_Qh
    return(list(grad = grad, hessian = Q_h, cov = covl))
  }
  
  gc()
  return(list(grad = grad, hessian = Q_h, obj = mean(vMY^2)))
}