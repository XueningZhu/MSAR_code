MSAR.Lse.Sig1<-function(vecY, W, ww, N, p, q, D, Sige, Omee, 
                        S, tS, S1, m, OmeeS, tSe, SYX, IX, IXbeta)
{
  ### define frequently used matrices 
  In = Diagonal(N, x = 1)
  Inp = Diagonal(N*p, x = 1)
  
  
  ### by eigenvalue decomposition, we have \Sigma_e = Q Lambda t(Q)
  ### here Q = ee$vectors
  ee = eigen(Sige)
  Sige_half = Matrix(t(sqrt(ee$values)*t(ee$vectors))) # QLambda^{1/2}
  ISige_half = kronecker(Sige_half, In)
  tISige_half = kronecker(t(Sige_half), In)
  IOmee_half = kronecker(Omee%*%Sige_half, In)
  
  ### (I-t(D)\otimes W)^{-1}(X\beta)
  SIXbeta = S1%*%IXbeta
  
  ### matrices used
  mtSe = m*tSe
  m2tSe = m^2*tSe
  Sem2tSe = OmeeS%*%(m2tSe)
  
  ### first order: E(Q_{j1j2}^d Q_{k1k2}^d)
  A1 = list()
  A2 = list()
  A3 = list()
  A4 = list()
  G = matrix(0, nrow = N*p, ncol = p^2)
  dM = matrix(0, nrow = N*p, ncol = p^2)
  
  # mS1 = as.matrix(S1)
  tmp2 = tISige_half %*% tSe
  tmp3 = tISige_half %*% OmeeS
  #tmp4 = as.matrix(tmp3 %*% m2tSe)
  tmp4 = tmp3 %*% m2tSe
  
  right = S1%*%ISige_half
  
  # mS1 = as.matrix(S1)
  ### first order: E(Q_{j1j2}^d Q_beta)
  Sig1dx = matrix(0, p^2, p*q)
  SSX = -as.matrix(Sem2tSe) %*% S %*% as.matrix(m2tSe) %*% IX
  
  mOmeeS = as.matrix(OmeeS)
  mSem2tS = as.matrix(Sem2tSe)
  gc()
  #cat("now start computational\n")
  for (j1 in 1:p)
  {
    for (j2 in 1:p)
    {
      #cat(j1, j2, "\n")
      jj = (j2-1)*p+j1 # jj is the row
      Ij2j1 = Matrix(0, nrow = p, ncol = p)
      Ij2j1[j2, j1] = 1
      Ij1j2 = t(Ij2j1)
      
      ### matrices gradients
      Se_g = -kronecker(Omee%*%Ij2j1, W)
      S_g = -kronecker(Ij2j1, W)
      V_g = kronecker(Ij1j2%*%Omee%*%t(D)+D%*%Omee%*%Ij2j1, ww)
      m_g = -m^2*diag(V_g)
      
      ### calculate tr(Mj1j2 Mk1k2)
      A1[[jj]] = m*m_g*tS + m^2*t(S_g)
      A2[[jj]] = as.matrix(OmeeS%*%A1[[jj]])
      A3[[jj]] = m2tSe%*%S_g
      A4[[jj]] = as.matrix(mSem2tS %*% S_g)
      A1[[jj]] = as.matrix(A1[[jj]])
      A3[[jj]] = as.matrix(A3[[jj]])
      
      G1j = tmp3 %*% ((m*m_g)*tSe) %*% SYX
      G2j = tmp3 %*% ((m^2)*t(Se_g)) %*% SYX
      G3j = tISige_half %*% A4[[jj]] %*% vecY
      G[,jj] = G1j[,1] + G2j[,1] + G3j[,1]
      Sig1dx[jj,] = as.numeric(t(SIXbeta)%*%t(S_g)%*%SSX)
      #Sig1dx[jj,] = as.numeric(t(SSX)%*%S_g%*%vecY)
      
      ### calculate diag(Mj1j2)
      left = tmp4 %*% S_g
      dM[,jj] = rowSums(left*t(right)) + rowSums((tISige_half%*%A2[[jj]])* t(IOmee_half))
      gc()
    }
  }
  wdel = as.vector(kronecker((solve(Sige_half)), In)%*%SYX)
  del4 = mean(wdel^4) -  3*mean(wdel^2)^2
  GG = crossprod(G)
  dM2 = crossprod(dM)
  #cat("now no computational\n")
  ### calculate 
  Sig1d1 = matrix(0, p^2, p^2)
  # Sig1d2 = matrix(0, p^2, p^2)
  # Sig1d3 = matrix(0, p^2, p^2)
  # Sig1d4 = matrix(0, p^2, p^2)
  for (j1 in 1:p)
  {
    for (j2 in 1:p)
    {
      for (k1 in 1:p)
      {
        for (k2 in 1:p)
        {
          #cat(j1,j2,k1,k2,'\n')
          jj = (j2-1)*p+j1 # jj is the row
          kk = (k2-1)*p+k1 # kk is the column
          if (kk>=jj)
          {
            #cat(kk, jj, "\n")
            ### E(Q_{j1j2}^dQ_{k1k2}^d)
            ### the same as Sig1d[jj,kk] = as.numeric(tr(M_g[[jj]]%*%t(M_g[[kk]])) + tr(M_g[[jj]]%*%M_g[[kk]])+ t(U_g[[jj]])%*%U_g[[kk]])
            
            ### calculate tr(Mj1j2 Mk1k2)
            Sig1d1[jj,kk] = sum(A2[[jj]]*t(A2[[kk]]))+sum(A4[[jj]]*t(A1[[kk]]))+
              sum(A4[[kk]]*t(A1[[jj]])) + sum(A3[[jj]]*t(A3[[kk]]))
            
            # Sig1d2[jj,kk] = sum(M_g[[jj]] * (t(M_g[[kk]])))
            # 
            # Sig1d3[jj,kk] = sum(M_g[[jj]] * ( M_g[[kk]])) +
            #   sum(U_g[[jj]]*U_g[[kk]]) 
            
            # Sig1d4[jj,kk] = del4*sum(diag(M_g[[jj]])*diag(M_g[[kk]]))
            
            # Sig1d[jj,kk] = sum(M_g[[jj]] * (t(M_g[[kk]]) + M_g[[kk]])) +
            #   sum(U_g[[jj]]*U_g[[kk]]) +
            #   del4*sum(diag(M_g[[jj]])*diag(M_g[[kk]]))
            
          }
          gc()
        }
      }
    }
  }
  Sig1d1[lower.tri(Sig1d1)] = t(Sig1d1)[lower.tri(Sig1d1)]
  Sig1d = Sig1d1 + GG + dM2*del4
  
  ### first order: E(Q_beta Q_beta^\top)
  Sig1x = t(IX)%*%Sem2tSe%*%S%*%(m2tSe)%*%IX
  
  
  ### the exact form of Sig1 = E{Q(theta)Q(theta)^\top}
  Sig1 = matrix(0, nrow = p^2+p*q, ncol = p^2+p*q)
  Sig1[1:p^2, 1:p^2] = Sig1d
  Sig1[1:p^2, (p^2+1):(p^2+p*q)] = Sig1dx
  Sig1[(p^2+1):(p^2+p*q), 1:p^2] = t(Sig1dx)
  Sig1[(p^2+1):(p^2+p*q),(p^2+1):(p^2+p*q)] = as.matrix(Sig1x)
  Sig1 = 4*Sig1
  gc()
  
  return(Sig1)
}