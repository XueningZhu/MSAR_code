library(MASS)
source("LSE_estimation.R")
source("LSE_grad_hessian.R")
source("LSE_inference.R")
source("simulator.R")

### Parameter setup
D1vec = c(0.3,-0.2,0,0.1); D1 = matrix(D1vec, nrow = 2)
p = nrow(D1)
B1 = cbind(c(-0.5, 1), c(1.3, 0.3))
q = nrow(B1)
ZSigma = 0.5^abs(outer(1:2, 1:2,"-")) 
Sige1 = Diagonal(n = p, x = c(0.4,0.6))
Sige1[1,2] = 0.1
Sige1[2,1] = 0.1
theta = c(as.vector(D1), as.vector(B1))
net_type = 3 # use powerlaw network
is_nor = F # use t-distribution

### Initial value setup
D0 = Matrix(0, p, p)
B0 = Matrix(0, q, p)
Sige0 = Diagonal(n =2, x = c(0.4,0.5))
Sige0[1,2] = 0.15
Sige0[2,1] = 0.15

# MSAR demo
set.seed(1234)
N = 500 # set the network size

# generate the network
if (net_type=="1")
  W = getDyadW(N, N1 = 10, delta = 1.2, normalize = T)
if (net_type=="2")
  W = getBlockW(N, Nblock = 20, normalize = T)
if (net_type=="3")
  W = getPowerLawW(N = N, alpha = 2.5) 

# generate dataset
X = mvrnorm(n = N, mu = rep(0,nrow(ZSigma)), Sigma = ZSigma) # generate X matrix
Ymat = getYmat(D1, W, N, X, B1, Sige1, is_nor) # generate Ymat by MSAR model

# conduct LSE estimation
lse_res = MSAR.Lse(Ymat, X, W, D = D0, B = B0, Sige = Sige0) # LSE estimation
Cov = lse_res$cov # asymptotic covariance
lse_res$D
lse_res$B
theta = as.vector(lse_res$theta) # estimated parameter (vec(D), vec(B))
times = lse_res$time # time

# statistical inference (95% confidence interval)
CI_low = theta-1.96*sqrt(diag(Cov))
CI_up = theta + 1.96*sqrt(diag(Cov))
