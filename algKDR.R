algKDR <- function(dat, p=2, max_loop=20){
  # Code by Debo Cheng
  # The function of CESD algorithm
  # dat includes, x: covariates, Tr: treatment variable, Y: outcome
  # p #target reduced dimensions
  # max_loop     #number of iterations in kdr method
  require(KDRcpp)
  ##install package devtools if necessary
  #devtools::install_github('aschmu/KDRcpp')
  
  
m <- nrow(dat)
ns <- ncol(dat)
# data.frame transfer to matrix
X <- data.matrix(dat[,-c(ns-1,ns)])
## data processing 
Tr <- as.matrix(dat$Tr,ncol=1,nrow=m);
Y<-as.matrix(dat$Y,ncol=1,nrow=m);

starttime<-proc.time()
# p <- 2 #target reduced dimension
  
Xs <- scale(X)
sx <- 5
sy <- 1.4
eps <- 0.0001     #regularization parameter for matrix inversion	
eta <-10.0        #range of the golden ratio search
anl	<- 4         #maximum value for anealing
eta <- 10
verbose <- TRUE   #print the optimization process info?
#Gaussian kernels are used. Deviation parameters are set by the median of
#mutual distances. In the anealling, sigma chages to 2*median to
#0.5*median
sgx <- 0.5*sx
sgy <- sy  #Y is discrete, tuning is not necessary.
# Ky <- RBFdot(Y, Y, .5)
B <- kdr_trace_cpp(X = Xs, Y = Y, K = p, 
                   max_loop = max_loop,
                   sigmax0 = sx*sqrt(p/m),
                   sigmay0 = sy, eps = eps,
                   eta = eta, anl = anl, 
                   verbose = verbose,
                   tol = 1e-9)
# cputime <- system.time()
rx <- Xs%*%B
runtime <- proc.time() - starttime
Runtime <- runtime[1]

retu<-list(rx,B,Runtime)
names(retu) <- c("rx","B","Runtime")

return(retu)
}
