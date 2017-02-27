#Hierarchical GP
library(TruncatedNormal)
library(dtmvnorm)
library(dplyr)



GP_cov <- function(d,lambda,ell,nugget = 0.){

  out_mat <- lambda^(-1)*exp(-as.matrix(dist(d))^2/ell^2) + nugget*diag(length(d))
  
  return(out_mat)
}

GP_cross_cov <- function(d,d_star,lambda,ell){
  inds <- 1:length(d)
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(c(d,d)))[inds,-inds]^2/ell^2)
  
  return(out_mat)
}

#Data stuff
read_hist_data <- function(hist_loc, q_file){
  q_vals <- read.table(q_file)[,1]
  
  Y_mat <- lapply(as.list(q_vals),function(x){
    read.table(file = paste0(hist_loc,"q_",as.character(x),".dat"))[,2]
  }) %>%
    do.call(rbind,.)
  
  bin_left = read.table(file = paste0(hist_loc,"q_",as.character(q_vals[1]),".dat"))[,1]
  
  colnames(Y_mat) <- as.character(bin_left)
  
  return(Y_mat)
  
}

jet_path <- "/Users/Jake/Dropbox/Research/JETSCAPE/JETSCAPE-STAT/"
hist_folder <- "q_dat_100k/"

q_vals_file <- paste0(jet_path,"qhat_vals_100k.dat")


#Data
num_counts <- 100000
Y <- read_hist_data(hist_loc = paste0(jet_path,hist_folder),
                    q_file = q_vals_file) *num_counts
d <- read.table(q_vals_file)[,1]

I <- dim(Y)[1]
J <- dim(Y)[2]
alpha <- c(as.numeric(colnames(Y)),1)
jitter = FALSE
if(jitter){
  alpha[-c(1,J+1)] <- alpha[-c(1,J+1)] + rnorm(J-1,0,1E-3)
}
b <- (alpha - lag(alpha))[-1]

#Parameter things
r <- 0.5
N <- (J-1)-(J-1)%%2 #Closest even number less than the number of bins J
R <- diag(rep(r^(1:(N/2)),2))

make_coefs <- function(){
  A <- B <- matrix(0,N/2,J)
  for(n in 1:(N/2)){
    for(j in 2:(J+1)){
      A[n,j-1] <- (sin(2*pi*n*alpha[j]) - sin(2*pi*n*alpha[j-1]))/(2*pi*n)
      B[n,j-1] <- (cos(2*pi*n*alpha[j-1]) - cos(2*pi*n*alpha[j]))/(2*pi*n)
    }
  }
  
  return(list(A=A,B=B))
}

A <- make_coefs()$A
B <- make_coefs()$B
C <- cbind(t(A),t(B))

#Priors
lam_a <- rep(1,N)
lam_b <- rep(1,N)

ell_a <- rep(1,N)
ell_b <- rep(1,N)



hier_gp_mh <- function(iters = 1E4, burnin_prop = 0.1,
                       delta_lam = rep(0.3,N),
                       delta_ell = rep(0.3,N),
                       X_kap = replicate(N,rep(1E-1,I)), #column n is diagonal of proposal
                                                        #for GP n
                       verbose = FALSE
                       ){
  burnin <- iters*burnin_prop
  
  X_array <- array(0,dim = c(iters,I,N)) #Columns are GPs, rows are histograms
  lam_mat <- ell_mat <- matrix(0,iters,N)

  #Current values
  X_cur <- matrix(0,I,N)
  ell_cur <- rgamma(N,ell_a,ell_b)
  lam_cur <- rgamma(N,lam_a,lam_b)
  
  p_cur <- t(replicate(I,b))
  l_cur <- sum(Y*log(p_cur))
  
  ell_acc <- lam_acc <- x_acc <- numeric(N)
  
  time_start <- proc.time()
  for(t in 1:(burnin + iters)){
    for(n in 1:N){
      if(!t%%100){
        if(verbose){
          flush.console()
          cat("\r t = ", t, "Elapsed Time: ",as.numeric((proc.time() - time_start)[3])) 
          #print(paste0("t = ", t, ", n = ", n))
        }
      }
      
      #Update X_n
      auto_reject = FALSE
      x_cur_n <- X_cur[,n]
      
      #Need to find upper and lower bounds for each histogram
      l <- -Inf*(numeric(I) + 1)
      u <- Inf*(numeric(I) + 1)
      for(i in 1:I){
        
        
        #All the constraints
        constr <- (-b - sqrt(2)*C[,-n]%*%R[-n,-n]%*%X_cur[i,-n])/(sqrt(2)*R[n,n]*C[,n])
        
        #Lower constraints
        (a_constr <- constr[which(round(C[,n],10)>0)])
        
        #Upper constraints
        (b_constr <- constr[which(round(C[,n],10)<0)])
        
        #Bad constraints
        (bad_constr <- which(round(C[,n],10)==0))
        
        #If any bad constraints
        if(length(bad_constr)){
          for(j in bad_constr){
            
            #For each bad constraint, check if the rest of X[i,] take care of it
            #If they don't auto-reject the whole vector
            if(b[j] + sqrt(2)*C[j,-n]%*%R[-n,-n]%*%X_cur[i,-n]<=0){
              print(paste("Bad times: i=",i,"n =",n,"j =",j))
              auto_reject = TRUE
              l_star = -Inf
            }
          }
        }
        
        l[i] <- max(a_constr)
        u[i] <- min(b_constr)
        
      }# i loop
      
      #Draw truncated normal proposal
      #See help file for instruction
      x_star_n <- mvrandn(l = l - x_cur_n,
                          u = u - x_cur_n, 
                          Sig = diag(X_kap[,n]),
                          n = 1) + 
        x_cur_n
      
      #Calculate proposed likelihood based on proposed x_star_n
      X_star <- X_cur
      X_star[,n] <- x_star_n
      
      
      p_star <- t(sqrt(2)*C%*%R%*%t(X_star) + replicate(I,b))
      if(sum(p_star<0)){
        print(paste("j = ",j,"i = ",i,"t = ",t,"n = ",n))
        stop("Dammit this shouldn't happen")
      }
      
      l_star <- sum(Y*log(p_star))
      
      #adjusted <- dtmvnorm(x_cur_n,mean = x_star_n,sigma = diag(X_kap[,n]),lower = l,upper = u,log = TRUE) -
       # dtmvnorm(x_star_n,mean = x_cur_n,sigma = diag(X_kap[,n]),lower = l,upper = u,log = TRUE)
      #Ratio
      ratio <- l_star - l_cur +
        
        #Priors
        dmvnorm(x=x_star_n,sigma = GP_cov(d,lambda = lam_cur[n], ell = ell_cur[n], nugget = 1E-8),
                log = TRUE) - 
        dmvnorm(x=x_cur_n,sigma = GP_cov(d,lambda = lam_cur[n], ell = ell_cur[n], nugget = 1E-8),
                log = TRUE)  #+
        
        #Proposals
        #Arrggg Botev doesn't include pdfs
    #    adjusted
      
#       if(n %in% c(1,10)){
#         print(paste('xst',x_star_n))
#         print(paste('xcur',x_cur_n))
#         print(paste('l',l))
#         print(paste('u',u))
#       }
      
      
      if(runif(1)<exp(ratio)){
        x_acc[n] <- x_acc[n] + 1
        l_cur <- l_star
        X_cur[,n] <- x_star_n
      }
      
      #############
      #Update ell_n
      #############
      z_ell <- rnorm(1)
      ell_st <- ell_cur[n]*exp(delta_ell[n]*z_ell)
      adjust_ratio <- delta_ell[n]*z_ell
      
      #MH Ratio
        
      ratio <- 
        #"Likelihood" - X prior
        dmvnorm(X_cur[,n],sigma = GP_cov(d,lam_cur[n],ell_st,nugget = 1E-8),log = TRUE) -
        dmvnorm(X_cur[,n],sigma = GP_cov(d,lam_cur[n],ell_cur[n],nugget = 1E-8), log = TRUE) +
        
        #Priors
        dgamma(ell_st,ell_a[n],ell_b[n],log = TRUE) -
        dgamma(ell_cur[n],ell_a[n],ell_b[n],log = TRUE) +
        
        #Proposal Adjustment
        adjust_ratio
        
      if(ratio + rexp(1)>0){
        ell_acc[n] <- ell_acc[n] + 1
        ell_cur[n] <- ell_st
      }
      
      ############
      #Update lam_n
      ############
      z_lam <- rnorm(1)
      lam_st <- lam_cur[n]*exp(delta_lam[n]*z_lam)
      adjust_ratio <- delta_lam[n]*z_lam
      
      ratio <- 
        #"Likelihood" - X prior
        dmvnorm(X_cur[,n],sigma = GP_cov(d,lam_st,ell_cur[n],nugget = 1E-8),log = TRUE) -
        dmvnorm(X_cur[,n],sigma = GP_cov(d,lam_cur[n],ell_cur[n],nugget = 1E-8), log = TRUE) +
        
        #Priors
        dgamma(lam_st,lam_a[n],lam_b[n],log = TRUE)-
        dgamma(lam_cur[n],lam_a[n],lam_b[n],log = TRUE) +
        
        #Proposal Adjustment
        adjust_ratio
      
      
      if(ratio + rexp(1)>0){
        lam_acc[n] <- lam_acc[n] + 1
        lam_cur[n] <- lam_st
      }
      
      
    }# n loop
    
    if(t > burnin){
      X_array[t-burnin,,] <- X_cur
      lam_mat[t-burnin,] <- lam_cur
      ell_mat[t-burnin,] <- ell_cur
    }
  }# t loop
  
  return(list(X = X_array,lam = lam_mat,ell = ell_mat,
              x_acc = x_acc/(iters+burnin),
              lam_acc = lam_acc/(iters+burnin),
              ell_acc = ell_acc/(iters+burnin)))
}

Rprof('Hope_against_hope.out')
hope <- hier_gp_mh(iters = 4000,verbose = TRUE,
                   X_kap = replicate(N,rep(1E-2,I)))
Rprof()

summaryRprof('Hope_against_hope.out')





