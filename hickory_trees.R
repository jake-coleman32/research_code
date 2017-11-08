library(gplots)

x <- c(31.4,112.4,56.2,46.6,126.2,64.7,112.4,120,11.1,24.8,113.9,18.2,44.3,114,125.4,
       33.5,113.4,12.2,46.7,126.7,83.2,115.5,10.2,12.5,32.8,117.2,31.6,51.2,114.6,
       38.7,11.6,21.5,72.9,10.3,19.2,118,10.2,18.5,38.3,106.,34.3,58,117,
       48.7,56.3,33.6,117.4,34,104.9,119.2,10.8,23.3,79.4,14.4,37.4,99.9,117.6)

y <- c(10.3,12.3,32.2,51.1,73.6,98.6,104.4,103.5,120.2,110.6,119.1,123.7,127.5,120.8,122.9,
       18.4,12.4,47.7,51.1,89.8,99.9,104.1,115.8,118,112.5,110.2,120.2,126.8,126,
       14.1,23.6,41.4,76.4,97.6,109.8,101.6,115.9,119.6,114.2,121.5,120.7,129.6,123.8,
       17.3,23.4,43.2,77.5,92.3,109.2,104.8,112.6,118.7,113.5,121.3,122.7,122.6,124)

plot(x,y, pch = 19, xlim = c(0,140),ylim = c(0,140),cex = .5,
     main = 'Hickory Trees Within Boundary')
arrows(x0 = c(10,130,10,10),
       y0 = c(10,10,130,10),
       x1 = c(10,130,130,130),
       y1 = c(130,130,130,10),
       length = 0)

grays = rgb(red = 0:255/255, blue = 0:255/255, green = 0:255/255)
make_counts <- function(p = 2,make_plot = FALSE,...){
  size = 120/p
  bins <- centers_x <- centers_y <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      bins[i,j] <- length(which(x>((j-1)*size+10) & x<(j*size+10) & 
                                  y>((p+1-i-1)*size+10) & y<((p+1-i)*size)+10))
      centers_x[i,j] <- 10+((j-1)*size + j*size)/2
      centers_y[i,j] <- 10+((p+1-i-1)*size + (p+1-i)*size)/2
    }
  }
  
  if(make_plot){
    heatmap(bins,Rowv=NA,Colv=NA, 
              col=grays[length(grays):1], scale = "none",
            revC=TRUE, labRow = NA, labCol = NA,...)
  }
  return(list(counts = bins, centers_x = centers_x,centers_y = centers_y, size=size))
}

#dat_2 <- make_counts(p=2,make_plot=FALSE)
dat_4 <- make_counts(p=4,make_plot= TRUE,main = "Aggregated Counts", margins = c(1.5,1.5))

#heatmap.2(dat_4$counts,Rowv=FALSE,Colv=FALSE, dendrogram = 'none',
 #         col=grays[length(grays):1], scale = "none",
  #        revC=TRUE, labRow = NA, labCol = NA,symm = TRUE,
   #       density.info = 'none')

kernel_mat <- function(centers_x,centers_y,theta_2 = 1,print_dist = FALSE){
  
  ker_mat <- dist_mat <- matrix(0,dim(centers_x)[1]^2,dim(centers_x)[2]^2)
  v_c_x <- c(centers_x)
  v_c_y <- c(centers_y)
  for(i in 1:dim(ker_mat)[1]){
    for(j in 1:dim(ker_mat)[2]){
      (dist_mat[i,j] <- sqrt((v_c_x[i] - v_c_x[j])^2 +(v_c_y[i] - v_c_y[j])^2))
      #dist_mat[i,j] <- abs(v_c_x[i] - v_c_x[j]) + abs(v_c_y[i]- v_c_y[j])
      (ker_mat[i,j] <- exp(-(dist_mat[i,j]/exp(theta_2))^2))
    }
  }
  if(print_dist){
    print("Distances:")
    print(dist_mat)
  }
  return(ker_mat)
}


ker_test <- kernel_mat(dat_4$centers_x,dat_4$centers_y,theta_2=4)
ker_test

mu_theta <- c(-0.383,-5.19,4)
sig_theta <- c(1,1,1)

gibbs_hick <- function(iters = 1E4,burn = 0.3,data_list = dat_4,
                      theta_kap = rep(1,3)){
  
  data = data_list$counts
  centers_x = data_list$centers_x
  centers_y = data_list$centers_y
  sizes = rep(data_list$size,length(c(data)))
  #sizes <- rep(1,length(c(data)))
  burnin <- floor(iters*burn)
  theta_acc = 0
  
  #Set up containers
  (J <- prod(dim(data)))
  alpha_mat <- beta_mat <- gamma_rf_mat <- matrix(0,iters,J)
  theta_mat <- matrix(0,iters,3)
  
  
  
  #Get initial values from priors
  (theta_cur <- rnorm(3,mu_theta,sig_theta))
  (alpha_cur <- rep(exp(theta_cur[1]),J))
  (beta_cur <- rep(exp(-theta_cur[2]),J))
  
  ker_cur <- kernel_mat(centers_x,centers_y,theta_2 = theta_cur[3])
  
  (gamma_rf <- rgamma(J,alpha_cur,scale=beta_cur^(-1)))#Need to double-check beta vs beta^(-1)
  lambda <- sweep(ker_cur,2,gamma_rf,FUN="*")
  lambda_i_plus <- apply(lambda,1,sum)
  (p <- sweep(lambda,1,lambda_i_plus,FUN="/"))
  
  N_ij <- matrix(0,dim(ker_cur)[1],dim(ker_cur)[2])
  for(i in 1:dim(N_ij)[1]){
    N_ij[i,] <- rmultinom(1,size=c(data)[i],prob = p[i,])
  }
  
  
  for(t in 1:(iters+burnin)){
    
    #Print out every 10% progress
    if(!(t%%(0.1*iters))){
      flush.console()
      cat("\r t = ", t,', Theta acc rate:"',theta_acc/t)
    }
    
    
    #Update Gamma RF
    (alpha_t <- alpha_cur + apply(N_ij,2,sum))
    ker_plus_j_cur <- apply(sweep(ker_cur,1,sizes,FUN="*"),2,sum)
    (beta_t <- beta_cur + ker_plus_j_cur)
    
    (gamma_rf <- rgamma(J,alpha_t,scale=beta_t^(-1)))
    (lambda <- sweep(ker_cur,2,gamma_rf,FUN="*"))
    (lambda_i_plus <- apply(lambda,1,sum))
    (p <- sweep(lambda,1,lambda_i_plus,FUN="/"))
    
    #Update augmentation pts
    N_ij <- matrix(0,dim(ker_cur)[1],dim(ker_cur)[2])
    for(i in 1:dim(N_ij)[1]){
      N_ij[i,] <- rmultinom(1,size=c(data)[i],prob = p[i,])
    }
    
    #Update theta with MH step
    (theta_st <- rnorm(3,theta_cur,theta_kap))
    (ker_st <- kernel_mat(centers_x,centers_y,theta_st[3]))
    ker_plus_j_st <- apply(sweep(ker_st,1,sizes,FUN="*"),2,sum)
    
    (alpha_st <- rep(exp(theta_st[1]),J))
    (beta_st <- rep(exp(-theta_st[2]),J))
    
    long_ker_cur <- c(ker_cur)
    long_ker_st = c(ker_st)
    long_n = c(N_ij)
    
    (log_ratio <- sum(dnorm(theta_st,mu_theta,sig_theta,log = TRUE)) -
        sum(dnorm(theta_cur,mu_theta,sig_theta,log = TRUE)) +
        sum(long_n*(log(long_ker_st) - log(long_ker_cur))) +
        sum((alpha_st)*log(gamma_rf) +
              alpha_st*log(beta_st) + log(gamma(alpha_cur)) -
              
              (alpha_cur)*log(gamma_rf) -
              alpha_cur*log(beta_cur) - log(gamma(alpha_st))) +
        -1*sum((beta_st - beta_cur + ker_plus_j_st - ker_plus_j_cur)*gamma_rf))
    
    err_count = 0
    tryCatch(
      if(runif(1)<exp(log_ratio)){
        theta_cur <- theta_st
        alpha_cur <- rep(exp(theta_st[1]),J)
        beta_cur <- rep(exp(-theta_st[2]),J)
        ker_cur <- ker_st
        theta_acc <- theta_acc + 1
      },
      error = function(e){print(paste('caught an error:',e))
        err_count = err_count + 1})
    
    
    if(t > burnin){
      theta_mat[t-burnin,] <- theta_cur
      alpha_mat[t-burnin,] <- alpha_cur
      beta_mat[t-burnin,] <- beta_cur
      gamma_rf_mat[t-burnin,] <- gamma_rf
    }
    
    
  }#t loop
  print(err_count)
  return(list(theta = theta_mat,alpha = alpha_mat, beta = beta_mat,gamma = gamma_rf_mat))
}

lego <- gibbs_hick(iters = 1E4,
                   theta_kap = c(.2,.2,.2))

theta_out <- lego$theta
thin = 100
iters <- 1:dim(theta_out)[1]
theta_thin <- theta_out[!iters%%thin,]
acf(theta_thin[,3])

hist(theta_thin[,2], xlab = '',
     main = expression('Histogram for'~theta[0]),
     cex.axis = 1.5)


traceplot(as.mcmc(lego$theta))

