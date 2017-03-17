#Trying to show regular way fails for histograms


#Read in data
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

jet_path_comp <- "/home/grad/jrc71/Documents/Research/Computer_Emulation/JETSCAPE/JETSCAPE-STAT"
jet_path_lap <- "/Users/Jake/Dropbox/Research/JETSCAPE/JETSCAPE-STAT/"
hist_folder <- "q_dat_100k/"

current_path = jet_path_lap

q_vals_file <- paste0(current_path,"qhat_vals_100k.dat")

num_counts <- 100000
holdout <- 7
Y <- read_hist_data(hist_loc = paste0(current_path,hist_folder),
                    q_file = q_vals_file) *num_counts
Y_new <- Y[holdout,]
Y_test <- Y[-holdout,]
d <- read.table(q_vals_file)[-holdout,1]
d_new <-read.table(q_vals_file)[holdout,1]
(d_scaled <- (d-min(d))/(max(d)-min(d)))
d_new_s <- (d_new-min(d))/(max(d)-min(d))
t_star <- 0.5
trunc_cols <- which(as.numeric(colnames(Y))<t_star)
Y_trunc <- Y_test[,trunc_cols]
Y_new_trunc <- Y_new[trunc_cols]

Y_final <- scale(Y_trunc,center=TRUE,scale = TRUE)

I <- dim(Y_final)[1]
J <- dim(Y_final)[2]

#PCA on data
scaled_pca = TRUE

Y_svd <- svd(Y_final)

cumsum(Y_svd$d^2)/sum(Y_svd$d^2)
q <- 5 #Number of components, greater than 99%

V <- Y_svd$v[,1:q]
S <- diag(Y_svd$d[1:q])

Z <- Y_final%*%V
if(scaled_pca){
  Z1 <- Z%*%solve(S)*(sqrt(I-1))
}else{
  Z1 = Z
}
var(Z1[,1])

#Mengyang's package
#Running the GPs
D1 <- d_scaled
D1 <- cbind(rep(1,length(D1)),D1)
X1 <- matrix(d_new_s)

#Mengyang's Robust GaSP
library(RobustGaSP)

run_emulator <- function(Z,D,X,...){

  start_t <- proc.time()
  
  rob_mod <- lapply(data.frame(Z),rgasp,design = D[,-1],...)
  rob_pred <- lapply(rob_mod,RobustGaSP::predict,testing_input = X)
  Z_st = do.call(cbind,lapply(rob_pred,function(x){x$mean}))
  Z_st_err = do.call(cbind,lapply(rob_pred,function(x){x$sd}))
  time_len = proc.time()-start_t
  return(list(rob_mod = rob_mod, Z_st=Z_st,Z_st_err = Z_st_err,
              time = time_len[["elapsed"]]/60))
}

run <- run_emulator(Z=Z1,D=D1,X = X1)

if(scaled_pca){
  Y_st <- (I-1)^(-0.5)*run$Z_st%*%S%*%t(V)
  Y_st_err = (I-1)^(-1)*run$Z_st_err%*%S^2%*%t(V^2)
}else{
  Y_st <- run$Z_st%*%t(V)
  Y_st_err = run$Z_st_err%*%t(V^2)
}


Y_comp <- sweep(Y_st,2,apply(Y_trunc,2,sd),FUN = "*") %>%
  sweep(2,apply(Y_trunc,2,mean),FUN = "+")/num_counts

Y_comp_err <- sweep(Y_st_err,2,apply(Y_trunc,2,var),FUN = "*")/num_counts^2
                    

##################
##Summary Plots of Prediction
##################

meeting_parent <- '/Users/Jake/Dropbox/Research/Computer_Emulation/meetings/2017/'
meeting_folder <- 'meeting_3_16/'
path <- paste0(meeting_parent,meeting_folder)
save_pics = TRUE
suffix = '_strawman'


##Prediction Comparison
if(save_pics) pdf(paste0(path,'pred_comp',suffix,'.pdf'))
plot(as.numeric(colnames(Y_trunc)),Y_comp[1,],
     cex = 0.9, ylim = c(min(Y_comp[1,]-2*sqrt(Y_comp_err[1,])),
                         max(Y_comp[1,]+2*sqrt(Y_comp_err[1,]))),main = 'Predicted Bin Probabilities - Strawman',
     xlab = 'Bin',ylab = 'Predicted Probability',col = rgb(1,0,0,0.3),pch = 19,
     las = 1)

arrow_bottoms <- Y_comp[1,]-2*sqrt(Y_comp_err[1,])
arrow_tops <- Y_comp[1,]+2*sqrt(Y_comp_err[1,])

arrows(as.numeric(colnames(Y_trunc)), arrow_bottoms,
       as.numeric(colnames(Y_trunc)), arrow_tops,
       length=0.05, angle=90, code=3)
points(as.numeric(colnames(Y_trunc)),Y_new_trunc/num_counts,
       col = rgb(0,0,1,0.3),cex = 0.8,pch = 19)
legend('bottomright',c('Predicted','Truth'),col = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),pch = 19,
       pt.cex = c(0.9,0.8))
if(save_pics) dev.off()

#Prediction Residuals
if(save_pics) pdf(paste0(path,'pred_resid',suffix,'.pdf'))
plot(as.numeric(colnames(Y_trunc)),
     Y_comp - Y_new_trunc/num_counts,
     cex = 0.9, 
     ylim = c(min(arrow_bottoms - Y_new_trunc/num_counts),
              max(arrow_tops - Y_new_trunc/num_counts)),
     main = bquote(A[j]~'Bin Probability Residuals - Strawman'),
     xlab = bquote(A[j]~'Bin'),ylab = 'Residual',col = rgb(1,0,0,0.7),pch = 19,
     las = 1, yaxt = 'n')
axis(2, at=c(-0.004,0,0.004),las = 1,labels = FALSE)
text(y=c(-0.004,0,0.004),par("usr")[1],labels = c(-0.004,0,0.004),pos = 2,xpd = TRUE)

arrows(as.numeric(colnames(Y_trunc)),
       arrow_bottoms - Y_new_trunc/num_counts,
       as.numeric(colnames(Y_trunc)),
       arrow_tops - Y_new_trunc/num_counts, length=0.05, angle=90, code=3)
abline(h = 0, col = rgb(0,0,1,0.7))
if(save_pics) dev.off()



cols = rainbow(I)

#Visualizing prediction compared to other runs
plot(as.numeric(colnames(Y_trunc)),
     Y_trunc[1,]/num_counts,col = cols[1],type = 'l',
     xlab = 'Aj',
     ylab = 'Fraction of Bin',
     main = 'Data Across Input')
for(i in 2:I){
  lines(as.numeric(colnames(Y_trunc)),
        Y_trunc[i,]/num_counts,col = cols[i],type = 'l')
}

plot(as.numeric(colnames(Y_trunc)),Y_new_trunc/num_counts,lwd = 2,type = 'l')
lines(as.numeric(colnames(Y_trunc)),Y_comp, lwd = 2)

legend('bottomright','Out of Sample Histgram',lwd = 2)

