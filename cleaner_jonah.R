library(ggplot2)
library(ggExtra)
library(reshape2)
library(dplyr)

setwd("/Users/Jake")

params <- c("Norm","alpha","tau","etas","kpi")
design_m <- read.table("~/Dropbox/Research/JETSCAPE/training_data/main/design.dat")
design_v <- read.table("~/Dropbox/Research/JETSCAPE/training_data/validation/design.dat")
colnames(design_m) <- colnames(design_v) <- params


#Getting main output
cent_names_m <- as.character(read.table("~/Dropbox/Research/JETSCAPE/training_data/main/cent.dat"))
out_names <- c("mult","v2","v3")
out_type <- c("","_err")
path_m <- "Dropbox/Research/JETSCAPE/training_data/main/"

out_m <- lapply(seq_along(1:length(out_type)), function(j){  
  temp <- lapply(seq_along(1:length(out_names)),function(i){
    out <- read.table(paste0(path_m,out_names[i],out_type[j],".dat"))
    names(out) <- cent_names_m
    return(out)
  })
  names(temp) <- out_names
  return(temp)
})
names(out_m) <- c("mean","sd")

#Getting validation output
path_v <- "Dropbox/Research/JETSCAPE/training_data/validation/"
cent_names_v <- as.character(read.table("Dropbox/Research/JETSCAPE/training_data/validation/cent.dat"))

out_v <- lapply(seq_along(1:length(out_type)), function(j){  
  temp <- lapply(seq_along(1:length(out_names)),function(i){
    out <- read.table(paste0(path_v,out_names[i],out_type[j],".dat"))
    names(out) <- cent_names_v
    return(out)
  })
  return(temp)
})
names(out_v) <- c("mean","sd")

#ALICE points
alice <- read.csv("Dropbox/Research/JETSCAPE/mtd-paper/data/exp/alice_means.csv")[,-1]
alice$cent <- as.factor(cent_names_m)

###################
####Plotting#######
###################
#Figure 2, design points
p <- ggplot(design_m,aes(x = etas,y = tau)) + geom_point() +
  xlab(expression(paste(eta,"/s"))) + ylab(expression(tau[0])) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 25,face = "bold"))
ggExtra::ggMarginal(p, type = "histogram",bins = 20)


#Figure 3, model outputs


#Reshaping it
mult_m_melted <- melt(out_m$mean$mult)
mult_m_melted$rowid <- 1:dim(out_m$mean$mult)[1]

ggplot(mult_m_melted,aes(colour = Legend)) + 
  geom_line(aes(x=variable,y = value, group = factor(rowid),colour = "Computer Model"),size = .2) + 
  xlab("Centrality %") + ylab("Glauber") +
  ggtitle(expression("{"~N[ch]~"}")) + 
  geom_point(data = alice, aes(x = cent, y = mult,colour = "Experiment")) + 
  scale_colour_manual(name = "Legend", values = c("gray","black"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid","blank"),
                        shape = c(NA,16)))) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 25,face = "bold"),
        legend.position="none")

v2_m_melted <- melt(out_m$mean$v2)
v2_m_melted$rowid <- 1:dim(out_m$mean$v2)[1]

ggplot(v2_m_melted) + 
  geom_line(aes(x=variable,y = value, group = factor(rowid),colour = "Computer Model"),size = .2) + 
  xlab("Centrality %") + ylab("Glauber") +
  ggtitle(expression(v[2]~"{2}")) +
  geom_point(data = alice, aes(x = cent, y = v2,colour = "Experiment")) + 
  scale_colour_manual(name = "Legend",values = c("gray","black"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid","blank"),
                        shape = c(NA,16)))) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 25,face = "bold"),
        legend.position="none")


v3_m_melted <- melt(out_m$mean$v3,aes(colour = Legend))
v3_m_melted$rowid <- 1:dim(out_m$mean$v3)[1]
ggplot(v3_m_melted) + 
  geom_line(aes(x=variable,y = value, group = factor(rowid),colour = "Computer Model"),size = .2) +
  xlab("Centrality %") + ylab("Glauber") +
  ggtitle(expression(v[3]~"{2}")) +
  geom_point(data = alice, aes(x = cent, y = v3,colour = "Experiment")) + 
  scale_colour_manual(name = "Legend",values = c("gray","black"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid","blank"),
                        shape = c(NA,16)))) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 25,face = "bold"))


#Fig 4 - comparison of two outputs
fig_4_data <- data.frame(mult_m = out_m$mean$mult$`22.5`, v2_m= out_m$mean$v2$`22.5`)[-c(1,4),]
p <- ggplot(fig_4_data, aes(x = sqrt(mult_m),y = v2_m)) + geom_point() + 
  xlab(expression(sqrt(N[ch]))) + ylab(expression(v[2]~"{2}")) + 
  annotate("text",x = 45,y = 0.04,label = "Glauber 20-25%",size = 6)
ggExtra::ggMarginal(p, type = "histogram",bins = 20)

#####################
#Principal components
#####################
#Center and standardize
##Note: Jonah instructs to "Divide each observable by its corresponding experimental
##      value from ALICE. This converts evertying to unitless quatities of order one."

#Cleaning observables, as specified on page 6 of Jonah's paper

mult_m_sqrt <- sqrt(out_m$mean$mult)#because sqrt is more normal
Y <- data.frame(mult = mult_m_sqrt,v2 = out_m$mean$v2, v3 = out_m$mean$v3)

alice$mult_sqrt <- sqrt(alice$mult) #because sqrt is more normal
alice_m <- c(alice$mult_sqrt,alice$v2,alice$v3)
alice_mat <- t(matrix(rep(alice_m,dim(Y)[1]),ncol=dim(Y)[1]))

weights_m <- c(rep(1.2,6),rep(1,6),rep(0.6,6))
weights_mat <- t(matrix(rep(weights_m,dim(Y)[1]),ncol=dim(Y)[1]))

# #Multiplying by weights, standardizing by ALICE values
Y2 = Y*weights_mat/alice_mat
Y2 <- Y2[-c(1,4),] #design points removed

Y3 <- (Y2 - lapply(Y2,mean))
#Y3 <- (Y[-c(1,4),] - lapply(Y[-c(1,4),],mean))/lapply(Y[-c(1,4),],sd)

Y_svd <- svd(Y3)
m <- dim(Y3)[1]
Z_reg <- sqrt(m)*as.matrix(Y3)%*%Y_svd$v[,1:5]#This is our observation data matrix
Z <- sqrt(m)*Y_svd$u[,1:5] # what jonah trains on
#Some PCA plots
eigs <- Y_svd$d^2
V_q <- cumsum(eigs)/sum(eigs)
plot(V_q[1:6],type = 'o',xlab = "Number of PC q", main = "Fraction of Variance Explained",
     ylab = "Explained variance V(q)",pch = 19)
legend('bottomright',"Glauber",pch = 19, lwd =1)

V = Y_svd$v[,1:5]
S = diag(Y_svd$d[1:5])

#############
##Emulation##
#############
#Might need to scale linearly?

#Scaling design matrix
scaled_m = lapply(design_m,function(x){
  (x-min(x))/(max(x) - min(x))
})

design_m_scaled = do.call(cbind,scaled_m)

scaled_v = lapply(design_v,function(x){
  (x-min(x))/(max(x) - min(x))
})
design_v_scaled = do.call(cbind,scaled_v)


#Running the GPs
D1 <- design_m_scaled[-c(1,4),]
D1 <- cbind(rep(1,dim(D1)[1]),D1)
yd1 <- Z[,1]
X1 <- design_v_scaled

#Mengyang's Robust GaSP
library(RobustGaSP)

run_emulator <- function(author = "Mengyang",Z,D,X,...){
  start_t <- proc.time()
  if(author=="Mengyang"){
    print("Performing Mengyang's Robust Emulation")
    rob_mod <- lapply(data.frame(Z),rgasp,design = D[,-1],...)
    rob_pred <- lapply(rob_mod,RobustGaSP::predict,testing_input = X)
    Z_st = do.call(cbind,lapply(rob_pred,function(x){x$mean}))
    Z_st_err = do.call(cbind,lapply(rob_pred,function(x){x$sd}))
    time_len = proc.time()-start_t
    return(list(rob_mod = rob_mod, Z_st=Z_st,Z_st_err = Z_st_err,
                time = time_len[["elapsed"]]/60))
  }else if(author=="me"){
    print("Performing my slow emulation")
    res <- lapply(X=data.frame(Z),FUN=GP_emulator,D = D,pred_X = X,...)
    Z_st = do.call(cbind,lapply(res,function(x){x$mu_predict}))
    Z_st_err = do.call(cbind,lapply(res,function(x){diag(x$cov_predict)}))
    ells = do.call(cbind,lapply(res,function(x){x$ell_hat}))
    time_len = proc.time()-start_t
    return(list(Z_st=Z_st,Z_st_err = Z_st_err,time = time_len[["elapsed"]]/60,ells = ells))
  }else{
    stop("Wrong authors")
  }
}

run_robust <- run_emulator(author ="Mengyang",Z = Z_reg,D = D1, X = X1, kernel_type = "pow_exp",
                           alpha = rep(2,dim(as.matrix(D1))[2]))
run_me <- run_emulator("me",Z = Z, D = D1, X = X1, R_type = "exponential",alp = 2)


#Getting predictions and 95% confidence intervals
run <- run_robust #which one?
#Y_st = m^(-.5)*run$Z_st%*%S%*%t(V)
Y_st = m^(-.5)*run$Z_st%*%t(V)
#The variance is element-wise b/c we're only taking the diagonal,
##and treating that as how we calculate the 95% interval
#See meeting_6_29
#Y_st_err = m^(-1)*run$Z_st_err%*%S^2%*%t(V^2)
Y_st_err = m^(-1)*run$Z_st_err%*%t(V^2)


##################
##REGULAR CENTER AND SCALE
#######################
#  Y_st_c <- data.frame(Y_st)*lapply(Y[-c(1,4),],sd) + lapply(Y[-c(1,4),],mean)
#  Y_st_err_c <- data.frame(Y_st_err)*lapply(Y[-c(1,4),],var)
#####################
#Reverse steps to get on original scale
means <- as.numeric(lapply(Y2,mean))
mean_mat <- t(matrix(rep(means,dim(Y_st)[1]),ncol=dim(Y_st)[1]))


Y_st_mean <- Y_st + mean_mat

alice_mat_st <- t(matrix(rep(alice_m,dim(Y_st)[1]),ncol=dim(Y_st)[1]))
weights_mat_st <- t(matrix(rep(weights_m,dim(Y_st)[1]),ncol=dim(Y_st)[1]))

Y_st_c <- Y_st_mean*alice_mat_st/weights_mat_st
Y_st_err_c <- Y_st_err*alice_mat_st^2/weights_mat_st^2 #errors
##########################

Y_st_c[,1:6] <- Y_st_c[,1:6]^2 #Squaring mult
Y_st_err_c[,1:6] <- Y_st_err_c[,1:6]^4*(2+4*Y_st_c[,1:6]/Y_st_err_c[,1:6])

colnames(Y_st_c) <- colnames(Y_st_err_c)  <-  paste0(rep(c("mult_","v2_","v3_"),each=6),
                                                     rep(as.character(cent_names_m),3))
Yst_comp <- data.frame(Y_st_c) %>% 
  dplyr::select(ends_with("_2.5"), ends_with("22.5"),ends_with("42.5"))

Yst_err_comp <- data.frame(Y_st_err_c) %>%
  dplyr::select(ends_with("_2.5"), ends_with("22.5"),ends_with("42.5"))

#Concatenate responses
Y_v <- do.call(cbind,out_v$mean)
Y_v_err <- do.call(cbind,out_v$sd)#errors from computer simulation

names(Y_v) <- names(Y_v_err)  <- paste0(rep(c("mult_","v2_","v3_"),each=3),
                                        rep(as.character(cent_names_v),3))



Y_graph = data.frame(mult_val = c(Y_v$mult_2.5,Y_v$mult_22.5,Y_v$mult_42.5),
                     mult_val_err = c(Y_v_err$mult_2.5,Y_v_err$mult_22.5,
                                      Y_v_err$mult_42.5),
                     
                     mult_pred = c(Yst_comp$mult_2.5,Yst_comp$mult_22.5,Yst_comp$mult_42.5),
                     mult_pred_err = c(Yst_err_comp$mult_2.5,Yst_err_comp$mult_22.5,
                                       Yst_err_comp$mult_42.5),
                     
                     v2_val = c(Y_v$v2_2.5,Y_v$v2_22.5,Y_v$v2_42.5),
                     v2_val_err = c(Y_v_err$v2_2.5,Y_v_err$v2_22.5,
                                    Y_v_err$v2_42.5),
                     v2_pred = c(Yst_comp$v2_2.5,Yst_comp$v2_22.5,Yst_comp$v2_42.5),
                     v2_pred_err = c(Yst_err_comp$v2_2.5,Yst_err_comp$v2_22.5,
                                     Yst_err_comp$v2_42.5),
                     
                     
                     v3_val = c(Y_v$v3_2.5,Y_v$v3_22.5,Y_v$v3_42.5),
                     v3_val_err = c(Y_v_err$v3_2.5,Y_v_err$v3_22.5,
                                    Y_v_err$v3_42.5),
                     v3_pred = c(Yst_comp$v3_2.5,Yst_comp$v3_22.5,Yst_comp$v3_42.5),
                     v3_pred_err = c(Yst_err_comp$v3_2.5,Yst_err_comp$v3_22.5,
                                     Yst_err_comp$v3_42.5),
                     
                     cent = c(rep("0-5%",dim(Yst_comp)[1]),rep("20-25%",dim(Yst_comp)[1]),
                              rep("40-45%",dim(Yst_comp)[1])))


ggplot(Y_graph, aes(x = mult_pred,y = mult_val)) + geom_point(aes(shape = cent,col = cent))+
  geom_abline() + 
  geom_pointrange(aes(ymin=mult_val-mult_val_err, ymax=mult_val+mult_val_err,
                      shape = cent,col = cent)) +
  geom_errorbarh(aes(xmin = mult_pred - 2*sqrt(mult_pred_err),
                     xmax = mult_pred + 2*sqrt(mult_pred_err),
                     col = cent))+
  ggtitle(expression("{"~N[ch]~"}")) +
  xlab("Emulator Prediction") + 
  ylab("Model Validation Output") + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 25,face = "bold"),
        legend.position="none")

ggplot(Y_graph, aes(x = v2_pred,y = v2_val)) + geom_point(aes(shape = cent,col = cent)) + 
  geom_abline()+
  geom_pointrange(aes(ymin=v2_val-v2_val_err, ymax=v2_val+v2_val_err,
                      shape = cent,col = cent)) +
  geom_errorbarh(aes(xmin = v2_pred - 2*sqrt(v2_pred_err),
                     xmax = v2_pred + 2*sqrt(v2_pred_err),
                     col = cent))+
  ggtitle(expression(v[2]~"{2}")) +
  xlab("Emulator Prediction") + 
  ylab("Model Validation Output") + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 25,face = "bold"),
        legend.position="none")


ggplot(Y_graph, aes(x = v3_pred,y = v3_val)) + geom_point(aes(shape = cent,col = cent)) + 
  geom_abline()+
  geom_pointrange(aes(ymin=v3_val-v3_val_err, ymax=v3_val+v3_val_err,
                      shape = cent,col = cent)) +
  geom_errorbarh(aes(xmin = v3_pred - 2*sqrt(v3_pred_err),
                     xmax = v3_pred + 2*sqrt(v3_pred_err),
                     col = cent))+
  ggtitle(expression(v[3]~"{2}")) +
  xlab("Emulator Prediction") + 
  ylab("Model Validation Output") + 
  scale_shape_discrete(name = "Centrality Bin") + 
  scale_colour_discrete(name = "Centrality Bin") + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 25,face = "bold"))

