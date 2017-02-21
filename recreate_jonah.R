library(ggplot2)
library(ggExtra)
library(reshape2)
library(dplyr)

params <- c("Norm","alpha","tau","etas","kpi")

#Reading in main
cent_m <- as.numeric(read.table("Dropbox/Research/JETSCAPE/training_data/main/cent.dat"))
design_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/design.dat")
colnames(design_m) <- params
mult_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/mult.dat")
mult_err_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/mult_err.dat")
v2_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/v2.dat")
v2_err_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/v2_err.dat")
v3_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/v3.dat")
v3_err_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/v3_err.dat")
colnames(mult_m) <- colnames(mult_err_m) <- colnames(v2_m) <-
  colnames(v2_err_m) <- colnames(v3_m) <- colnames(v3_err_m) <- as.character(cent_m)

#Reading in validation
cent_v <- read.table("Dropbox/Research/JETSCAPE/training_data/validation/cent.dat")
design_v <- read.table("Dropbox/Research/JETSCAPE/training_data/validation/design.dat")
colnames(design_v) <- params
mult_v <- read.table("Dropbox/Research/JETSCAPE/training_data/validation/mult.dat")
mult_err_v <- read.table("Dropbox/Research/JETSCAPE/training_data/validation/mult_err.dat")
v2_v <- read.table("Dropbox/Research/JETSCAPE/training_data/validation/v2.dat")
v2_err_v <- read.table("Dropbox/Research/JETSCAPE/training_data/validation/v2_err.dat")
v3_v <- read.table("Dropbox/Research/JETSCAPE/training_data/validation/v3.dat")
v3_err_v <- read.table("Dropbox/Research/JETSCAPE/training_data/validation/v3_err.dat")
colnames(mult_v) <- colnames(mult_err_v) <- colnames(v2_v) <-
  colnames(v2_err_v) <- colnames(v3_v) <- colnames(v3_err_v) <- as.character(cent_v)

#Figure 2, design points
p <- ggplot(design_m,aes(x = etas,y = tau)) + geom_point() +
  xlab(expression(paste(eta,"/s"))) + ylab(expression(tau[0])) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 25,face = "bold"))
ggExtra::ggMarginal(p, type = "histogram",bins = 20)


#Figure 3, model outputs
#ALICE points
alice <- read.csv("Dropbox/Research/JETSCAPE/mtd-paper/data/exp/alice_means.csv")[,-1]
alice$cent <- as.factor(cent_m)

#Reshaping it
mult_m_melted <- melt(mult_m)
mult_m_melted$rowid <- 1:dim(mult_m)[1]

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

v2_m_melted <- melt(v2_m)
v2_m_melted$rowid <- 1:dim(v2_m)[1]

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

  
v3_m_melted <- melt(v3_m,aes(colour = Legend))
v3_m_melted$rowid <- 1:dim(v3_m)[1]
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
fig_4_data <- data.frame(mult_m = mult_m$`22.5`, v2_m= v2_m$`22.5`)[-c(1,4),]
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

mult_m_sqrt <- sqrt(mult_m)#because sqrt is more normal
alice$mult_sqrt <- sqrt(alice$mult) #because sqrt is more normal

alice_m <- c(alice$mult_sqrt,alice$v2,alice$v3)
weights_m <- c(rep(1.2,6),rep(1,6),rep(0.6,6))

Y <- data.frame(mult = mult_m_sqrt,v2 = v2_m, v3 = v3_m)

clean_Y <- function(Ydat,weight_vec=weights_m,alice_vec = alice_m,trim = FALSE){
  Y2 <- Ydat
  for(i in 1:dim(Ydat)[1]){
    for(j in 1:dim(Ydat)[2]){
      Y2[i,j] <- weight_vec[j]*Ydat[i,j]/alice_vec[j]
    }
  }
  if(trim==TRUE){
    Y2 <- Y2[-c(1,4),]
  }
  return(Y2)
}

Y2 <- clean_Y(Y,trim = TRUE)
Y3 <- Y2 - lapply(Y2,mean)

Y_svd <- svd(Y3)
m <- dim(Y3)[1]
Z <- sqrt(m)*as.matrix(Y3)%*%Y_svd$v[,1:5]#This is our observation data matrix

#Some PCA plots
eigs <- Y_svd$d^2
V_q <- cumsum(eigs)/sum(eigs)
plot(V_q[1:6],type = 'o',xlab = "Number of PC q", main = "Fraction of Variance Explained",
     ylab = "Explained variance V(q)",pch = 19)
legend('bottomright',"Glauber",pch = 19, lwd =1)

U <- Y_svd$v[,1:5]
rownames(U) <- c(rep("N_ch",6),rep("v2",6),rep("v3",6))
View(U) #Different from Fig. 6, but can't figure out why. May be normalization

#############
##Emulation##
#############
#Might need to scale linearly?

#Scaling design matrix
Norm_d <- with(design_m,((Norm-min(Norm))/(max(Norm) - min(Norm))))
alpha_d <- with(design_m,((alpha-min(alpha))/(max(alpha) - min(alpha))))
tau_d <- with(design_m,((tau-min(tau))/(max(tau) - min(tau))))
etas_d <- with(design_m,((etas-min(etas))/(max(etas) - min(etas))))
kpi_d <- with(design_m,((kpi-min(kpi))/(max(kpi) - min(kpi))))

#Scaling validation matrix
Norm_x <- with(design_v,((Norm-min(Norm))/(max(Norm) - min(Norm))))
alpha_x <- with(design_v,((alpha-min(alpha))/(max(alpha) - min(alpha))))
tau_x <- with(design_v,((tau-min(tau))/(max(tau) - min(tau))))
etas_x <- with(design_v,((etas-min(etas))/(max(etas) - min(etas))))
kpi_x <- with(design_v,((kpi-min(kpi))/(max(kpi) - min(kpi))))

design_m_scaled <- cbind(Norm_d,alpha_d,tau_d,etas_d,kpi_d)
design_v_scaled <- cbind(Norm_x,alpha_x,tau_x,etas_x,kpi_x)

#Running the GPs
D1 <- design_m_scaled[-c(1,4),]
D1 <- cbind(rep(1,dim(D1)[1]),D1)
yd1 <- Z[,1]
X1 <- design_v_scaled

start_t <- proc.time()
res1 <- GP_emulator(D=D1,yd = yd1,X=X1)
res2 <- GP_emulator(D=D1,yd = Z[,2],X=X1)
res3 <- GP_emulator(D=D1,yd = Z[,3],X=X1)
res4 <- GP_emulator(D=D1,yd = Z[,4],X=X1)
res5 <- GP_emulator(D=D1,yd = Z[,5],X=X1)
stop_t <- proc.time()
stop_t-start_t

start_t_m <- proc.time()
res1_m <- GP_emulator(D=D1,yd = yd1,X=X1,R_type = "matern")
res2_m <- GP_emulator(D=D1,yd = Z[,2],X=X1,R_type = "matern")
res3_m <- GP_emulator(D=D1,yd = Z[,3],X=X1,R_type = "matern")
res4_m <- GP_emulator(D=D1,yd = Z[,4],X=X1,R_type = "matern")
res5_m <- GP_emulator(D=D1,yd = Z[,5],X=X1,R_type = "matern")
stop_t_m <- proc.time()
stop_t_m-start_t_m

#Getting predictions and 95% confidence intervals
Z_st = cbind(res1$mu_predict,res2$mu_predict,res3$mu_predict,
             res4$mu_predict,res5$mu_predict)
Z_st_err <- cbind(diag(res1$cov_predict),diag(res2$cov_predict),diag(res3$cov_predict),
                  diag(res4$cov_predict),diag(res5$cov_predict))

Y_st = m^(-.5)*Z_st%*%t(U)
Y_st_err = m^(-1)*Z_st_err%*%t(U^2)

##Reverse steps to get on original scale
means <- as.numeric(lapply(Y2,mean))
mean_mat <- t(matrix(rep(means,dim(Y_st)[1]),ncol=dim(Y_st)[1]))
alice_mat <- t(matrix(rep(alice_m,dim(Y_st)[1]),ncol=dim(Y_st)[1]))
weights_mat <- t(matrix(rep(weights_m,dim(Y_st)[1]),ncol=dim(Y_st)[1]))

Y_st_mean <- Y_st + mean_mat

Y_st_c <- Y_st_mean*alice_mat/weights_mat
Y_st_err_c <- Y_st_err*alice_mat^2/weights_mat^2 #errors

Y_st_c[,1:6] <- Y_st_c[,1:6]^2 #Squaring mult
Y_st_err_c[,1:6] <- Y_st_err_c[,1:6]^4*(2+4*Y_st_c[,1:6]/Y_st_err_c[,1:6])

colnames(Y_st_c) <- colnames(Y_st_err_c)  <-  paste0(rep(c("mult_","v2_","v3_"),each=6),
                           rep(as.character(cent_m),3))
Yst_comp <- data.frame(Y_st_c) %>% 
  select(ends_with("_2.5"), ends_with("22.5"),ends_with("42.5"))

Yst_err_comp <- data.frame(Y_st_err_c) %>%
  select(ends_with("_2.5"), ends_with("22.5"),ends_with("42.5"))

#Concatenate responses
Y_v <- cbind(mult_v,v2_v,v3_v)
Y_v_err <- cbind(mult_err_v,v2_err_v,v3_err_v)#errors from the computer validation

names(Y_v) <- names(Y_v_err)  <- paste0(rep(c("mult_","v2_","v3_"),each=3),
                    rep(as.character(cent_v),3))



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