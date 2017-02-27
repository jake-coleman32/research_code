#LHC for Abhijit's code
library(lhs)
setwd("/Users/Jake/Dropbox/Research/JETSCAPE/JETSCAPE-STAT/")
set.seed(48)
qhats <- maximinLHS(10,1)
(qhats <- round(sort(qhats*4.75 + 0.25),2))
qhats <- c(0.25,qhats,5)
write.table(qhats, file = 'qhat_vals_100k.dat', row.names=FALSE, col.names = FALSE)

#(qhats <- c(0.25,read.table(file = 'qhat_vals.dat')[[1]],5))
#read.table(file = 'A_J_dist.dat')[,2]*100000

#Reading output
qs <- which(qhats<3.35)
qs <- 1:length(qhats)
paths <- paste0("q_new/q_",as.character(qhats[qs]),".dat")
dists <- lapply(seq_along(1:length(paths)),function(x){
  read.table(file=paths[x])
})

for(i in 1:length(paths)){
  plot(dists[[i]][,1],dists[[i]][,2],main = paste("Hist of",qhats[i]),
       xlab = "Aj",ylab = "Frequency")
}
