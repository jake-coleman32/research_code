
t = c(.11,.432,.754,1.077,1.399,1.721,2.043,2.366,2.688,3.010)
y = c(4.73,4.72,4.234,
      3.177,2.966,3.653,
      1.970,2.267,2.084,
      2.079,2.409,2.371,
      1.908,1.665,1.685,
      1.773,1.603,1.922,
      1.370,1.661,1.757,
      1.868,1.505,1.638,
      1.390,1.275,1.679,
      1.461,1.157,1.530)

t_rep <- rep(t,each = 3)

plot(t_rep,y,pch = 1, main = "MLE Fit Only",xlab = 't',
     ylab = 'Chemical Concentration', cex = .8,ylim = c(0,max(y)+.5))
x <- seq(0,3,length.out = 50)
lines(x,5*exp(-0.63*x), lwd = 2)


plot(t_rep,y,pch = 1, main = "Two-Model Calibration Toy",xlab = 't',
     ylab = 'Chemical Concentration',ylim = c(0,max(y)+.5))
x <- seq(0,3,length.out = 50)
lines(x,5*exp(-0.63*x), lwd = 2)
lines(x,2.7*exp(-0.24*x),col = 'red',lwd = 2)
lines(x,3.5*exp(-1.7*x)+1.5,col = 'blue',lwd = 2)
legend('topright',c("Model 1 (MLE)","Model 2"),lwd=2, col = c('black','red'))


