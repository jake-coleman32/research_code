#Design matrix - main
params <- c("Norm","alpha","tau","etas","kpi")
design_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/design.dat")
colnames(design_m) <- params


cent_names <- as.character(read.table("Dropbox/Research/JETSCAPE/training_data/main/cent.dat"))
out_names <- c("mult","v2","v3")
out_type <- c("","_err")


path <- "Dropbox/Research/JETSCAPE/training_data/main/"

out_m <- sapply(out_type, function(x){  
  temp <- lapply(seq_along(1:length(out_names)),function(i){
    out <- read.table(paste0(path,out_names[i],x,".dat"))
    names(out) <- cent_names
    return(out)
  })
  return(temp)
})


out_m <- lapply(seq_along(1:length(out_m)),function(i){
  out <- read.table(paste0(path,out_names[i],".dat"))
  names(out) <- cent_names
  out
})



out_m_err <- lapply(seq_along(1:length(out_m_err)),function(i){
  out <- read.table(paste0(path,out_names[i],".dat"))
  names(out) <- cent_names
  out
})


mult_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/mult.dat")
mult_err_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/mult_err.dat")
v2_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/v2.dat")
v2_err_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/v2_err.dat")
v3_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/v3.dat")
v3_err_m <- read.table("Dropbox/Research/JETSCAPE/training_data/main/v3_err.dat")
colnames(mult_m) <- colnames(mult_err_m) <- colnames(v2_m) <-
  colnames(v2_err_m) <- colnames(v3_m) <- colnames(v3_err_m) <- as.character(cent_m)