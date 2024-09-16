# Created 08-07-2024

# First analysis of the potential 2D data bases for the 2nd paper
options(error=NULL)
library(jpeg)
library(imager)
library(stringr)
library(pracma)
library(VDPO)
library(splines)
library(plotly)
library(Epi)
library(caret)
library(ks)


# This is for upload the response variable y={0,1}
data=read.csv("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/2do paper/codigo/Paper 2 Classification/2D/2D-Partially-observed-data/Data sets/archive/Dataset_for_AQI_Classification/Dataset_for_AQI_Classification/Dimapur_AQI_All_Info.csv")
# 

N=length(data$Filename)

X=array(dim = c(224,224,N))

both_path="C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/2do paper/codigo/Paper 2 Classification/2D/2D-Partially-observed-data/Dimapur/Both"

for (i in 1:N) {
  
  image=list.files(both_path,full.names = TRUE)[i]
  
  img <- readJPEG(image)
  
  bar <- img[,,1] + img[,,2] + img[,,3]  # Combine color channels
  bar <- bar / max(bar)  # Normalize the grayscale values
  
  X[,,i]=bar
  
}

# plot(1:2, type='n')
# rasterImage(X[1:70,1:200,350], 1.2, 1.27, 1.8, 1.73)
# # 
# plot_ly(z=X[1:50,,20], type="surface")

X_list <- lapply(1:dim(X)[3], function(i) X[1:70,1:200,i]) # 

data_pollution=data.frame(y=c(rep(0,table(data$AQI_Class)[1]),rep(1,table(data$AQI_Class)[2])))
data_pollution[["X"]]=X_list

miss_coord=c(21,63) # c(25,26) c(36,37) c(44,45) c(51,52) c(57,58)

data_miss = add_miss_2(data_pollution$X, n_missing = 1, min_distance_x = miss_coord[1], min_distance_y = miss_coord[2])

# case=20
# 
# plot_ly(z=X[1:60,1:220,case], type="surface")
# plot_ly(z=data_miss$X_miss[[case]], type="surface")

#add_miss_2(data_pollution$X, n_missing = 0, min_distance_x = 1, min_distance_y = 1)

data_pollution[["X_miss"]] <-  data_miss$X_miss
data_pollution[["miss_points"]] <-  data_miss$miss_points
data_pollution[["missing_points"]] <-  data_miss$missing_points

# all.equal(data_miss$X_miss,data_pollution$X) # SI NO HAY MISSING ESTO ES TRUE

aux=15

c1=aux
c2=aux
c1_beta=aux
c2_beta=aux

formula <- y ~ ffpo_2d(X = X, miss_points = miss_points, missing_points = missing_points, nbasis=c(c1,c2,c1_beta,c2_beta))
res <- VDPO(formula = formula, data = data_pollution, family = stats::binomial())
# model_B=ffpo_2d(X = data_pollution$X, miss_points = data_pollution$miss_points, missing_points = data_pollution$missing_points, nbasis=c(c1,c2,c1_beta,c2_beta))

# res$fit$fitted.values

miss_error=response=rep(0,N)

aux=ROC(as.vector(res$fit$fitted.values),as.vector(data_pollution$y),plot = NULL)

best_pos_12=which.max(aux$res[,1]+aux$res[,2]) # which.max(aux$res[,1]+aux$res[,2])
optimal_cut_12=aux$res[best_pos_12,5]

AUC=aux$AUC

for (i in 1:N) {
  
  if (res$fit$fitted.values[i]<optimal_cut_12) {
    
    response[i]=0
  }else{
    response[i]=1
  }
  
  miss_error[i]  = as.double(response[i]!=data_pollution$y[i])
  
}

sum(miss_error)

sum(miss_error)/N

miss_table=confusionMatrix(as.factor(response),as.factor(data_pollution$y))

# Beta_VDPO=invvec(model_B$Phi_ffpo2d %*% res$theta_ffpo2d)

# res$fit$fitted.values
