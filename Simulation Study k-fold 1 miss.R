## SIMULATION STUDY

options(error=NULL)
library(ks)
library(tictoc)
library(refund)
library(fda)
library(mgcv)
library(SOP)
library(expm)
library(writexl)
library(pracma)
library(splines)
library(rgl)
library(plot3D)
library(plotly)
library(devtools)
library(ReconstPoFD)
library(Tplyr)
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(MASS)
library(caret)
library(disdat)
library(dismo)

# DEFINING THE GENERAL PARAMETERS   

R=100
N=100
k=10

Rsq = 0.95

d1=8
d2=8

c1=15
c2=15

c1_beta=15
c2_beta=15

####### DEFINING SOME EMPTY LISTS, ARRAYS, VECTORS AND MATRICES

error_y_test=array(dim = c(k,N/k,R))

error_Beta=Error_y_mean=array(dim = R)

nu=response=err_mean=err=array(0,dim=c(R,N))

temp=rep(0,R)

miss_points=vector(mode="list", length = N)

A=array(0,dim = c(N, N*c1*c2, R))

A_train=array(0,dim = c(N*((k-1)/k), (N*((k-1)/k))*c1*c2, R))
A_test=array(0,dim = c(N/k, (N/k)*c1*c2, R))

X_hat=X_hat_train=X_hat_test=Real_X_train=Real_X_test=Real_X=vector(mode = "list", length = R)
missing_points_test=missing_points_train=missing_points=model=aux_B=Beta_2d=vector(mode = "list", length = R)

####### DEFINING THE OBSERVATION GRID

x_b=20 #round(runif(1,10,30))
y_b=20 #round(runif(1,10,20))

x_observations=seq(0,1,length.out=x_b)
y_observations=seq(0,1,length.out=y_b)

observation_points=expand.grid(x_observations,y_observations)

####### DEFINING THE PARAMETERS OF THE INNER PRODUCT   

Inner_matrix=array(dim = c(N*c1*c2,c1_beta*c2_beta))
Inner_matrix_train=array(dim = c((N*((k-1)/k))*c1*c2,c1_beta*c2_beta))
Inner_matrix_test=array(dim = c(N/k*c1*c2,c1_beta*c2_beta))

sub=500
n=2*sub # THIS IS THE NUMBER OF INTERVALS FOR THE SIMPSON METHOD
m=2*sub

W_delta=array(dim = (n+1)*(m+1))

simp_w=rep(1,n+1)
even=seq(2,n+1-1,2)
odd=seq(3,n+1-1,2)

simp_w[even]=2
simp_w[odd]=4

h=(x_observations[x_b]-x_observations[1])/n  # if y_b and y_a are functions of x
HX=(y_observations[y_b]-y_observations[1])/m # this line should go inside the next for loop (the int_i lopp or outer loop)

Sim_w_x=(h/3)*simp_w

x = seq(x_observations[1],x_observations[x_b],h)
y = seq(y_observations[1],y_observations[y_b],HX)

W_x=(HX/3)*Sim_w_x
W_x_even=2*W_x
W_x_odd=4*W_x

for (aux in 1:(m+1)) {
  
  # print(c("aux = ", aux))
  
  if (aux==1 || aux==(m+1)) {
    
    W_delta[((n+1)*(aux-1)+1):((n+1)*aux)]= W_x
    
  }else{
    
    if (aux%%2==0) {
      
      W_delta[((n+1)*(aux-1)+1):((n+1)*aux)]= W_x_even
      
    }else{
      
      W_delta[((n+1)*(aux-1)+1):((n+1)*aux)]= W_x_odd
    }}
  
}

W_delta_save=W_delta


####### DEFINING THE B-spline BASIS USED IN THE MODEL (only depend on the observation grid and the number
# of basis c) 

aux_model_k=bspline(x_observations, x_observations[1]-1e-04, x_observations[x_b]+1e-04, c1-3, 3) # x_observations and y_observation have to be the complete domain
aux_2_model_k=bspline(y_observations, y_observations[1]-1e-04, y_observations[y_b]+1e-04, c2-3, 3) # this domain may not by fulfill by any of the surfaces 

knots_x=aux_model_k$knots
knots_y=aux_2_model_k$knots

fx=spline.des(knots_x, x, 3+1, 0*x)$design
fy=spline.des(knots_y, y, 3+1, 0*y)$design

aux_model=aux_model_k$B
aux_2_model=aux_2_model_k$B

B_kron_model=(kronecker(aux_2_model,(aux_model)))

#

aux_model_beta=bspline(x_observations, x_observations[1]-1e-04, x_observations[x_b]+1e-04, c1_beta-3, 3) # x_observations and y_observation have to be the complete domain
aux_2_model_beta=bspline(y_observations, y_observations[1]-1e-04, y_observations[y_b]+1e-04, c2_beta-3, 3) # this domain may not by fulfill by any of the surfaces 

knots_x_beta=aux_model_beta$knots
knots_y_beta=aux_2_model_beta$knots

fx_beta=spline.des(knots_x_beta, x, 3+1, 0*x)$design
fy_beta=spline.des(knots_y_beta, y, 3+1, 0*y)$design

#########

for (ind in 1:R) {
  
  set.seed(1000+ind)
  
  print(c("ind = ", ind))
  
  ######## GENERATING THE TWO DIMENSIONAL BASIS USING bifd 
  
  # Beta_2d[[ind]]=Beta_fun(x_observations, y_observations)
  
  # Beta_2d[[ind]]=Beta_H_exp(x_observations, y_observations)
  
  Beta_2d[[ind]]=Beta_H_saddle(x_observations, y_observations)
  
  ######## GENERATING THE TWO DIMENSIONAL BASIS USING bifd 
  
  Real_X[[ind]]=vector(mode = "list", length = N)
  
  # breaks_x=bspline(x_observations, x_observations[1]-1e-04, x_observations[x_b]+1e-04, d1-3, 3)$knots
  # breaks_y=bspline(y_observations, y_observations[1]-1e-04, y_observations[y_b]+1e-04, d2-3, 3)$knots
  # 
  # dife=diff(breaks_x)[1]
  # breaks_x=c(breaks_x[1]-dife,breaks_x, breaks_x[length(breaks_x)]+dife)
  # breaks_x=c(breaks_x[1]-dife,breaks_x, breaks_x[length(breaks_x)]+dife)
  # n_x=length(breaks_x)
  # 
  # dife=diff(breaks_y)[1]
  # breaks_y=c(breaks_y[1]-dife,breaks_y, breaks_y[length(breaks_y)]+dife)
  # breaks_y=c(breaks_y[1]-dife,breaks_y, breaks_y[length(breaks_y)]+dife)
  # n_y=length(breaks_y)
  # 
  # x_basis=create.bspline.basis(breaks=breaks_x,norder=bdeg[1]+1,dropin=c(1:5,(n_x-2):(n_x+2)))
  # y_basis=create.bspline.basis(breaks=breaks_y,norder=bdeg[1]+1,dropin=c(1:5,(n_y-2):(n_y+2)))
  
  missing_points[[ind]]=vector(mode = "list", length = N)
  
  for (j in 1:N) {
    
    # coef_xy=matrix(0,nrow = x_basis$nbasis, ncol=y_basis$nbasis)
    # 
    # coef_short=matrix(rnorm(d1*d2), nrow = d1, ncol = d2)
    # 
    # coef_xy[(length(x_basis$dropind)/2)+1:d1,(length(y_basis$dropind)/2)+1:d2] =coef_short
    # 
    # X_bifd=bifd(coef_xy, x_basis, y_basis)
    # 
    # X_bifd$coefs=coef_short
    # 
    # Real_X[[ind]][[j]]=eval.bifd(x_observations, y_observations, X_bifd) # sbasismat %*% coef %*% t(tbasismat)
    # 
    # nu[ind,j] = response_int_bifd(X_bifd,Beta_H_saddle, x_observations, y_observations)
    
    
    Real_X[[ind]][[j]]=Data_H(x_observations, y_observations)$DATA_N # GENERATING THE TWO DIMENSIONAL BASIS USING
    # Harold's THESIS
    
    # nu[ind,j]= response_int_H(Data_H,Beta_H_exp, x_observations, y_observations, N=FALSE)
    
    nu[ind,j]= response_int_H(Data_H,Beta_H_saddle, x_observations, y_observations, N=FALSE)
    
    # nu[ind,j]= response_int_H(Data_H,Beta_fun, x_observations, y_observations, N=FALSE)
    
    ######## MISSING FROM BORDERS AND INSIDE (SUB-grid VERSION)
    
    n_missing = 1 #sample(2:3,1)
    
    min_distance_x=9
    min_distance_y=9

    x_missing=sort(sample(1:(x_b-min_distance_x),n_missing))
    y_missing=sort(sample(1:(y_b-min_distance_y),n_missing))
    
    x_pos=x_missing:(x_missing+min_distance_x-1)
    y_pos=y_missing:(y_missing+min_distance_y-1)

    missing_points[[ind]][[j]]=aux=expand.grid(x_pos,y_pos)

    # boundries_x=matrix(nrow=n_missing, ncol=2)
    # boundries_y=matrix(nrow=n_missing, ncol=2)
    #
    # size_miss= 1 #round(runif(1,1,3)) #8

    # MORE THAN ONE MISSING SPOT

    # boundries_x[,1]=x_missing-(min_distance_x/size_miss)
    # boundries_x[,2]=x_missing+(min_distance_x/size_miss)
    # boundries_y[,1]=y_missing-(min_distance_y/size_miss)
    # boundries_y[,2]=y_missing+(min_distance_y/size_miss)

    # boundries_x[,2][which(boundries_x[,2]>x_observations[x_b])]=x_observations[x_b]
    # boundries_y[,2][which(boundries_y[,2]>y_observations[y_b])]=y_observations[y_b]
    #
    # boundries_x[,1][which(boundries_x[,1]<x_observations[1])]=x_observations[1]
    # boundries_y[,1][which(boundries_y[,1]<y_observations[1])]=y_observations[1]

    ###

    # if (which(x_missing==x_observations)-size_miss < 1) {
    #   miss_x_left=1
    # }else{
    #   miss_x_left=which(x_missing==x_observations)-size_miss
    # }
    #
    # if (which(x_missing==x_observations)+size_miss > x_b) {
    #   miss_x_right=x_b
    # }else{
    #   miss_x_right=which(x_missing==x_observations)+size_miss
    # }
    #
    # if (which(y_missing==y_observations)-size_miss < 1) {
    #   miss_y_left=1
    # }else{
    #   miss_y_left=which(y_missing==y_observations)-size_miss
    # }
    #
    # if (which(y_missing==y_observations)+size_miss > y_b) {
    #   miss_y_right=y_b
    # }else{
    #   miss_y_right=which(y_missing==y_observations)+size_miss
    # }


    # boundries_x[,1]=x_observations[miss_x_left]
    # boundries_x[,2]=x_observations[miss_x_right]
    # boundries_y[,1]=y_observations[miss_y_left]
    # boundries_y[,2]=y_observations[miss_y_right]
    #
    #
    # aux=NULL
    #
    # for (ind_miss in 1:n_missing) {
    #
    #   x_pos=which(x_observations>=boundries_x[ind_miss,1] & x_observations<=boundries_x[ind_miss,2])
    #   y_pos=which(y_observations>=boundries_y[ind_miss,1] & y_observations<=boundries_y[ind_miss,2])
    #
    #   aux=rbind(aux,expand.grid(x_pos,y_pos))
    #
    # }
    #
    # missing_points[[ind]][[j]]=unique(as.matrix(aux,ncol=2))
    # missing_points[[ind]][[j]]=missing_points[[ind]][[j]][order(missing_points[[ind]][[j]][,2]),]
    
    # missing_points[[ind]][[j]]=c(which(x_observations==x_missing),which(y_observations==y_missing))
    ###########
    
    #### FROM THIS POINT ON THE MISSING POINTS HAVE TO BE DEFINED
    
    for (add_miss in 1:dim(missing_points[[ind]][[j]])[1]) {

      Real_X[[ind]][[j]][missing_points[[ind]][[j]][add_miss,1],missing_points[[ind]][[j]][add_miss,2]]=NA
    }
    
    # Real_X[[ind]][[j]][missing_points[[ind]][[j]][1],missing_points[[ind]][[j]][2]]=NA
    
    miss_points[[j]]=vector(mode="list", length = ncol(Real_X[[ind]][[j]]))
    
    for (i in 1:ncol(Real_X[[ind]][[j]])) {
      
      miss_spots=NULL
      
      for (j_row in 1:nrow(Real_X[[ind]][[j]])) {
        
        if (is.na(Real_X[[ind]][[j]][j_row,i])) {
          
          miss_spots=c(miss_spots,j_row)
          
        }
        
      }
      
      if (!is.null(miss_spots)) {
        
        miss_points[[j]][[i]]=miss_spots
        
      }
      
    }
  }
  
  for (ind_miss_points in 1:N) {
    
    if (length(miss_points[[ind_miss_points]])!=ncol(Real_X[[ind]][[j]])) {
      stop(c("Revisa que el length de miss point en ", ind_miss_points, "sea", ncol(Real_X[[ind]][[j]])),call. = FALSE)
    }
    
  }
  
  var_e <-  (1/Rsq - 1) * var(nu[ind,]) #(1-Rsq)*var(nu[ind,])
  
  response[ind,]=nu[ind,]+rnorm(N,sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL #rnorm(nu[ind,])
  
  # response[ind,]=rbinom(N,1,exp(nu[ind,])/(1+exp(nu[ind,])))
  
  # y=rpois(N,exp(nu)) # POISSON MODEL
  
  ########### PRE-MULTIPLICATION OF THE KNRONECKER MATRIX BY AN IDENTITY MATRIX WITH ZEROS IN THE DIAG
  
  groups=kfold(Real_X[[ind]],k=k)
  
  tic()
  
  X_hat_train[[ind]]=vector(mode = "list", length = (N*((k-1)/k)))
  
  for (fold in 1:k) {
   
    print(c("fold = ", fold))
    
    test  = which(groups==fold)
    train = which(groups!=fold)

    response_train = response[train]
    response_test = response[test]
    
    Real_X_train[[ind]]=Real_X[[ind]][train]
    Real_X_test[[ind]]=Real_X[[ind]][test]
    
    missing_points_train[[ind]]=missing_points[[ind]][train]
    
    for (j in 1:(N*((k-1)/k))) {

    print(c("j = ", j))
    
    NA_ind=NULL
    for (ind_i in 1:dim(Real_X_train[[ind]][[j]])[1]) {
      for (ind_j in 1:dim(Real_X_train[[ind]][[j]])[2]) {
        if(is.na(Real_X_train[[ind]][[j]][ind_i,ind_j])){
          NA_ind=sort(c(NA_ind,(dim(Real_X_train[[ind]][[j]])[1]*(ind_j-1))+ind_i))
        }}}
    
    
    I=diag(dim(B_kron_model)[1])
    I[NA_ind,]=0
    
    aux_22=B2XZG(I%*%B_kron_model, pord = c(2,2), c = c(c1,c2))
    
    X_response_train=c(Real_X_train[[ind]][[j]])
    
    aux_33=XZG2theta(X = aux_22$X,Z = aux_22$Z,G = aux_22$G,T = aux_22$T,y= X_response_train)#, weights = w)
    
    A_train[j,((c2*c1*(j-1))+1):(j*c1*c2),ind]=aux_33$theta
    
    X_hat_train[[ind]][[j]]= matrix(t(aux_33$theta) %*% t(B_kron_model), nrow=x_b, ncol = y_b)
    
    for (add_miss in 1:dim(missing_points_train[[ind]][[j]])[1]) {
      X_hat_train[[ind]][[j]][missing_points_train[[ind]][[j]][add_miss,1],missing_points_train[[ind]][[j]][add_miss,2]]=NA
    }
    
    # X_hat_train[[ind]][[j]][missing_points_train[[ind]][[j]][1],missing_points_train[[ind]][[j]][2]]=NA
    
    # err[ind,j]=max(abs(Real_X_train[[ind]][[j]]-X_hat_train[[ind]][[j]]),na.rm=TRUE)
    # err_mean[ind,j]=mean(abs(Real_X_train[[ind]][[j]]-X_hat_train[[ind]][[j]]),na.rm=TRUE)
    
    # print(c(err[ind,j]))
    # print(c(err_mean[ind,j]))
    
    # HERE BEGINS THE DOUBLE INTEGRAL
    
    W_delta=W_delta_save
    
    for (int_i in 1:length(y)) {
      
      na_x=NULL
      
      # print(c("int_i = ", int_i))
      
      if (length(y)!=m+1) {
        stop("length(y) has to be equal to m+1",call.=FALSE)
      }
      
      prev_obs=max(which(y[int_i]>=y_observations))
      
      next_obs=min(which(y[int_i]<=y_observations))
      
      na_22=unique(append(miss_points[[j]][[prev_obs]],miss_points[[j]][[next_obs]]))
      
      observed_points=which(!seq(1:nrow(Real_X_train[[ind]][[j]])) %in% na_22)
      
      for (aux_o in observed_points) {
        
        aux_prev=aux_o-1
        aux_next=aux_o+1
        
        if ((aux_prev %in% na_22) && (aux_next %in% na_22)) {
          
          na_22=sort(c(na_22,aux_o))
          observed_points=observed_points[-which(aux_o==observed_points)]
        }
        
      }
      
      where_diff=which(diff(observed_points) !=1)
      
      if (!isempty(where_diff)) {
        
        observed_points_final=NULL
        
        observed_points_1=range(observed_points[1:where_diff[1]])
        
        if (length(where_diff)>1) {
          
          for (aux_diff in seq(length(where_diff)-1)) {
            
            observed_points_final=c(observed_points_final,
                                    range(observed_points[(where_diff[aux_diff]+1):where_diff[aux_diff+1]]))
            
            
          }}
        
        observed_points_last=range(observed_points[(where_diff[length(where_diff)]+1):length(observed_points)])
        
        observed_points_group=c(observed_points_1, observed_points_final, observed_points_last)
      }else{
        
        observed_points_group=range(observed_points)
        
      }
      
      
      ###########
      
      where_diff_na_22=which(diff(na_22) !=1)
      
      if (!isempty(where_diff_na_22)) {
        
        na_22_final=NULL
        
        na_22_1=range(na_22[1:(where_diff_na_22[1])])
        na_22_1[2]=na_22_1[2]+1
        
        if (na_22[1]!=1) {
          na_22_1[1]=na_22_1[1]-1
        }
        
        if (length(where_diff_na_22)>1) {
          
          for (aux_diff in seq(length(where_diff_na_22)-1)) {
            
            aux=range(na_22[(where_diff_na_22[aux_diff]+1):where_diff_na_22[aux_diff+1]])
            aux[1]=aux[1]-1
            aux[2]=aux[2]+1
            
            na_22_final=c(na_22_final,aux)
            
            
          }}
        
        na_22_last=c(na_22[length(na_22)]-1,na_22[length(na_22)])
        
        if (na_22[length(na_22)]!=nrow(Real_X_train[[ind]][[j]])) {
          na_22_last[length(na_22_last)]=na_22_last[length(na_22_last)]+1
        }
        
        na_22_group=c(na_22_1, na_22_final, na_22_last)
      }else{
        
        if (!isempty(na_22)) {
          
          na_22_group=range(na_22)
          
          if (na_22_group[2]!=nrow(Real_X_train[[ind]][[j]])) {
            na_22_group[2]=na_22_group[2]+1
          }
          
          if (na_22_group[1]!=1) {
            na_22_group[1]=na_22_group[1]-1
          }
          
        }else{
          
          na_22_group=NULL
        }
      }
      
      if (!is.null(na_22_group)) {
        
        for(tag_x in 1:(length(na_22_group)/2)){
          
          na_x=c(na_x,which(x>=x_observations[na_22_group[(2*tag_x-1)]] & x<=x_observations[na_22_group[2*tag_x]]))
        }
        
        W_delta[((n+1)*(int_i-1)+1):((n+1)*int_i)][na_x]=0
        
      }
      
    } # for in int_i
    
    W=invvec(W_delta, ncol = m+1, nrow = n+1)
    
    # check=10
    # all.equal(as.matrix(W_delta[((n+1)*(check-1)+1):(check*(n+1))]),as.matrix(W[,check]))
    
    aux_GLAM=RH(t(Rten2((fy_beta),fy)),RH(t(Rten2((fx_beta),fx)),W))
    dim(aux_GLAM)=c(c1,c1_beta,c2,c2_beta)
    aux_GLAM_apperm=matrix(aperm(aux_GLAM,c(1,3,2,4)),nrow = c1*c2)
    
    # all.equal(aux,t(aux_GLAM))
    
    # all.equal(Inner_matrix_WO[(c1*c2*(j-1)+1):(c1*c2*j),],(aux_apperm))
    # all.equal(Inner_matrix_WO[(c1*c2*(j-1)+1):(c1*c2*j),],(aux_GLAM_apperm))
    
    Inner_matrix_train[(c1*c2*(j-1)+1):(c1*c2*j),]=aux_GLAM_apperm
    
    
  } # for in j
  
  aux_B[[ind]]=B2XZG(A_train[,,ind]%*%Inner_matrix_train, pord = c(2,2), c = c(c1_beta,c2_beta))
  
  model[[ind]]=XZG2theta(X = aux_B[[ind]]$X,Z = aux_B[[ind]]$Z,G = aux_B[[ind]]$G, T = aux_B[[ind]]$T,y=response_train)#, family = binomial() )
  
  ########## CALCULATING THE PREDICCTION ERROR IN THE TEST PARTITION

  X_hat_test[[ind]]=vector(mode = "list", length = (N/k))
  
  missing_points_test[[ind]]=missing_points[[ind]][test]
  
  for (j_test in 1:(N/k)) {
    
    print(c("j_test = ", j_test))
    
    NA_ind=NULL
    for (ind_i in 1:dim(Real_X_test[[ind]][[j_test]])[1]) {
      for (ind_j in 1:dim(Real_X_test[[ind]][[j_test]])[2]) {
        if(is.na(Real_X_test[[ind]][[j_test]][ind_i,ind_j])){
          NA_ind=sort(c(NA_ind,(dim(Real_X_test[[ind]][[j_test]])[1]*(ind_j-1))+ind_i))
        }}}
    
    
    I=diag(dim(B_kron_model)[1])
    I[NA_ind,]=0
    
    aux_22=B2XZG(I%*%B_kron_model, pord = c(2,2), c = c(c1,c2))
    
    X_response_test=c(Real_X_test[[ind]][[j_test]])
    
    aux_33=XZG2theta(X = aux_22$X,Z = aux_22$Z,G = aux_22$G,T = aux_22$T,y= X_response_test)#, weights = w)
    
    A_test[j_test,((c2*c1*(j_test-1))+1):(j_test*c1*c2),ind]=aux_33$theta
    
    X_hat_test[[ind]][[j_test]]= matrix(t(aux_33$theta) %*% t(B_kron_model), nrow=x_b, ncol = y_b)
    
    for (add_miss in 1:dim(missing_points_test[[ind]][[j_test]])[1]) {
      X_hat_test[[ind]][[j_test]][missing_points_test[[ind]][[j_test]][add_miss,1],missing_points_test[[ind]][[j_test]][add_miss,2]]=NA
    }
    
    # X_hat_test[[ind]][[j_test]][missing_points_test[[ind]][[j_test]][1],missing_points_test[[ind]][[j_test]][2]]=NA
    
    # HERE BEGINS THE DOUBLE INTEGRAL
    
    for (int_i in 1:length(y)) {
      
      na_x=NULL
      
      # print(c("int_i = ", int_i))
      
      if (length(y)!=m+1) {
        stop("length(y) has to be equal to m+1",call.=FALSE)
      }
      
      prev_obs=max(which(y[int_i]>=y_observations))
      
      next_obs=min(which(y[int_i]<=y_observations))
      
      na_22=unique(append(miss_points[[j_test]][[prev_obs]],miss_points[[j_test]][[next_obs]]))
      
      observed_points=which(!seq(1:nrow(Real_X_test[[ind]][[j_test]])) %in% na_22)
      
      for (aux_o in observed_points) {
        
        aux_prev=aux_o-1
        aux_next=aux_o+1
        
        if ((aux_prev %in% na_22) && (aux_next %in% na_22)) {
          
          na_22=sort(c(na_22,aux_o))
          observed_points=observed_points[-which(aux_o==observed_points)]
        }
        
      }
      
      where_diff=which(diff(observed_points) !=1)
      
      if (!isempty(where_diff)) {
        
        observed_points_final=NULL
        
        observed_points_1=range(observed_points[1:where_diff[1]])
        
        if (length(where_diff)>1) {
          
          for (aux_diff in seq(length(where_diff)-1)) {
            
            observed_points_final=c(observed_points_final,
                                    range(observed_points[(where_diff[aux_diff]+1):where_diff[aux_diff+1]]))
            
            
          }}
        
        observed_points_last=range(observed_points[(where_diff[length(where_diff)]+1):length(observed_points)])
        
        observed_points_group=c(observed_points_1, observed_points_final, observed_points_last)
      }else{
        
        observed_points_group=range(observed_points)
        
      }
      
      
      ###########
      
      where_diff_na_22=which(diff(na_22) !=1)
      
      if (!isempty(where_diff_na_22)) {
        
        na_22_final=NULL
        
        na_22_1=range(na_22[1:(where_diff_na_22[1])])
        na_22_1[2]=na_22_1[2]+1
        
        if (na_22[1]!=1) {
          na_22_1[1]=na_22_1[1]-1
        }
        
        if (length(where_diff_na_22)>1) {
          
          for (aux_diff in seq(length(where_diff_na_22)-1)) {
            
            aux=range(na_22[(where_diff_na_22[aux_diff]+1):where_diff_na_22[aux_diff+1]])
            aux[1]=aux[1]-1
            aux[2]=aux[2]+1
            
            na_22_final=c(na_22_final,aux)
            
            
          }}
        
        na_22_last=c(na_22[length(na_22)]-1,na_22[length(na_22)])
        
        if (na_22[length(na_22)]!=nrow(Real_X_test[[ind]][[j_test]])) {
          na_22_last[length(na_22_last)]=na_22_last[length(na_22_last)]+1
        }
        
        na_22_group=c(na_22_1, na_22_final, na_22_last)
      }else{
        
        if (!isempty(na_22)) {
          
          na_22_group=range(na_22)
          
          if (na_22_group[2]!=nrow(Real_X_test[[ind]][[j_test]])) {
            na_22_group[2]=na_22_group[2]+1
          }
          
          if (na_22_group[1]!=1) {
            na_22_group[1]=na_22_group[1]-1
          }
          
        }else{
          
          na_22_group=NULL
        }
      }
      
      if (!is.null(na_22_group)) {
        
        for(tag_x in 1:(length(na_22_group)/2)){
          
          na_x=c(na_x,which(x>=x_observations[na_22_group[(2*tag_x-1)]] & x<=x_observations[na_22_group[2*tag_x]]))
        }
        
        W_delta[((n+1)*(int_i-1)+1):((n+1)*int_i)][na_x]=0
        
      }
      
    } # for in int_i
    
    W=invvec(W_delta, ncol = m+1, nrow = n+1)
 
    aux_GLAM=RH(t(Rten2((fy_beta),fy)),RH(t(Rten2((fx_beta),fx)),W))
    dim(aux_GLAM)=c(c1,c1_beta,c2,c2_beta)
    aux_GLAM_apperm=matrix(aperm(aux_GLAM,c(1,3,2,4)),nrow = c1*c2)
    
    Inner_matrix_test[(c1*c2*(j_test-1)+1):(c1*c2*j_test),]=aux_GLAM_apperm
    
    
  } # for in j_test
  
  y_hat_test=A_test[,,ind]%*%Inner_matrix_test%*%model[[ind]]$theta
  
  error_y_test[fold,,ind]=abs(y_hat_test-(response_test))
  
  } # for in fold
  
  temp[ind]=toc(echo=0)
  
  ##### CALCULATING THE INTEGRATED MEAN SQUARED ERROR FOR THE BETA ESTIMATOR

  print(c("Calculating Beta error for iteration number", ind))
  
  X_hat[[ind]]=vector(mode = "list", length = N)
  
    for (j in 1:N) {
      
      NA_ind=NULL
      for (ind_i in 1:dim(Real_X[[ind]][[j]])[1]) {
        for (ind_j in 1:dim(Real_X[[ind]][[j]])[2]) {
          if(is.na(Real_X[[ind]][[j]][ind_i,ind_j])){
            NA_ind=sort(c(NA_ind,(dim(Real_X[[ind]][[j]])[1]*(ind_j-1))+ind_i))
          }}}
      
      
      I=diag(dim(B_kron_model)[1])
      I[NA_ind,]=0
      
      aux_22=B2XZG(I%*%B_kron_model, pord = c(2,2), c = c(c1,c2))

      X_response=c(Real_X[[ind]][[j]])
      
      aux_33=XZG2theta(X = aux_22$X,Z = aux_22$Z,G = aux_22$G,T = aux_22$T,y= X_response)#, weights = w)
      
      A[j,((c2*c1*(j-1))+1):(j*c1*c2),ind]=aux_33$theta
      
      X_hat[[ind]][[j]]= matrix(t(aux_33$theta) %*% t(B_kron_model), nrow=x_b, ncol = y_b)
      
      for (add_miss in 1:dim(missing_points[[ind]][[j]])[1]) {
        X_hat[[ind]][[j]][missing_points[[ind]][[j]][add_miss,1],missing_points[[ind]][[j]][add_miss,2]]=NA
      }

      # X_hat[[ind]][[j]][missing_points[[ind]][[j]][1],missing_points[[ind]][[j]][2]]=NA
      
      # HERE BEGINS THE DOUBLE INTEGRAL
      
      W_delta=W_delta_save
      
      for (int_i in 1:length(y)) {
        
        na_x=NULL
   
        if (length(y)!=m+1) {
          stop("length(y) has to be equal to m+1",call.=FALSE)
        }
        
        prev_obs=max(which(y[int_i]>=y_observations))
        
        next_obs=min(which(y[int_i]<=y_observations))
        
        na_22=unique(append(miss_points[[j]][[prev_obs]],miss_points[[j]][[next_obs]]))
        
        observed_points=which(!seq(1:nrow(Real_X[[ind]][[j]])) %in% na_22)
        
        for (aux_o in observed_points) {
          
          aux_prev=aux_o-1
          aux_next=aux_o+1
          
          if ((aux_prev %in% na_22) && (aux_next %in% na_22)) {
            
            na_22=sort(c(na_22,aux_o))
            observed_points=observed_points[-which(aux_o==observed_points)]
          }
          
        }
        
        where_diff=which(diff(observed_points) !=1)
        
        if (!isempty(where_diff)) {
          
          observed_points_final=NULL
          
          observed_points_1=range(observed_points[1:where_diff[1]])
          
          if (length(where_diff)>1) {
            
            for (aux_diff in seq(length(where_diff)-1)) {
              
              observed_points_final=c(observed_points_final,
                                      range(observed_points[(where_diff[aux_diff]+1):where_diff[aux_diff+1]]))
              
              
            }}
          
          observed_points_last=range(observed_points[(where_diff[length(where_diff)]+1):length(observed_points)])
          
          observed_points_group=c(observed_points_1, observed_points_final, observed_points_last)
        }else{
          
          observed_points_group=range(observed_points)
          
        }
        
        
        ###########
        
        where_diff_na_22=which(diff(na_22) !=1)
        
        if (!isempty(where_diff_na_22)) {
          
          na_22_final=NULL
          
          na_22_1=range(na_22[1:(where_diff_na_22[1])])
          na_22_1[2]=na_22_1[2]+1
          
          if (na_22[1]!=1) {
            na_22_1[1]=na_22_1[1]-1
          }
          
          if (length(where_diff_na_22)>1) {
            
            for (aux_diff in seq(length(where_diff_na_22)-1)) {
              
              aux=range(na_22[(where_diff_na_22[aux_diff]+1):where_diff_na_22[aux_diff+1]])
              aux[1]=aux[1]-1
              aux[2]=aux[2]+1
              
              na_22_final=c(na_22_final,aux)
              
              
            }}
          
          na_22_last=c(na_22[length(na_22)]-1,na_22[length(na_22)])
          
          if (na_22[length(na_22)]!=nrow(Real_X[[ind]][[j]])) {
            na_22_last[length(na_22_last)]=na_22_last[length(na_22_last)]+1
          }
          
          na_22_group=c(na_22_1, na_22_final, na_22_last)
        }else{
          
          if (!isempty(na_22)) {
            
            na_22_group=range(na_22)
            
            if (na_22_group[2]!=nrow(Real_X[[ind]][[j]])) {
              na_22_group[2]=na_22_group[2]+1
            }
            
            if (na_22_group[1]!=1) {
              na_22_group[1]=na_22_group[1]-1
            }
            
          }else{
            
            na_22_group=NULL
          }
        }
        
        if (!is.null(na_22_group)) {
          
          for(tag_x in 1:(length(na_22_group)/2)){
            
            na_x=c(na_x,which(x>=x_observations[na_22_group[(2*tag_x-1)]] & x<=x_observations[na_22_group[2*tag_x]]))
          }
          
          W_delta[((n+1)*(int_i-1)+1):((n+1)*int_i)][na_x]=0
          
        }
        
      } # for in int_i
      
      W=invvec(W_delta, ncol = m+1, nrow = n+1)
  
      aux_GLAM=RH(t(Rten2((fy_beta),fy)),RH(t(Rten2((fx_beta),fx)),W))
      dim(aux_GLAM)=c(c1,c1_beta,c2,c2_beta)
      aux_GLAM_apperm=matrix(aperm(aux_GLAM,c(1,3,2,4)),nrow = c1*c2)
 
      Inner_matrix[(c1*c2*(j-1)+1):(c1*c2*j),]=aux_GLAM_apperm
      
      
    } # for in j for the Beta IMSE error
    
    aux_B_Beta=B2XZG(A[,,ind]%*%Inner_matrix, pord = c(2,2), c = c(c1_beta,c2_beta))
    
    model_Beta=XZG2theta(X = aux_B_Beta$X,Z = aux_B_Beta$Z,G = aux_B_Beta$G, T = aux_B_Beta$T,y=response[ind,])#, family = binomial() )
    
  
  sub_IMSE=50
  n_IMSE=2*sub_IMSE # THIS IS THE NUMBER OF INTERVALS FOR THE SIMPSON METHOD
  m_IMSE=2*sub_IMSE
  
  W_delta=array(dim = (n_IMSE+1)*(m_IMSE+1))
  
  simp_w=rep(1,n_IMSE+1)
  even=seq(2,n_IMSE+1-1,2)
  odd=seq(3,n_IMSE+1-1,2)
  
  simp_w[even]=2
  simp_w[odd]=4
  
  h=(x_observations[x_b]-x_observations[1])/n_IMSE  # if y_b and y_a are functions of x
  HX=(y_observations[y_b]-y_observations[1])/m_IMSE # this line should go inside the next for loop (the int_i lopp or outer loop)
  
  Sim_w_x=(h/3)*simp_w
  
  x_IMSE = seq(x_observations[1],x_observations[x_b],h)
  y_IMSE = seq(y_observations[1],y_observations[y_b],HX)
  
  W_x=(HX/3)*Sim_w_x
  W_x_even=2*W_x
  W_x_odd=4*W_x
  
  for (aux in 1:(m_IMSE+1)) {
    
    # print(c("aux = ", aux))
    
    if (aux==1 || aux==(m_IMSE+1)) {
      
      W_delta[((n_IMSE+1)*(aux-1)+1):((n_IMSE+1)*aux)]= W_x
      
    }else{
      
      if (aux%%2==0) {
        
        W_delta[((n_IMSE+1)*(aux-1)+1):((n_IMSE+1)*aux)]= W_x_even
        
      }else{
        
        W_delta[((n_IMSE+1)*(aux-1)+1):((n_IMSE+1)*aux)]= W_x_odd
      }}
    
  }
  
  fx_beta_IMSE=spline.des(knots_x_beta, x_IMSE, 3+1, 0*x_IMSE)$design
  fy_beta_IMSE=spline.des(knots_y_beta, y_IMSE, 3+1, 0*y_IMSE)$design
  
  B_kron_beta=kronecker(fy_beta_IMSE,fx_beta_IMSE)
  
  Beta_hat = B_kron_beta %*% model_Beta$theta #matrix(B_kron_beta %*% model[[ind]]$theta, nrow=length(x_IMSE), ncol = length(y_IMSE))
  
  Beta_true=vec(Beta_H_saddle(x_IMSE, y_IMSE))
  
  error_Beta[ind]=imse_beta_2d(Beta_hat, Beta_true, diag(W_delta), 1)
  
} # for in ind

for (r in 1:R) {

  Error_y_mean[r]=mean(rowSums(error_y_test[,,r])/(N/k))
    
}

mean(Error_y_mean)

mean(error_Beta)

case=2

plot_ly(z = (Real_X[[ind]][[case]]),type="surface")

plot_ly(z = X_hat[[ind]][[case]],type="surface")

B_kron_beta=kronecker(aux_2_model_beta$B,aux_model_beta$B)

# Beta_hat=matrix(B_kron_beta %*% model[[ind]]$theta, nrow=x_b, ncol = y_b)

Beta_hat=matrix(B_kron_beta %*% model_Beta$theta, nrow=x_b, ncol = y_b)

err_Beta=max(abs(Beta_2d[[ind]]-Beta_hat),na.rm=TRUE)
err_Beta_mean=mean(abs(Beta_2d[[ind]]-Beta_hat),na.rm=TRUE)

plot_ly(z = Beta_2d[[ind]],type="surface")
plot_ly(z = (Beta_hat) ,type="surface")

y_hat=A[,,ind]%*%Inner_matrix%*%model[[ind]]$theta

error_y_abs=abs(y_hat-t(response))
error_y_rel=abs((y_hat-t(response))/t(abs(response)))


