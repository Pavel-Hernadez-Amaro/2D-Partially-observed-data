######### JUST SOME ELEMENTS OF A ROW OR COLUMN ARE MISSING
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

R=1
N=100

err_mean=err=temp=array(0,dim=c(R,N))

Rsq=0.9

d1=8
d2=8

c1=20
c2=20
 
c1_beta=50
c2_beta=50

bdeg=c(3,3)

sub=1000
n=2*sub # THIS IS THE NUMBER OF INTERVALS FOR THE SIMPSON METHOD
m=2*sub

W_delta=array(dim = (n+1)*(m+1))
# WO=matrix(0,nrow = (n+1)*(m+1),ncol = (n+1)*(m+1))

aux_B_WO=model_WO=model=aux_B=B_kron_beta=list()

Inner_matrix=Inner_matrix_WO=array(dim = c(N*c1*c2,c1_beta*c2_beta))

na_x=na_y=I=X_hat=missing_points=list(list(),list())

Real_X=vector(mode = "list", length = R)

simp_w=rep(1,n+1)
even=seq(2,n+1-1,2)
odd=seq(3,n+1-1,2)

simp_w[even]=2
simp_w[odd]=4

Beta_x=Beta_y=Beta_2d=list()

# array(dim = c(1,y_b))  
nu=response=array(dim = c(R,N))

A=array(0,dim = c(N, N*c1*c2, R))

for (ind in 1:R) {
  
  set.seed(1000+ind)
  
  print(c("ind = ", ind))
  
  x_b=round(runif(1,10,30))
  y_b=round(runif(1,10,20))
  # n_points=50
  
  Real_X=lapply(Real_X, function(x) vector(mode = "list", length = N))
  
  # Real_X=lapply(Real_X, function(x) matrix(nrow = x_b, ncol = y_b))
  
  x_observations=seq(0,1,length.out=x_b)
  y_observations=seq(0,1,length.out=y_b)
  
  observation_points=expand.grid(x_observations,y_observations)
  
  h=(x_observations[x_b]-x_observations[1])/n  # if y_b and y_a are functions of x
  HX=(y_observations[y_b]-y_observations[1])/m # this line should go inside the next for loop (the int_i lopp or outer loop)
  
  # Sim_w_x=(h/3)*diag(simp_w)
  # Sim_w_y=(HX/3)*diag(simp_w)
  
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
  
  W=invvec(W_delta, ncol = m+1, nrow = n+1)
  
  # aux_GLAM=RH(t(Rten(),W))
  # dim(aux_GLAM)=c(c1,c1_beta,c2,c2_beta)
  # aux_GLAM_apperm=matrix(aperm(aux_GLAM,c(1,3,2,4)),nrow = c1*c2)

  ################## CREATING THE FUNCTIONAL COEFFICIENT

  # Beta_2d[[ind]]=Beta_fun(x_observations, y_observations)
  
  Beta_2d[[ind]]=Beta_H_exp(x_observations, y_observations)
  
  ######## GENERATING THE TWO DIMENSIONAL BASIS USING bifd 

  breaks_x=bspline(x_observations, x_observations[1]-1e-04, x_observations[x_b]+1e-04, d1-3, 3)$knots
  breaks_y=bspline(y_observations, y_observations[1]-1e-04, y_observations[y_b]+1e-04, d2-3, 3)$knots

  dife=diff(breaks_x)[1]
  breaks_x=c(breaks_x[1]-dife,breaks_x, breaks_x[length(breaks_x)]+dife)
  breaks_x=c(breaks_x[1]-dife,breaks_x, breaks_x[length(breaks_x)]+dife)
  n_x=length(breaks_x)

  dife=diff(breaks_y)[1]
  breaks_y=c(breaks_y[1]-dife,breaks_y, breaks_y[length(breaks_y)]+dife)
  breaks_y=c(breaks_y[1]-dife,breaks_y, breaks_y[length(breaks_y)]+dife)
  n_y=length(breaks_y)

  x_basis=create.bspline.basis(breaks=breaks_x,norder=bdeg[1]+1,dropin=c(1:5,(n_x-2):(n_x+2)))
  y_basis=create.bspline.basis(breaks=breaks_y,norder=bdeg[1]+1,dropin=c(1:5,(n_y-2):(n_y+2)))

  for (j in 1:N) {
    
    # coef_xy=matrix(0,nrow = x_basis$nbasis, ncol=y_basis$nbasis)
    # 
    # coef_short=matrix(rnorm(d1*d2), nrow = d1, ncol = d2)
    # 
    # coef_xy[(length(x_basis$dropind)/2)+1:d1,(length(y_basis$dropind)/2)+1:d2] =coef_short
    # 
    # X_bifd=bifd(coef_xy, x_basis, y_basis)
    # 
    # # coef_xy=matrix(rnorm(d1*d2),nrow = d1, ncol = d2)
    # X_bifd$coefs=coef_short
    # 
    # Real_X[[ind]][[j]]=eval.bifd(x_observations, y_observations, X_bifd) # sbasismat %*% coef %*% t(tbasismat)
    # 
    # nu[ind,j]= response_int_bifd(X_bifd,Beta_H_exp, x_observations, y_observations)
    
    Real_X[[ind]][[j]]=Data_H(x_observations, y_observations) # GENERATING THE TWO DIMENSIONAL BASIS USING
                                                              # Harold's THESIS

    nu[ind,j]= response_int_H(Data_H,Beta_H_exp, x_observations, y_observations)
      
      
  }
  
  var_e <-  (1/Rsq - 1) * var(nu[ind,]) #(1-Rsq)*var(nu[ind,])
  
  response[ind,]=nu[ind,]+rnorm(N,sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL #rnorm(nu[ind,])
  
  # response[ind,]=rbinom(N,1,exp(nu[ind,])/(1+exp(nu[ind,])))
  
  # response[ind,]=rpois(N,exp(nu)) # POISSON MODEL
  
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
  
  for (j in 1:N) {
    
    # tic()
    
    print(c("j = ", j))
    
    aux_22=B2XZG(B_kron_model, pord = c(2,2), c = c(c1,c2))
    
    X_response=c(Real_X[[ind]][[j]])
    
    aux_33=XZG2theta(X = aux_22$X,Z = aux_22$Z,G = aux_22$G,T = aux_22$T,y= X_response)
    
    # temp_aux=toc(echo =0)
    # temp[index]=as.double(temp_aux$toc-temp_aux$tic)

    A[j,((c2*c1*(j-1))+1):(j*c1*c2),ind]=aux_33$theta
    
    X_hat[[ind]][[j]]= matrix(t(aux_33$theta) %*% t(B_kron_model), nrow=x_b, ncol = y_b)
    
    # for (add_miss in 1:dim(missing_points[[ind]][[j]])[1]) {
    #   X_hat[[ind]][[j]][missing_points[[ind]][[j]][add_miss,1],missing_points[[ind]][[j]][add_miss,2]]=NA
    # }
    
    
    # max(abs(coef_x-aux_33$theta))
    
    err[ind,j]=max(abs(Real_X[[ind]][[j]]-X_hat[[ind]][[j]]),na.rm=TRUE)
    err_mean[ind,j]=mean(abs(Real_X[[ind]][[j]]-X_hat[[ind]][[j]]),na.rm=TRUE)
    
    print(c(err[ind,j]))
    print(c(err_mean[ind,j]))
    
    # HERE BEGINS THE DOUBLE INTEGRAL
    
    # check=10
    # all.equal(as.matrix(W_delta[((n+1)*(check-1)+1):(check*(n+1))]),as.matrix(W[,check]))
    
    
    # aux_GLAM=RH(t(Rten(fy)),RH(t(Rten(fx)),W))
    aux_GLAM=RH(t(Rten2((fy_beta),fy)),RH(t(Rten2((fx_beta),fx)),W))
    dim(aux_GLAM)=c(c1,c1_beta,c2,c2_beta)
    aux_GLAM_apperm=matrix(aperm(aux_GLAM,c(1,3,2,4)),nrow = c1*c2)
    
    Inner_matrix[(c1*c2*(j-1)+1):(c1*c2*j),]=aux_GLAM_apperm
    
    
  } # for in j
  
  # I_x=kronecker(fy,fx)
  # I_beta=kronecker(fy_beta,fx_beta)
  # 
  # WO=diag(W_delta)
  # 
  # aux_product=t(I_x) %*% WO %*% I_beta

  # dim(aux_GLAM)
  # dim(aux_GLAM_apperm)
  # dim(aux_product)
  # all.equal(aux_GLAM_apperm, aux_GLAM)
  # all.equal(aux_GLAM_apperm, aux_product)

  # for (aux_j in 1:N) {
  #   Inner_matrix_WO[(c1*c2*(aux_j-1)+1):(c1*c2*aux_j),]=aux_product
  # }
  
  # all.equal(Inner_matrix_WO[(c1*c2*(j-1)+1):(c1*c2*j),],(aux_apperm))
  # all.equal(Inner_matrix_WO[(c1*c2*(j-1)+1):(c1*c2*j),],(aux_GLAM_apperm))
  
  aux_B[[ind]]=B2XZG(A[,,ind]%*%Inner_matrix, pord = c(2,2), c = c(c1_beta,c2_beta))
  
  model[[ind]]=XZG2theta(X = aux_B[[ind]]$X,Z = aux_B[[ind]]$Z,G = aux_B[[ind]]$G, T = aux_B[[ind]]$T,y=response[ind,], family = gaussian() )

  # aux_B_WO[[ind]]=B2XZG(A[,,ind]%*%Inner_matrix_WO, pord = c(2,2), c = c(c1_beta,c2_beta))
  # 
  # model_WO[[ind]]=XZG2theta(X = aux_B_WO[[ind]]$X,Z = aux_B_WO[[ind]]$Z,G = aux_B_WO[[ind]]$G, T = aux_B_WO[[ind]]$T,y=response[ind,], family = gaussian() )
  
  # temp[ind,j]=toc(echo=0)
  
  B_kron_beta[[ind]]=kronecker(aux_2_model_beta$B,aux_model_beta$B)
  
  
} # for in ind

# all.equal(Inner_matrix_WO,Inner_matrix)

case=13

plot_ly(z = (Real_X[[ind]][[case]]),type="surface")

plot_ly(z = X_hat[[ind]][[case]],type="surface")

Beta_hat=matrix(B_kron_beta[[ind]] %*% model[[ind]]$theta, nrow=x_b, ncol = y_b)
err_Beta=max(abs(Beta_2d[[ind]]-Beta_hat),na.rm=TRUE)
err_Beta_mean=mean(abs(Beta_2d[[ind]]-Beta_hat),na.rm=TRUE)
err_Beta_rel=max(abs(Beta_2d[[ind]]-Beta_hat)/abs(Beta_2d[[ind]]),na.rm=TRUE)

# Beta_hat_WO=matrix(B_kron_beta[[ind]] %*% model_WO[[ind]]$theta, nrow=x_b, ncol = y_b)
# err_Beta_WO=max(abs(Beta_2d[[ind]]-Beta_hat_WO),na.rm=TRUE)

plot_ly(z = Beta_2d[[ind]],type="surface")
plot_ly(z = (Beta_hat) ,type="surface")
# plot_ly(z = (Beta_hat_WO) ,type="surface")

y_hat=A[,,ind]%*%Inner_matrix%*%model[[ind]]$theta
 
# y_hat_WO=A[,,ind]%*%Inner_matrix_WO%*%model_WO[[ind]]$theta

# all.equal(y_hat_WO,y_hat)

# all.equal(model_WO[[ind]]$theta,model[[ind]]$theta)

error_y_abs=abs(y_hat-t(response))
error_y_rel=abs((y_hat-t(response))/t(abs(response)))

mean(error_y_abs)


mean(abs(t(nu)-t(response)))
