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

Rsq = 0.95

d1=8
d2=8

c1=25
c2=25

c1_beta=50
c2_beta=50

sub=1000
n=2*sub # THIS IS THE NUMBER OF INTERVALS FOR THE SIMPSON METHOD
m=2*sub

W_delta=array(dim = (n+1)*(m+1))

miss_points=vector(mode="list", length = N)

model=aux_B=list()

Inner_matrix=array(dim = c(N*c1*c2,c1_beta*c2_beta))

na_x=na_y=Real_X=I=X_hat=missing_points=list(list(),list())

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
  
  # Beta_x[[ind]]=((10*x_observations/x_b)-5)/50
  # 
  # # Beta_x[[ind]]=-1*(5-40*((x_observations/x_b)-0.5)^2)/50
  # Beta_y[[ind]]=-1*(5-40*((y_observations/y_b)-0.5)^2)/50
  # #
  # Beta_2d[[ind]]=kronecker(t(Beta_y[[ind]]),Beta_x[[ind]])
  
  # Beta_2d[[ind]]=Beta_fun(x_observations, y_observations) #Beta_spline
  
  # Beta_2d[[ind]]=as.matrix(sin(2*pi*x_observations))%*%(2*y_observations-1)

  Beta_2d[[ind]]=Beta_H(x_observations, y_observations)
  
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
    
    coef_xy=matrix(0,nrow = x_basis$nbasis, ncol=y_basis$nbasis)
    
    coef_short=matrix(rnorm(d1*d2), nrow = d1, ncol = d2)
    
    coef_xy[(length(x_basis$dropind)/2)+1:d1,(length(y_basis$dropind)/2)+1:d2] =coef_short 
    
    X_bifd=bifd(coef_xy, x_basis, y_basis)
    
    # coef_xy=matrix(rnorm(d1*d2),nrow = d1, ncol = d2)
    X_bifd$coefs=coef_short
    
    Real_X[[ind]][[j]]=eval.bifd(x_observations, y_observations, X_bifd) # sbasismat %*% coef %*% t(tbasismat)
    
    # Real_X[[ind]][[j]]=matrix((as.matrix(coef_x*exp(-((observation_points[,1])/10)-((observation_points[,2])/10)))),
    #                           nrow=x_b, ncol = y_b)
    
    # nu[ind,j]=sum(Real_X[[ind]][[j]]*Beta_2d[[ind]])/(length(x_observations)+length(y_observations))
    
    nu[ind,j] = response_int_bifd(X_bifd,Beta_H, x_observations, y_observations)
    
    ######## MISSING FROM BORDERS AND INSIDE (CLOUD VERSION)
    
    n_missing=sample(2:8,1)
    
    x_missing=sort(sample(x_observations,n_missing))
    y_missing=sort(sample(y_observations,n_missing))
    
    miss_point_up=cbind(x_missing,y_missing)+min(x_observations[x_b],y_observations[y_b])/4
    miss_point_down=cbind(x_missing,y_missing)-min(x_observations[x_b],y_observations[y_b])/4
    
    miss_point_up[which(miss_point_up[,1]>x_observations[x_b]),1]=x_observations[x_b]
    miss_point_up[which(miss_point_up[,2]>y_observations[y_b]),2]=y_observations[y_b]
    
    miss_point_down[which(miss_point_down[,1]<x_observations[1]),1]=x_observations[1]
    miss_point_down[which(miss_point_down[,2]<y_observations[1]),2]=y_observations[1]
    
    aux=NULL
    
    for (ind_miss in 1:n_missing) {
      
      x_pos=which(x_observations>=miss_point_down[ind_miss,1] & x_observations<=miss_point_up[ind_miss,1])
      y_pos=which(y_observations>=miss_point_down[ind_miss,2] & y_observations<=miss_point_up[ind_miss,2])
      
      aux_x=sort(sample(x_pos,sample(2:length(x_pos),1)))
      aux_y=sort(sample(y_pos,sample(2:length(y_pos),1)))
      
      aux=rbind(aux,expand.grid(aux_x,aux_y))
      
    }
    
    missing_points[[ind]][[j]]=unique(as.matrix(aux,ncol=2))
    missing_points[[ind]][[j]]=missing_points[[ind]][[j]][order(missing_points[[ind]][[j]][,2]),]
    
    
    ###########
    
    #### FROM THIS POINT ON THE MISSING POINTS HAVE TO BE DEFINED
    
    for (add_miss in 1:dim(missing_points[[ind]][[j]])[1]) {
      
      Real_X[[ind]][[j]][missing_points[[ind]][[j]][add_miss,1],missing_points[[ind]][[j]][add_miss,2]]=NA
    }
    
    
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
  
  # if (length(miss_zones)==ncol(Real_X[[ind]][[j]])-1) {
  #   miss_zones[[ncol(Real_X[[ind]][[j]])]]=0
  #   miss_zones[[ncol(Real_X[[ind]][[j]])]]=list(NULL)
  # }
  
  var_e <-  (1/Rsq - 1) * var(nu[ind,]) #(1-Rsq)*var(nu[ind,])
  
  response[ind,]=nu[ind,]+rnorm(N,sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL #rnorm(nu[ind,])
  
  # response[ind,]=rbinom(N,1,exp(nu[ind,])/(1+exp(nu[ind,])))
  
  # y=rpois(N,exp(nu)) # POISSON MODEL
  
  ########### PRE-MULTIPLICATION OF THE KNRONECKER MATRIX BY AN IDENTITY MATRIX WITH ZEROS IN THE DIAG
  
  # c1=30
  # c2=30
  
  #
  
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
  
  for (j in 1:N) {
    
    tic()
    
    print(c("j = ", j))
    
    NA_ind=NULL
    for (ind_i in 1:dim(Real_X[[ind]][[j]])[1]) {
      for (ind_j in 1:dim(Real_X[[ind]][[j]])[2]) {
        if(is.na(Real_X[[ind]][[j]][ind_i,ind_j])){
          NA_ind=sort(c(NA_ind,(dim(Real_X[[ind]][[j]])[1]*(ind_j-1))+ind_i))
        }}}
    
    
    I[[ind]][[j]]=diag(dim(B_kron_model)[1])
    I[[ind]][[j]][NA_ind,]=0
    
    aux_22=B2XZG(I[[ind]][[j]]%*%B_kron_model, pord = c(2,2), c = c(c1,c2))
    
    w=rep(1,length(c(Real_X[[ind]][[j]])))
    w[is.na(c(Real_X[[ind]][[j]]))]=0
    X_response=c(Real_X[[ind]][[j]])
    X_response[is.na(X_response)]=0
    
    aux_33=XZG2theta(X = aux_22$X,Z = aux_22$Z,G = aux_22$G,T = aux_22$T,y= X_response, weights = w)
    
    # temp_aux=toc(echo =0)
    # temp[index]=as.double(temp_aux$toc-temp_aux$tic)
    
    temp[ind,j]=toc(echo=0)
    
    A[j,((c2*c1*(j-1))+1):(j*c1*c2),ind]=aux_33$theta
    
    X_hat[[ind]][[j]]= matrix(t(aux_33$theta) %*% t(B_kron_model), nrow=x_b, ncol = y_b)
    
    for (add_miss in 1:dim(missing_points[[ind]][[j]])[1]) {
      X_hat[[ind]][[j]][missing_points[[ind]][[j]][add_miss,1],missing_points[[ind]][[j]][add_miss,2]]=NA
    }
    
    
    # max(abs(coef_x-aux_33$theta))
    
    err[ind,j]=max(abs(Real_X[[ind]][[j]]-X_hat[[ind]][[j]]),na.rm=TRUE)
    err_mean[ind,j]=mean(abs(Real_X[[ind]][[j]]-X_hat[[ind]][[j]]),na.rm=TRUE)
    
    print(c(err[ind,j]))
    print(c(err_mean[ind,j]))
    
    # HERE BEGINS THE DOUBLE INTEGRAL
    
    for (int_i in 1:length(y)) {
      
      na_x=NULL
      
      # print(c("int_i = ", int_i))
      
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
    
    # check=10
    # all.equal(as.matrix(W_delta[((n+1)*(check-1)+1):(check*(n+1))]),as.matrix(W[,check]))
    

    # aux_GLAM=RH(t(Rten(fy)),RH(t(Rten(fx)),W))
    aux_GLAM=RH(t(Rten2((fy_beta),fy)),RH(t(Rten2((fx_beta),fx)),W))
    dim(aux_GLAM)=c(c1,c1_beta,c2,c2_beta)
    aux_GLAM_apperm=matrix(aperm(aux_GLAM,c(1,3,2,4)),nrow = c1*c2)
    
    # all.equal(aux,t(aux_GLAM))
    
    # all.equal(Inner_matrix_WO[(c1*c2*(j-1)+1):(c1*c2*j),],(aux_apperm))
    # all.equal(Inner_matrix_WO[(c1*c2*(j-1)+1):(c1*c2*j),],(aux_GLAM_apperm))
    
    Inner_matrix[(c1*c2*(j-1)+1):(c1*c2*j),]=aux_GLAM_apperm
    
      
  } # for in j
  
  aux_B[[ind]]=B2XZG(A[,,ind]%*%Inner_matrix, pord = c(2,2), c = c(c1_beta,c2_beta))
  
  model[[ind]]=XZG2theta(X = aux_B[[ind]]$X,Z = aux_B[[ind]]$Z,G = aux_B[[ind]]$G, T = aux_B[[ind]]$T,y=response[ind,])#, family = binomial() )
  
  
  
} # for in ind

case=19

plot_ly(z = (Real_X[[ind]][[case]]),type="surface")

plot_ly(z = X_hat[[ind]][[case]],type="surface")

B_kron_beta=kronecker(aux_2_model_beta$B,aux_model_beta$B)

Beta_hat=matrix(B_kron_beta %*% model[[ind]]$theta, nrow=x_b, ncol = y_b)

err_Beta=max(abs(Beta_2d[[ind]]-Beta_hat),na.rm=TRUE)
err_Beta_mean=mean(abs(Beta_2d[[ind]]-Beta_hat),na.rm=TRUE)

plot_ly(z = Beta_2d[[ind]],type="surface")
plot_ly(z = (Beta_hat) ,type="surface")

y_hat=A[,,ind]%*%Inner_matrix%*%model[[ind]]$theta

error_y_abs=abs(y_hat-t(response))
error_y_rel=abs((y_hat-t(response))/t(abs(response)))


