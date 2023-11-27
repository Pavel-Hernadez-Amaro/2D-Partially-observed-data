Beta_fun=function(x_observations, y_observations){

  Beta_x=-1*(5-40*((x_observations/x_b)-0.5)^2)/50
  Beta_y=-1*(5-40*((y_observations/y_b)-0.5)^2)/50
  
  res=matrix(0, nrow = length(x_observations), ncol = length(y_observations))
  res=kronecker(t(Beta_y),Beta_x)
  
  # for (i in 1:length(x_observations)) {
  #   for (j in 1:length(y_observations)) {
  #    
  #     res[i,j]=sin(2*pi*(x_observations[i]^2+y_observations[j]^2))
  #      
  #   }}
  
  # sin(log((x_observations[i]^2)))+log(y_observations[j]^2)
  
  # sin(2*pi*(x_observations^2))+cos(2*pi*y_observations^2)
  # sin(2*pi*(x_observations[i]^2+y_observations[j]^2))
  
  # Beta_x[[ind]]=((10*x_observations/x_b)-5)/50
  
  # Beta_x[[ind]]=-1*(5-40*((x_observations/x_b)-0.5)^2)/50
  # Beta_y[[ind]]=-1*(5-40*((y_observations/y_b)-0.5)^2)/50
  # Beta_2d=kronecker(t(Beta_y[[ind]]),Beta_x[[ind]])
  
  
  
  return(res)
  
}

Beta_2d=Beta_fun(x_observations, y_observations)

plot_ly(z = Beta_2d,type="surface")

#########################

d1=10
d2=10

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

coef_xy=matrix(0,nrow = x_basis$nbasis, ncol=y_basis$nbasis)

coef_xy[(length(x_basis$dropind)/2)+1:d1,(length(y_basis$dropind)/2)+1:d2]=matrix(rnorm(d1*d2),nrow = d1, ncol = d2)


aux=bifd(coef_xy, x_basis, y_basis)

coef_xy=matrix(rnorm(d1*d2),nrow = d1, ncol = d2)
aux$coefs=coef_xy

Beta_spline=eval.bifd(x_observations, y_observations, aux) # sbasismat %*% coef %*% t(tbasismat)

plot_ly(z = Beta_spline,type="surface")

#######

Beta_H=function(x_observations,y_observations,case=1){
  
  aux_beta=matrix(nrow = length(x_observations), ncol = length(y_observations))
  
  for (i in 1:length(x_observations)) {
    for (h in 1:length(y_observations)) {
      
      if (case==2) {
        
        aux_beta[i,h] = 5*exp(-8*((x_observations[i]-0.75)^2 + (y_observations[h] -0.75)^2)) + 5*exp(-8*((x_observations[i]-0.1)^2 + (y_observations[h] -0.1)^2))
      }else{

        aux_beta[i,h] = (4*x_observations[i]-2)^3 -3*(4*x_observations[i]-2)*(4*y_observations[h]-2)^2
        
      }
      
      
      
    }
  }
  aux_beta
  }

# aux_beta=matrix(nrow = x_b, ncol = y_b)





