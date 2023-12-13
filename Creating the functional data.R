### DIFFERENT WAYS TO GENERATE THE FUNCTIONAL DATA (SURFACE)

d1=10

aux=bspline(x_observations, x_a-1e-04, x_b+1e-04, d1-3, 3)$B

d2=10

aux_2=bspline(y_observations, y_a-1e-04, y_b+1e-04, d2-3, 3)$B

B_kron=(kronecker(aux_2,(aux)))

coef_x=rnorm(d1*d2)

# Real_X[[ind]][[j]]= matrix(coef_x %*% t(B_kron), nrow=x_b, ncol = y_b)

##############

d1=21
d2=21

breaks_x=bspline(x_observations, x_a-1e-04, x_b+1e-04, d1-3, 3)$knots
breaks_y=bspline(y_observations, y_a-1e-04, y_b+1e-04, d2-3, 3)$knots

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

eval.bifd(x_observations, y_observations, aux) # sbasismat %*% coef %*% t(tbasismat)

##########

coef_x=runif(x_b*y_b,0.5,2)

Real_X[[ind]][[j]]=matrix(as.matrix(coef_x*(sin((observation_points[,1])^2/100)*cos((observation_points[,2])^2/100))),
                          nrow=x_b, ncol = y_b)


########################

Data_H=function(x_observations, y_observations, epsilon_1=0.2, epsilon_2=0.2, epsilon_data=0.005){

x_b=length(x_observations)
y_b=length(y_observations)

DATA=matrix(nrow = x_b, ncol=y_b)

a1=rnorm(1,0,epsilon_1)
a2=rnorm(1,0,epsilon_2)

for (i in 1:x_b) {
  for (j in 1:y_b) {

    
    DATA[i,j] = a1*cos(2*pi*x_observations[i]) + a2*cos(2*pi*y_observations[j]) + 1 + rnorm(1,0,epsilon_data)
  }
}

DATA

}

data_example=Data_H(x_observations, y_observations, epsilon_1 = 0.2, epsilon_2 = 0.2, epsilon_data = 0.005)

plot_ly(z = data_example, type="surface")

