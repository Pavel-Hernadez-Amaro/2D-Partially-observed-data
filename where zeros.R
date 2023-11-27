###############
# posicion de los ceros en la base de kronecker cuando faltan filas o columnas enteras

M_1=magic(30)#matrix(1:30,nrow = 30)
M_2=magic(20)#matrix(1:20,nrow = 20)

n_1=15:30  ### menor que dim(aux)[1]==30
n_2=5:20  ### menor que dim(aux_2)[1]==20

# ESTOS SON CEROS DE FILAS Y COLUMNAS ENTERAS
M_1[n_1,]=0
M_2[n_2,]=0

M_k=kronecker(M_2,M_1)

NA_ind_1=NULL
for ( ind in n_1) {
for (ind_i in 1:dim(M_2)[1]) {
  NA_ind_1=c(NA_ind_1,dim(M_1)[1]*(ind_i-1)+ind)
}}

NA_ind_2=NULL
for (ind_j in n_2) {
  NA_ind_2=c(NA_ind_2, (dim(M_1)[1]*(ind_j-1)+1):(dim(M_1)[1]*(ind_j)))
  
}

NA_ind=sort(unique(c(NA_ind_1,NA_ind_2)))

all.equal(length(which(rowSums(M_k)==0)),length(NA_ind))

NA_ind=NULL
for (ind_i in 1:length(n_1)) {
  for (ind_j in 1:length(n_2)) {
    NA_ind=c(NA_ind,(dim(M_2)[1])*(ind_i-1)+ind_j)
    
  }}

I=diag(dim(B_kron)[1])
I[NA_ind,]=0

###############
###############
###############

# posicion de los ceros en la base de kronecker cuando faltan solo una parte de una fila o columna

M_1=matrix(1:6,nrow = 1)#magic(5)#matrix(1:30,nrow = 30)
M_2=matrix(1:10,nrow = 5)#magic(6)#matrix(1:20,nrow = 20)

# n_1=15:30  ### menor que dim(aux)[1]==30
# n_2=5:20  ### menor que dim(aux_2)[1]==20

observed_data=kronecker(M_2,M_1)
observed_data[5,c(1,4)]=NA
observed_data[1:3,6]=NA
observed_data[2,3]=NA

c1=30
c2=30

aux_model=bspline(1:6, 1-1e-04, 6+1e-04, c1-3, 3)$B
aux_2_model=bspline(1:10, 1-1e-04, 10+1e-04, c2-3, 3)$B

B_kron_model=(kronecker(aux_2_model,(aux_model)))

NA_ind=NULL
for (ind_i in 1:dim(observed_data)[1]) {
  for (ind_j in 1:dim(observed_data)[2]) {
    if(is.na(observed_data[ind_i,ind_j])){
      NA_ind=sort(c(NA_ind,(dim(observed_data)[1]*(ind_j-1))+ind_i))
    }}}

I=diag(dim(B_kron_model)[1])
I[NA_ind,]=NA

B_kron_model[NA_ind,]=NA

Real_X= matrix(coef_x %*% t(I%*%B_kron_model), nrow=5, ncol = 12)
