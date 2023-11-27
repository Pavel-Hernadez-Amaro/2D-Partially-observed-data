######### COMPLETE ROWS AND COLUMNS

aux_l_x=round(length(1:round(x_b/3))/4)
aux_m_x=round(length(round(x_b/3):round(2*x_b/3))/4)
aux_r_x=round(length(round(2*x_b/3):x_b)/4)

miss_left_x=sort(round(runif(aux_l_x,1,x_b/3)))
miss_mid_x=sort(round(runif(aux_m_x,x_b/3,2*x_b/3)))
miss_right_x=sort(round(runif(aux_r_x,2*x_b/3,x_b)))

aux_l_y=round(length(1:round(x_b/3))/4)
aux_m_y=round(length(round(x_b/3):round(2*x_b/3))/4)
aux_r_y=round(length(round(2*x_b/3):x_b)/4)

miss_left_y=sort(round(runif(aux_l_y,1,y_b/3)))
miss_mid_y=sort(round(runif(aux_m_y,y_b/3,2*y_b/3)))
miss_right_y=sort(round(runif(aux_r_y,2*y_b/3,y_b)))

n_1[[ind]][[j]]=unique(c(miss_left_x,miss_mid_x,miss_right_x))
n_2[[ind]][[j]]=unique(c(miss_left_y,miss_mid_y,miss_right_y))

Real_X[[ind]][[j]][n_1[[ind]][[j]],]=NA
Real_X[[ind]][[j]][,n_2[[ind]][[j]]]=NA

#######################

############ JUST END ROWS AND COLUMNS

n_1[[ind]][[j]]=sort(round(runif(1,x_b/5,x_b)):x_b)
n_2[[ind]][[j]]=sort(round(runif(1,y_b/5,y_b)):y_b)

#####

######### JUST IN THE MIDDLE

n_1=n_2=missing_points=list(list(),list())

n_missing=sample(2:8,1)

x_missing=sort(sample(1:x_b,n_missing))
y_missing=sort(sample(1:y_b,n_missing))

miss_point_up=cbind(x_missing,y_missing)+round(min(x_b,y_b)/4)
miss_point_down=cbind(x_missing,y_missing)-round(min(x_b,y_b)/4)

miss_point_up[which(miss_point_up[,1]>x_b),1]=x_b
miss_point_up[which(miss_point_up[,2]>y_b),2]=y_b

miss_point_down[which(miss_point_down[,1]<x_a),1]=x_a
miss_point_down[which(miss_point_down[,2]<y_a),2]=y_a

aux=NULL
for (ind_miss in 1:n_missing) {

  aux_x=sort(sample(miss_point_down[ind_miss,1]:miss_point_up[ind_miss,1],sample(2:6,1)))
  aux_y=sort(sample(miss_point_down[ind_miss,2]:miss_point_up[ind_miss,2],sample(2:6,1)))
  
  aux=rbind(aux,expand.grid(aux_x,aux_y))
  
}

missing_points[[ind]][[j]]=aux

n_1[[ind]][[j]]=missing_points[[ind]][[j]][,1]
n_2[[ind]][[j]]=missing_points[[ind]][[j]][,2]


Real_X[[ind]][[j]][missing_points[[ind]][[j]][1,1],]



############ 166-184

aux=NULL

for (tag in 1:(length(n_1[[ind]][[j]])-1)) {
  
  aux=c(aux,which(x>=n_1[[ind]][[j]][tag] & x<n_1[[ind]][[j]][tag]+1))
}

na_x[[ind]][[j]]=sort(unique(aux))

aux=NULL

for (tag in 1:(length(n_2[[ind]][[j]])-1)) {
  aux=c(aux,which(y>=n_2[[ind]][[j]][tag] & y<n_2[[ind]][[j]][tag]+1))
}

na_y[[ind]][[j]]=sort(unique(aux))


