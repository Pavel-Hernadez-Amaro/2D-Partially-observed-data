response_int_bifd=function(f_X,f_Beta, x_observations, y_observations, sub_response=50){
  
  n_y=m_y=2*sub_response
  
  W_delta_y=array(dim = (n_y+1)*(m_y+1))
  
  h=(x_observations[length(x_observations)]-x_observations[1])/n_y  # if y_b and y_a are functions of x_observations
  HX=(y_observations[length(y_observations)]-y_observations[1])/m_y # this line should go inside the next for loop (the int_i lopp or outer loop)
  
  
  x = seq(x_observations[1],x_observations[length(x_observations)],h)
  y = seq(y_observations[1],y_observations[length(y_observations)],HX)
  
  
  simp_w_y=rep(1,n_y+1)
  even_y=seq(2,n_y+1-1,2)
  odd_y=seq(3,n_y+1-1,2)
  
  simp_w_y[even_y]=2
  simp_w_y[odd_y]=4
  
  Sim_w_x_y=(h/3)*simp_w_y
  
  W_x_y=(HX/3)*Sim_w_x_y
  W_x_even_y=2*W_x_y
  W_x_odd_y=4*W_x_y
  
  for (aux in 1:(m_y+1)) {
    
    # print(c("aux = ", aux))
    
    if (aux==1 || aux==(m_y+1)) {
      
      W_delta_y[((n_y+1)*(aux-1)+1):((n_y+1)*aux)]= W_x_y
      
    }else{
      
      if (aux%%2==0) {
        
        W_delta_y[((n_y+1)*(aux-1)+1):((n_y+1)*aux)]= W_x_even_y
        
      }else{
        
        W_delta_y[((n_y+1)*(aux-1)+1):((n_y+1)*aux)]= W_x_odd_y
      }}
    
  }
  
  X_eval=eval.bifd(x, y, f_X)
  
  Beta_eval=f_Beta(x, y)
  
  y_int=as.double(t(vec(X_eval)) %*% diag(W_delta_y) %*% vec(Beta_eval))
  
  y_int
  
  
  
}

response_int_H=function(f_X,f_Beta, x_observations, y_observations, sub_response=50){
  
  n_y=m_y=2*sub_response
  
  W_delta_y=array(dim = (n_y+1)*(m_y+1))
  
  h=(x_observations[length(x_observations)]-x_observations[1])/n_y  # if y_b and y_a are functions of x_observations
  HX=(y_observations[length(y_observations)]-y_observations[1])/m_y # this line should go inside the next for loop (the int_i lopp or outer loop)
  
  
  x = seq(x_observations[1],x_observations[length(x_observations)],h)
  y = seq(y_observations[1],y_observations[length(y_observations)],HX)
  
  
  simp_w_y=rep(1,n_y+1)
  even_y=seq(2,n_y+1-1,2)
  odd_y=seq(3,n_y+1-1,2)
  
  simp_w_y[even_y]=2
  simp_w_y[odd_y]=4
  
  Sim_w_x_y=(h/3)*simp_w_y
  
  W_x_y=(HX/3)*Sim_w_x_y
  W_x_even_y=2*W_x_y
  W_x_odd_y=4*W_x_y
  
  for (aux in 1:(m_y+1)) {
    
    # print(c("aux = ", aux))
    
    if (aux==1 || aux==(m_y+1)) {
      
      W_delta_y[((n_y+1)*(aux-1)+1):((n_y+1)*aux)]= W_x_y
      
    }else{
      
      if (aux%%2==0) {
        
        W_delta_y[((n_y+1)*(aux-1)+1):((n_y+1)*aux)]= W_x_even_y
        
      }else{
        
        W_delta_y[((n_y+1)*(aux-1)+1):((n_y+1)*aux)]= W_x_odd_y
      }}
    
  }
  
  X_eval=f_X(x, y)
  
  Beta_eval=f_Beta(x, y)
  
  y_int=as.double(t(vec(X_eval)) %*% diag(W_delta_y) %*% vec(Beta_eval))
  
  y_int
  
  
  
}
