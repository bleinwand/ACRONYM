 try2 <- function(code, silent = FALSE) {
tryCatch(code, error = function(c) {
    if (!silent) {NA}
    else{code}})}




estimate_local_pars_sigma_zero = function(local_matrix, func_form, maxiter = 100, alpha=NA, beta=NA, sigma=NA, rho=NA, Psi_u = NA, Psi_v = NA, currentiter = 0){
  
  num_nodes = nrow(local_matrix)
  num_nodes2 = ncol(local_matrix)
  #This function updates the value (ranging from 0 to 1) of the H function with incorporated noise. 
  #To get the estimated probability, multiply the output of this function by alpha and add beta
  #this function relies on Psi_u, Psi_v, and rho (all embedded in q2), as well as sigma, and the likelihood maximizing choice of epsilons, assuming those parameters are true. Having this available allows for easier optimizing of other parameters
  update_estimates3 = function(){
    entries<<- pnorm((q2 + (sigma_current_opt*eps_estimate))/sqrt(1+sigma_current_opt^2))
  }
  
  #The following functions allows us to update the whole-network parameters: alpha, beta, sigma and rho 
  #The idea is pretty basic: minimize the likelihood. 
  
  
  ab_update_opt = function(alpha_beta){
    if(func_form == "normal"){
      q3 = (qnorm(u2[,2]) + qnorm(u2[,3], sd = alpha_beta[3]))/sqrt(1+alpha_beta[3]^2)
    }
    if(func_form=="exponential"){
      if(alpha_beta[3]==1){
        q3 = qnorm(pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1))
      }else{
        q3 = qnorm(1 - (((alpha_beta[3]*exp(-(qexp(u2[,2]) + qexp(u2[,3], rate=alpha_beta[3])))) - exp(-alpha_beta[3]*((qexp(u2[,2]) + qexp(u2[,3], rate=alpha_beta[3])))))/(alpha_beta[3]-1)))
      }
    }
    
    if(func_form=="negexponential"){
      if(alpha_beta[3]==1){
        q3 = -qnorm(pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1))
      }else{
        q3 = -qnorm(1 - (((alpha_beta[3]*exp(-(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=alpha_beta[3])))) - exp(-alpha_beta[3]*((qexp(1-u2[,2]) + qexp(1-u2[,3], rate=alpha_beta[3])))))/(alpha_beta[3]-1)))
      }
    }
    
    if(func_form == "linear"){
      q3 = qnorm(Linear_H( x= u2[,2] + (alpha_beta[3]*u2[,3]), w= alpha_beta[3]))
      
    }
    
    if(func_form == "uniform"){
      #u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), times=num_nodes2)
      #u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes2),3],function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), each=num_nodes) 
      q3 = qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
    }
    if(func_form == "neguniform"){
      #u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), times=num_nodes2)
      #u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes2),3],function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), each=num_nodes) 
      q3 = -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
    }
    
    
    entries2 = pnorm((q3 + (alpha_beta[4] * eps_estimate))/sqrt(1+ alpha_beta[4]^2))
    opt_ests2 = (alpha_beta[1]*entries2)+alpha_beta[2]
    return((-sum((u2[,1]* log(as.vector(opt_ests2)))[!(is.na(u2[,1]))]) -sum(((1-u2[,1])*log(as.vector(1-opt_ests2)))[!(is.na(u2[,1]))]))
           *ifelse(alpha_beta[4]<.01 || alpha_beta[3]< .01 || alpha_beta[1]>.99 || alpha_beta[1]<.01 || alpha_beta[2]>.99 || alpha_beta[2]<.01 || (alpha_beta[1]+alpha_beta[2])>.999|| (alpha_beta[1]+alpha_beta[2])<.001   , 10000, 1))
  }
  
  
  
  
  ab_update_opt_sym = function(alpha_beta){
    if(func_form == "normal"){
      q3 = (qnorm(u2[,2]) + qnorm(u2[,3], sd = 1))/sqrt(1+1^2)
    }
    if(func_form=="exponential"){
      q3 = qnorm(pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1))
      
    }
    
    if(func_form=="negexponential"){
      q3 = -qnorm(pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1))
      
    }
    
    if(func_form == "linear"){
      q3 = qnorm(Linear_H(x= u2[,2] + u2[,3], w= 1))
    }
    
    if(func_form == "uniform"){
      #u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), times=num_nodes2)
      #u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), each=num_nodes) 
      q3 = qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
    }
    if(func_form == "neguniform"){
      #u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), times=num_nodes2)
      #u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), each=num_nodes) 
      q3 = -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
    }
    
    
    entries2 = pnorm((q3 + (alpha_beta[3] * eps_estimate))/sqrt(1+ alpha_beta[3]^2))
    opt_ests2 = (alpha_beta[1]*entries2)+alpha_beta[2]
    return((-sum((u2[,1]* log(as.vector(opt_ests2)))[!(is.na(u2[,1]))]) -sum(((1-u2[,1])*log(as.vector(1-opt_ests2)))[!(is.na(u2[,1]))]))
           *ifelse(alpha_beta[3]<0 || alpha_beta[1]>1 || alpha_beta[1]<0 || alpha_beta[2]>1 || alpha_beta[2]<0 || (alpha_beta[1]+alpha_beta[2])>.999|| (alpha_beta[1]+alpha_beta[2])<.001 , 10000, 1))
  }
  
  
  
  
  
  
  #we need a separate function to independently update sigma alone at the end after we take our estimated values of epsilon instead of random values.
  
  sigma_analytical_func3 = function(sigma){ 
    #opt_ests2 = (sqrt(1+sigma_^2)*(opt_ests - beta_current_opt)/sqrt(1+sigma_current_opt^2)) + beta_current_opt  
    if(func_form == "normal"){
      opt_ests2 = (alpha_current_opt* pnorm(  ( ((qnorm(u2[,2]) + (rho_current_opt*qnorm(u2[,3])))/sqrt(1+rho_current_opt^2))+ (sigma*eps_estimate))  / sqrt(1+(sigma^2)))) + beta_current_opt
    }
    if(func_form == "exponential"){
      if(rho_current_opt==1){
        opt_ests2 = beta_current_opt + alpha_current_opt*pnorm((qnorm(  pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1)) + (sigma*eps_estimate))/sqrt(1+sigma^2))
      }else{
        opt_ests2 = (alpha_current_opt* pnorm( (qnorm(1 - (((rho_current_opt*exp(-(qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1))) + (sigma*eps_estimate)) / sqrt(1+(sigma^2)))) + beta_current_opt 
      }
    }
    
    if(func_form == "negexponential"){
      if(rho_current_opt==1){
        opt_ests2 = beta_current_opt + alpha_current_opt*pnorm((-qnorm(  pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1)) + (sigma*eps_estimate))/sqrt(1+sigma^2))
      }else{
        opt_ests2 = (alpha_current_opt* pnorm( (-qnorm(1 - (((rho_current_opt*exp(-(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(1-u2[,2]) + qexp(1-u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1))) + (sigma*eps_estimate)) / sqrt(1+(sigma^2)))) + beta_current_opt 
      }
    }
    
    if(func_form == "linear"){
      opt_ests2 = beta_current_opt + alpha_current_opt*pnorm((qnorm(  Linear_H( x= u2[,2] + ( rho_current_opt *u2[,3]),w= rho_current_opt)) + (sigma*eps_estimate))/sqrt(1+sigma^2))
      
    }
    
    
    if(func_form == "uniform"){
      #u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), times=num_nodes2)
      #u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), each=num_nodes) 
      q3 = qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      opt_ests2 = (alpha_current_opt*pnorm((q3+(sigma*eps_estimate))/sqrt(1+sigma^2))) + beta_current_opt
    }
    if(func_form == "neguniform"){
      #u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), times=num_nodes2)
      #u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), each=num_nodes) 
      q3 = -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      opt_ests2 = (alpha_current_opt*pnorm((q3+(sigma*eps_estimate))/sqrt(1+sigma^2))) + beta_current_opt
      
    }
    
    return(-sum((u2[,1]* log(as.vector(opt_ests2)))[!(is.na(u2[,1]))]) -sum(((1-u2[,1])*log(as.vector(1-opt_ests2)))[!(is.na(u2[,1]))]))
  }
  
  ###the following update Psi_u and Psi_v estimates one at a time, so the node needs to be specified in the function. these will actually be estimated with an lapply function later
  #The idea is the same though, just trying to minimize the likelihood by choosing a sociability parameter
  #this names of these functions are different because they were built off an existing function, so normal differentiates from e.g. uniform. These may be subject to change. 
  
  Normal_version_func_Psi_u = function(theta, row_node){
    relevant_rows  = seq(from=row_node, to = length(local_matrix), by=length(Psi_v_current_opt) ) 
    
    if(func_form == "normal"){
      values = pnorm((((theta + (rho_current_opt*qnorm(Psi_v_current_opt))) /(sqrt((rho_current_opt^2) +1))) +(sigma_current_opt*eps_estimate[relevant_rows]))/sqrt(1+sigma_current_opt^2)) 
    }
    
    if(func_form == "exponential"){
      if(rho_current_opt==1){
        values = pnorm((qnorm(  pgamma(qexp(pnorm(theta)) + qexp(Psi_v_current_opt, rate=1), shape=2, rate=1)) + (sigma_current_opt*eps_estimate[relevant_rows]))/sqrt(1+sigma_current_opt^2))
      }else{
        values = pnorm( (qnorm(1 - (((rho_current_opt*exp(-(qexp(pnorm(theta)) + qexp(Psi_v_current_opt, rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(pnorm(theta)) + qexp(Psi_v_current_opt, rate=rho_current_opt)))))/(rho_current_opt-1))) + (sigma_current_opt*eps_estimate[relevant_rows])) / sqrt(1+(sigma_current_opt^2))) 
      }
    }
    
    if(func_form == "negexponential"){
      if(rho_current_opt==1){
        values = pnorm((-qnorm(  pgamma(qexp(1-pnorm(theta)) + qexp(1-Psi_v_current_opt, rate=1), shape=2, rate=1)) + (sigma_current_opt*eps_estimate[relevant_rows]))/sqrt(1+sigma_current_opt^2))
      }else{
        values = pnorm( (-qnorm(1 - (((rho_current_opt*exp(-(qexp(1-pnorm(theta)) + qexp(1-Psi_v_current_opt, rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(1-pnorm(theta)) + qexp(1-Psi_v_current_opt, rate=rho_current_opt)))))/(rho_current_opt-1))) + (sigma_current_opt*eps_estimate[relevant_rows])) / sqrt(1+(sigma_current_opt^2))) 
      }
    }
    
    if(func_form == "linear"){
      values = pnorm((qnorm(  Linear_H( x =pnorm(theta) + (rho_current_opt*Psi_v_current_opt), w=rho_current_opt)) + (sigma_current_opt*eps_estimate[relevant_rows]))/sqrt(1+sigma_current_opt^2))
      
    }
    
    
    if(func_form == "uniform"){
      u1_u = grid_points[which.min(abs(pnorm(theta)- integrated_value))]
      #u2_u =  unlist(lapply(Psi_v_current_opt, function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F) 
      values = pnorm((qnorm(qunif(p=exp(u1_u+u2_u[seq(1, length(u2_u), by = num_nodes)]), min=0, max = 1)) + (sigma_current_opt*eps_estimate[relevant_rows]))/ sqrt(1+(sigma_current_opt^2))) 
    }
    if(func_form == "neguniform"){
      u1_u = grid_points[which.min(abs((1-pnorm(theta))- integrated_value))]
      #u2_u =  unlist(lapply(Psi_v_current_opt, function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F)       
      values = pnorm((-qnorm(qunif(p=exp(u1_u+u2_u[seq(1, length(u2_u), by = num_nodes)]), min=0, max = 1)) + (sigma_current_opt*eps_estimate[relevant_rows]))/ sqrt(1+(sigma_current_opt^2)))
      
    }
    
    
    truths = local_matrix[row_node,]
    opt_ests2 = alpha_current_opt*values +  beta_current_opt
    return(-sum((truths*log(as.vector(opt_ests2)))[!(is.na(truths))]) -sum(((1-truths)*log(as.vector(1-opt_ests2)))[!(is.na(truths))]))
    
  }
  
  
  Normal_version_func_Psi_v = function(theta, col_node){
    relevant_cols  = (((col_node-1)*length(Psi_u_current_opt))+1):((col_node)*length(Psi_u_current_opt)) 
    if(func_form == "normal"){
      values = pnorm((((qnorm(Psi_u_current_opt) + (rho_current_opt*theta)) /(sqrt((rho_current_opt^2) +1))) +(sigma_current_opt*eps_estimate[relevant_cols]))/sqrt(1+sigma_current_opt^2)) 
    }
    
    if(func_form == "exponential"){
      if(rho_current_opt==1){
        values =  pnorm((qnorm(  pgamma(qexp(Psi_u_current_opt) + qexp(pnorm(theta), rate=1), shape=2, rate=1)) + (sigma_current_opt*eps_estimate[relevant_cols]))/sqrt(1+sigma_current_opt^2))
      }else{
        values = pnorm( (qnorm(1 - (((rho_current_opt*exp(-(qexp(Psi_u_current_opt) + qexp(pnorm(theta), rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(Psi_u_current_opt) + qexp(pnorm(theta), rate=rho_current_opt)))))/(rho_current_opt-1))) + (sigma_current_opt*eps_estimate[relevant_cols])) / sqrt(1+(sigma_current_opt^2))) 
      }
    }
    if(func_form == "negexponential"){
      if(rho_current_opt==1){
        values =  pnorm((-qnorm(  pgamma(qexp(1-Psi_u_current_opt) + qexp(1-pnorm(theta), rate=1), shape=2, rate=1)) + (sigma_current_opt*eps_estimate[relevant_cols]))/sqrt(1+sigma_current_opt^2))
      }else{
        values = pnorm( (-qnorm(1 - (((rho_current_opt*exp(-(qexp(1-Psi_u_current_opt) + qexp(1-pnorm(theta), rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(1-Psi_u_current_opt) + qexp(1-pnorm(theta), rate=rho_current_opt)))))/(rho_current_opt-1))) + (sigma_current_opt*eps_estimate[relevant_cols])) / sqrt(1+(sigma_current_opt^2))) 
      }
    }
    
    
    if(func_form == "linear"){
      values = pnorm((qnorm(  Linear_H( x=Psi_u_current_opt + (rho_current_opt*pnorm(theta)), w=rho_current_opt)) + (sigma_current_opt*eps_estimate[relevant_cols]))/sqrt(1+sigma_current_opt^2))
      
    }
    
    if(func_form == "uniform"){
      #u1_u =  unlist(lapply(Psi_u_current_opt, function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F) 
      u2_u = grid_points[which.min(abs(pnorm(theta)- integrated_value))]
      values = pnorm((qnorm(qunif(p=exp(u1_u[1:num_nodes2]+u2_u), min=0, max = 1)) + (sigma_current_opt*eps_estimate[relevant_cols]))/ sqrt(1+(sigma_current_opt^2))) 
    }
    if(func_form == "neguniform"){
      #u1_u =  unlist(lapply(Psi_u_current_opt, function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F)   
      u2_u = grid_points[which.min(abs((1-pnorm(theta))- integrated_value))]
      values = pnorm((-qnorm(qunif(p=exp(u1_u[1:num_nodes2]+u2_u), min=0, max = 1)) + (sigma_current_opt*eps_estimate[relevant_cols]))/ sqrt(1+(sigma_current_opt^2)))
      
    }
    
    
    truths = local_matrix[,col_node]
    opt_ests2 = alpha_current_opt*values +  beta_current_opt
    return(-sum((truths*log(as.vector(opt_ests2)))[!(is.na(truths))]) -sum(((1-truths)*log(as.vector(1-opt_ests2)))[!(is.na(truths))]))
  }
  
  
  # below are the 2 functions that allow for estimating all the Psi_u values with one line of code
  Psi_u_opt_lapply =  function(x){
    return(pnorm(optim(par = 0, fn= Normal_version_func_Psi_u, row_node=x, method = "Brent", upper=4, lower=-4)$par))
  }
  
  
  
  Psi_v_opt_lapply = function(x){
    return(pnorm(optim(par = qnorm(0), fn= Normal_version_func_Psi_v, col_node=x, method = "Brent", upper=4, lower=-4)$par))
  }
  
  
  
  
  
  
  
  
  
  
  #the following function sets the parameter estimates to a baseline before estimation
  
  reset_pars = function(){
    #Nothing clever here, just ranking each row node and column node by their degree
    Psi_u_current_opt <<- rank(rowSums(local_matrix, na.rm = T))/(num_nodes+1)
    Psi_v_current_opt <<- rank(colSums(local_matrix, na.rm = T))/(num_nodes2+1)
    #Psi_v_current_opt <<- (colMeans(local_matrix)-min(colMeans(local_matrix)))/diff(range(colMeans(local_matrix)))
    #Psi_v_current_opt <<- ((Psi_v_current_opt -.5)*num_nodes2/(num_nodes2+1))+.5
    #Psi_u_current_opt <<- (rowMeans(local_matrix)-min(rowMeans(local_matrix)))/diff(range(rowMeans(local_matrix)))
    #Psi_u_current_opt <<- ((Psi_u_current_opt -.5)*num_nodes/(num_nodes+1))+.5
    
    
    
    #Choose beta to be the minimum percentage of existing edges for any particular node
    beta_current_opt <<- min(c(rowSums(local_matrix, na.rm = T)/nrow(local_matrix), colSums(local_matrix, na.rm = T)/ncol(local_matrix)))
    #alpha is based on the expected density  density = alpha/2 +beta 
    alpha_current_opt <<- max(.01, min(.99,  2*(mean(local_matrix, na.rm = T) - beta_current_opt)))
    
    #Rho and sigma start at 1, "at balance" in a sense
    sigma_current_opt <<- 0 #1 CHANGESIGMA
    rho_current_opt <<- 1
    #valence<<-1
    
    #We're just initializing the likelihood maximizing values of epsilons to be all 0
    eps_estimate <<- rep(0, length(local_matrix))
    
    #u2 is a matrix containing vectorized edge presence, estimated row sociabilities, and estimated column sociabilities, allowing for faster computation
    u2 <<- cbind(as.vector(local_matrix), rep(Psi_u_current_opt, times= length(Psi_v_current_opt)), rep(Psi_v_current_opt, each=length(Psi_u_current_opt)))
    
    #q2 is the output of the H-function, with nothing else included
    if(func_form == "normal"){
      q2 <<- (qnorm(u2[,2]) + qnorm(u2[,3], sd = rho_current_opt))/sqrt(1+rho_current_opt^2)
    }
    if(func_form == "exponential"){
      #   valence <<- 1.05 - ifelse(sign(diff(range(rowSums(local_matrix)))/diff(range(colSums(local_matrix))))< 1, .1, 0)
      # rho_current_opt<<- valence
      
      q2 <<-  qnorm(pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1)) 
      
      # q2 <<- qnorm(1 - (((rho_current_opt*exp(-(qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1)))
      
    }
    if(func_form == "negexponential"){
      
      q2 <<-  -qnorm(pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1)) 
      
    }
    
    if(func_form == "linear"){
      q2 <<- qnorm(Linear_H( x=u2[,2] + u2[,3], w= 1))
      
    }
    
    if(func_form == "uniform"){
      u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), times=num_nodes2)
      u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), each=num_nodes) 
      q2 <<- qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      
    }
    if(func_form == "neguniform"){
      u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), times=num_nodes2)
      u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), each=num_nodes) 
      q2 <<- -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      
    }
    
    
  }
  
  #the following function is used to take our parameter estimates and convert them into edge estimates (without the alpha and beta part, but that will be included in a single line of code) . 
  
  #This function gets the estimates outside of alpha and beta (that is, assuming, alpha=1 and beta=0),which is the 3-d H-function incorporating both nodes and epsilon.
  getting_ests_opt =  function(group1_Psi,group2_Psi){
    if(func_form=="normal"){
      pnorm(qnorm(pnorm(qnorm(group1_Psi) + qnorm(group2_Psi, sd = rho_current_opt), sd = sqrt(1+rho_current_opt^2))), sd = sqrt(1+2*(sigma_current_opt^2))) 
    }else{
      
      if(func_form=="exponential"){
        if(rho_current_opt==1){
          pnorm((qnorm(  pgamma(qexp(group1_Psi) + qexp(group2_Psi, rate=1), shape=2, rate=1)) )/sqrt(1+2*sigma_current_opt^2))
        }else{
          pnorm(qnorm(1 - (((rho_current_opt*exp(-(qexp(group1_Psi) + qexp(group2_Psi, rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(group1_Psi) + qexp(group2_Psi, rate=rho_current_opt)))))/(rho_current_opt-1))), sd = sqrt(1+2*(sigma_current_opt^2))) 
        }
      }else{
        
        if(func_form=="negexponential"){
          if(rho_current_opt==1){
            pnorm((-qnorm(  pgamma(qexp(1-group1_Psi) + qexp(1-group2_Psi, rate=1), shape=2, rate=1)) )/sqrt(1+2*sigma_current_opt^2))
          }else{
            pnorm(-qnorm(1 - (((rho_current_opt*exp(-(qexp(1-group1_Psi) + qexp(1-group2_Psi, rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(1-group1_Psi) + qexp(1-group2_Psi, rate=rho_current_opt)))))/(rho_current_opt-1))), sd = sqrt(1+2*(sigma_current_opt^2))) 
          }
        }else{
          if(func_form == "linear"){
            pnorm((qnorm( Linear_H( x= u2[,2] + (rho_current_opt*u2[,3]), w= rho_current_opt)) )/sqrt(1+2*sigma_current_opt^2))
            
          }else{
            
            if(func_form == "uniform" ){
              #u1_u = unlist(lapply(group1_Psi, function(x){grid_points[which.min(abs(x- integrated_value))]}), use.names = F)
              #u2_u = unlist(lapply(group2_Psi, function(x){grid_points[which.min(abs(x- integrated_value))]}), use.names = F)
              pnorm(qnorm(exp(group1_Psi+group2_Psi))/sqrt(1+2*(sigma_current_opt^2)))
            }else{
              
              if(func_form == "neguniform"){
                #u1_u = unlist(lapply(group1_Psi, function(x){grid_points[which.min(abs((1-x)- integrated_value))]}), use.names = F)
                #u2_u = unlist(lapply(group2_Psi, function(x){grid_points[which.min(abs((1-x)- integrated_value))]}), use.names = F)
                pnorm(qnorm(exp(group1_Psi+group2_Psi))/sqrt(1+2*(sigma_current_opt^2)))
                
              }
              
            }
          }
        }
      }
    }
  }
  
  
  #the following function is only used to accept or reject changes to all Psi_u or Psi_v values in one fell swoop.  
  Psi_update_func = function(Psis){
    f =  rep(Psis[1:num_nodes], times= num_nodes2)
    g =  rep(Psis[(num_nodes+1):(num_nodes+num_nodes2)], each= num_nodes)
    if(func_form == "normal"){
      q3 = (qnorm(f) + qnorm(g, sd = rho_current_opt))/sqrt(1+rho_current_opt^2)
    }
    if(func_form == "exponential"){
      if(rho_current_opt==1){
        q3 = qnorm(  pgamma(qexp(f) + qexp(g, rate=1), shape=2, rate=1))
      }else{
        q3 = qnorm(1 - (((rho_current_opt*exp(-(qexp(f) + qexp(g, rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(f) + qexp(g, rate=rho_current_opt)))))/(rho_current_opt-1)))
        
      }
    }
    
    
    if(func_form == "negexponential"){
      if(rho_current_opt==1){
        q3 = -qnorm(  pgamma(qexp(1-f) + qexp(1-g, rate=1), shape=2, rate=1))
      }else{
        q3 = -qnorm(1 - (((rho_current_opt*exp(-(qexp(1-f) + qexp(1-g, rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(1-f) + qexp(1-g, rate=rho_current_opt)))))/(rho_current_opt-1)))
        
      }
    }
    
    if(func_form == "linear"){
      
      q3 = qnorm(Linear_H( x=f + (rho_current_opt*g), w= rho_current_opt))
      
    }
    
    if(func_form == "uniform"){
      u1_u = rep(unlist(lapply(f[1:num_nodes], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), times=num_nodes2)
      u2_u =  rep(unlist(lapply(g[seq(1, length(g), by=num_nodes)],function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), each=num_nodes) 
      q3 = qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
    }
    if(func_form == "neguniform"){
      u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), times=num_nodes2)
      u2_u =  rep(unlist(lapply(u2[seq(1, length(g), by=num_nodes),3],function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), each=num_nodes) 
      q3 = -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      
    }
    
    
    entries2 = pnorm((q3 + (sigma_current_opt * eps_estimate))/sqrt(1+ sigma_current_opt^2))
    opt_ests2 = (alpha_current_opt * entries2)+ beta_current_opt
    return( -sum(u2[,1]* log(as.vector(opt_ests2)), na.rm = T) -sum((1-u2[,1])*log(as.vector(1-opt_ests2)), na.rm = T))
  }
  
  
  
  epsFinder = function(epsilon, q, alpha, beta, sigma, u){
    
    sigRatio = sigma/sqrt(1+sigma^2)
    
    z = (q/sqrt(1+sigma^2)) + (sigRatio*epsilon)
    
    `if`(is.na(u), 0, `if`(u==1,
         (
           epsilon - 
             (
               (sigRatio * alpha * dnorm(z)) /
                 (alpha*pnorm(z)+beta)   
             )
         )^2 ,
         
         sq_error =(
           epsilon + 
             (
               (sigRatio * alpha * dnorm(z)) /
                 (1-(alpha*pnorm(z)+beta))   
             )
         )^2 
    )
    
    )
    
  }
  
  
  
  
  epsFinder_lapply =  function(x){
    return(optim(par = eps_estimate[x], fn= epsFinder, q = q2[x], alpha = alpha_current_opt, beta=beta_current_opt, sigma=sigma_current_opt, u = u2[x,1], method = "Brent", upper=4, lower=-4, control = list(abstol=.000025))$par)
  }
  
  
  
  update_eps = function(){
    #eps_estimate <<-   rnorm(length(eps_estimate), mean=
    # (u2[,1]*(alpha_current_opt*dnorm(qnorm(entries))*abs(sigma_current_opt/sqrt(1+sigma_current_opt^2)))/(alpha_current_opt*entries+beta_current_opt)) - ((1-u2[,1])*(alpha_current_opt*dnorm(qnorm(entries)))*abs(sigma_current_opt/sqrt(1+sigma_current_opt^2))/(1-(alpha_current_opt*entries+beta_current_opt)))
    #)
    
    eps_estimate <<-rnorm(length(eps_estimate), mean= unlist(lapply(1:length(eps_estimate), FUN = epsFinder_lapply)))
  }
  
  
  
  update_a_b_rho_sig = function(){
    #opt4 =  nlm(f= ab_update_opt, p = c(alpha_current_opt,beta_current_opt, 1, 1))
    opt4 = nlminb(start = c(alpha_current_opt,beta_current_opt, rho_current_opt, 1), objective = ab_update_opt, lower= c(.01, .01, .1, 0), upper = c(.99, .99, 10, 0 ) )#CHANGESIGMA
    # if(opt4$minimum <= current_likelihood){
    # 
    #   alpha_current_opt <<- opt4$estimate[1]
    #   beta_current_opt <<- opt4$estimate[2]
    #   rho_current_opt <<- opt4$estimate[3]
    #   sigma_current_opt <<- opt4$estimate[4]
    if(opt4$objective <= current_likelihood && is.numeric(opt4$par[3])){
      alpha_current_opt <<- opt4$par[1]
      beta_current_opt <<- opt4$par[2]
      rho_current_opt <<- opt4$par[3]
      sigma_current_opt <<- 0#CHANGESIGMA
      
      
      if(func_form == "normal"){
        q2 <<- (qnorm(u2[,2]) + qnorm(u2[,3], sd = rho_current_opt))/sqrt(1+rho_current_opt^2)
      }
      if(func_form == "exponential"){
        if(rho_current_opt==1){
          q2<<- qnorm(  pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1))
        }else{
          q2 <<- qnorm(1 - (((rho_current_opt*exp(-(qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1)))
        }
      }
      
      if(func_form == "negexponential"){
        if(rho_current_opt==1){
          q2 <<- -qnorm(  pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1))
        }else{
          q2 <<- -qnorm(1 - (((rho_current_opt*exp(-(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(1-u2[,2]) + qexp(1-u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1)))
        }
      }
      
      
      if(func_form == "linear"){
        
        q2 <<- qnorm(Linear_H( x=u2[,2] + (rho_current_opt*u2[,3]), w=rho_current_opt))
        
      }
      
      
      #    if(func_form == "uniform"){
      #    
      #     q2 <<- qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      # }
      # if(func_form == "neguniform"){
      #    
      #     q2 <<- -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      # 
      # }
      
      
      #update estimates
      update_estimates3()                           
      #current_likelihood <<- opt4$minimum
      current_likelihood = opt4$objective
    }
    
  }
  
  
  
  
  update_a_b_rho_sig_sym = function(){
    #opt4 =  nlm(f= ab_update_opt_sym, p = c(alpha_current_opt,beta_current_opt, 1))
    opt4 = nlminb(start = c(alpha_current_opt,beta_current_opt, 1), objective = ab_update_opt_sym, lower= c(.01, .01,  0), upper = c(.99, .99, 0 ) ) #CHANGESIGMA
    # if(opt4$minimum <= current_likelihood){
    # 
    #   alpha_current_opt <<- opt4$estimate[1]
    #   beta_current_opt <<- opt4$estimate[2]
    #   sigma_current_opt <<- opt4$estimate[3]
    #   
    if(opt4$objective <= current_likelihood){
      
      alpha_current_opt <<- opt4$par[1]
      beta_current_opt <<- opt4$par[2]
      sigma_current_opt <<- 0 #CHANGESIGMA opt4$par[3]
      
      # if(func_form == "normal"){
      # q2 <<- (qnorm(u2[,2]) + qnorm(u2[,3], sd = 1))/sqrt(1+1^2)
      #       }
      # if(func_form == "exponential"){
      #  qnorm(  pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1))
      # 
      # }
      #   
      #   if(func_form == "negexponential"){
      #  -qnorm(  pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1))
      # }
      # 
      #   
      #   if(func_form == "uniform"){
      #     #u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), times=num_nodes2)
      #     #u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), each=num_nodes) 
      #     qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      # }
      # if(func_form == "neguniform"){
      #    #u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), times=num_nodes2)
      #     #u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), each=num_nodes) 
      #      -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      #     
      # 
      # }
      
      
      
      #update estimates
      update_estimates3()                           
      #current_likelihood <<- opt4$minimum
      current_likelihood = opt4$objective
    }
    
  }
  
  
  
  
  
  
  
  
  update_Psi_u = function(){
    Psi_u_temp = unlist(lapply(1:num_nodes, FUN = Psi_u_opt_lapply), use.names = F)
    if(Psi_update_func(c(Psi_u_temp, Psi_v_current_opt))<= current_likelihood){
      Psi_u_current_opt <<- Psi_u_temp
      u2[,2] <<-  rep(Psi_u_current_opt, times= length(Psi_v_current_opt))
      if(func_form == "normal"){
        q2 <<- (qnorm(u2[,2]) + qnorm(u2[,3], sd = rho_current_opt))/sqrt(1+rho_current_opt^2)
      }
      if(func_form == "exponential"){
        if(rho_current_opt==1){
          q2<<- qnorm(  pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1)) 
        }else{
          q2 <<- qnorm(1 - (((rho_current_opt*exp(-(qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1)))
        }
      }
      
      if(func_form == "negexponential"){
        if(rho_current_opt==1){
          q2<<- -qnorm(  pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1)) 
        }else{
          q2 <<- -qnorm(1 - (((rho_current_opt*exp(-(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(1-u2[,2]) + qexp(1-u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1)))
        }
      }
      
      
      if(func_form == "linear"){
        
        q2 <<- qnorm(Linear_H( x=u2[,2] + (rho_current_opt*u2[,3]), w= rho_current_opt))
        
      }
      
      if(func_form == "uniform"){
        u1_u <<- rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), times=num_nodes2)
        u2_u <<-  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), each=num_nodes) 
        q2<<- qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
        
      }
      if(func_form == "neguniform"){
        u1_u <<- rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), times=num_nodes2)
        u2_u <<-  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), each=num_nodes) 
        q2<<- -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
        
      }
      
      
      update_estimates3()
      current_likelihood <<- Psi_update_func(c(Psi_u_current_opt, Psi_v_current_opt))
    }
  }
  
  
  
  update_Psi_v = function(){
    Psi_v_temp = unlist(lapply(1:num_nodes2, FUN = Psi_v_opt_lapply), use.names = F)
    if(Psi_update_func(c(Psi_u_current_opt, Psi_v_temp))<= current_likelihood){
      Psi_v_current_opt <<- Psi_v_temp
      u2[,3] <<-  rep(Psi_v_current_opt, each= length(Psi_u_current_opt))
      if(func_form == "normal"){
        q2 <<- (qnorm(u2[,2]) + qnorm(u2[,3], sd = rho_current_opt))/sqrt(1+rho_current_opt^2)
      }
      if(func_form == "exponential"){
        if(rho_current_opt==1){
          q2<<- qnorm(  pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1)) 
        }else{
          q2 <<- qnorm(1 - (((rho_current_opt*exp(-(qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1)))
        }
      }
      
      if(func_form == "negexponential"){
        if(rho_current_opt==1){
          q2<<- -qnorm(  pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1)) 
        }else{
          q2 <<- -qnorm(1 - (((rho_current_opt*exp(-(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(1-u2[,2]) + qexp(1-u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1)))
        }
      }
      
      if(func_form == "linear"){
        
        q2 <<- qnorm(Linear_H(x=u2[,2] + (rho_current_opt*u2[,3]), w=rho_current_opt))
        
      }
      
      if(func_form == "uniform"){
        u1_u <<- rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), times=num_nodes2)
        u2_u <<-  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), each=num_nodes) 
        q2<<- qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
        
      }
      if(func_form == "neguniform"){
        u1_u <<- rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), times=num_nodes2)
        u2_u <<-  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), each=num_nodes) 
        q2<<- -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
        
      }
      update_estimates3() 
      current_likelihood <<- Psi_update_func(c(Psi_u_current_opt, Psi_v_current_opt))
    }
  }
  
  
  
  
  
  update_Psi_u_sym = function(){
    Psi_u_temp = unlist(lapply(1:num_nodes, FUN = Psi_u_opt_lapply), use.names = F)
    if(Psi_update_func(c(Psi_u_temp, Psi_u_temp))<= current_likelihood){
      Psi_u_current_opt <<- Psi_u_temp
      Psi_v_current_opt <<- Psi_u_temp
      u2[,2] <<-  rep(Psi_u_current_opt, times= length(Psi_v_current_opt))
      u2[,3] <<-  rep(Psi_v_current_opt, each= length(Psi_u_current_opt))
      
      if(func_form == "normal"){
        q2 <<- (qnorm(u2[,2]) + qnorm(u2[,3], sd = 1))/sqrt(1+1^2)
      }
      if(func_form == "exponential"){
        q2<<- qnorm(  pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1)) 
        
      }
      
      if(func_form == "negexponential"){
        q2<<- -qnorm(  pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1)) 
      }
      
      if(func_form == "linear"){
        
        q2 <<- qnorm(Linear_H( x=u2[,2] +u2[,3], w= 1))
        
      }
      
      if(func_form == "uniform"){
        uu = unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F)
        u1_u <<- rep(uu, times=num_nodes2)
        u2_u <<- rep(uu, each=num_nodes) 
        q2<<- qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
        
      }
      if(func_form == "neguniform"){
        uu = unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F)
        u1_u <<- rep(uu, times=num_nodes2)
        u2_u <<-  rep(uu, each=num_nodes) 
        q2<<- -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
        
      }
      update_estimates3()
      current_likelihood <<- Psi_update_func(c(Psi_u_current_opt, Psi_v_current_opt))
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  #check if the matrix is symmetric. If so, we keep rho as 1, and we only have one set of nodes 
  symmetry_check =  isSymmetric(local_matrix) 
  
  #if we have some existing guesses (usually from an earlier analysis), use these. Otherwise we need more naive estimates
  if(sum(is.na(c(alpha, beta, sigma, rho, Psi_u, Psi_v)))==0){
    alpha_current_opt = alpha  
    beta_current_opt = beta
    sigma_current_opt = sigma
    rho_current_opt = ifelse(symmetry_check, 1, rho)
    Psi_u_current_opt = Psi_u
    Psi_v_current_opt = Psi_v
    
    Psi_u_best = Psi_u_current_opt
    Psi_v_best = Psi_v_current_opt
    alpha_best = alpha_current_opt
    beta_best = beta_current_opt
    rho_best = rho_current_opt
    sigma_best = sigma_current_opt
    
    
    eps_estimate = rep(0, length(local_matrix))
    
    #u2 is a matrix containing vectorized edge presence, estimated row sociabilities, and estimated column sociabilities, allowing for faster computation
    u2 = cbind(as.vector(local_matrix), rep(Psi_u_current_opt, times= length(Psi_v_current_opt)), rep(Psi_v_current_opt, each=length(Psi_u_current_opt)))
    
    #q2 is the output of the H-function, with nothing else included
    if(func_form == "normal"){
      q2 = (qnorm(u2[,2]) + qnorm(u2[,3], sd = rho_current_opt))/sqrt(1+rho_current_opt^2)
    }
    if(func_form == "exponential"){
      #   valence <<- 1.05 - ifelse(sign(diff(range(rowSums(local_matrix)))/diff(range(colSums(local_matrix))))< 1, .1, 0)
      # rho_current_opt<<- valence
      
      q2 =  qnorm(pgamma(qexp(u2[,2]) + qexp(u2[,3], rate=1), shape=2, rate=1)) 
      
      # q2 <<- qnorm(1 - (((rho_current_opt*exp(-(qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))) - exp(-rho_current_opt*((qexp(u2[,2]) + qexp(u2[,3], rate=rho_current_opt)))))/(rho_current_opt-1)))
      
    }
    if(func_form == "negexponential"){
      
      q2 =  -qnorm(pgamma(qexp(1-u2[,2]) + qexp(1-u2[,3], rate=1), shape=2, rate=1)) 
      
    }
    
    if(func_form == "linear"){
      
      q2 = qnorm(Linear_H( x = u2[,2] + u2[,3], w=1))
      
    }
    
    
    if(func_form == "uniform"){
      u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), times=num_nodes2)
      u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs(x- integrated_value))]}),use.names = F), each=num_nodes) 
      q2= qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      
    }
    if(func_form == "neguniform"){
      u1_u = rep(unlist(lapply(u2[1:num_nodes,2], function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), times=num_nodes2)
      u2_u =  rep(unlist(lapply(u2[seq(1, nrow(u2), by=num_nodes),3],function(x){grid_points[which.min(abs((1-x)- integrated_value))]}),use.names = F), each=num_nodes) 
      q2= -qnorm(qunif(p=exp(u1_u+u2_u), min=0, max = 1)) 
      
    }
  }else{
    reset_pars()
    if(min(rowMeans(local_matrix, na.rm = T))==0||max(rowMeans(local_matrix, na.rm = T))==1){
      local_matrix[c(which(rowMeans(local_matrix, na.rm = T)==0),which(rowMeans(local_matrix, na.rm = T)==1) ), ]=NA
    }
    if(min(colMeans(local_matrix, na.rm = T))==0||max(colMeans(local_matrix, na.rm = T))==1){
      local_matrix[, c(which(colMeans(local_matrix, na.rm = T)==0),which(colMeans(local_matrix, na.rm = T)==1) )]=NA
    }
    Psi_u_best = Psi_u_current_opt
    Psi_v_best = Psi_v_current_opt
    alpha_best = alpha_current_opt
    beta_best = beta_current_opt
    rho_best = rho_current_opt
    sigma_best = sigma_current_opt
  }
  
  
  
  
  update_estimates3()
  if(func_form %in% c("uniform", "neguniform")){
    opt_ests = (alpha_current_opt*outer(u1_u[1:num_nodes2], u2_u[seq(1, length(u2_u), by = num_nodes)], getting_ests_opt) ) + beta_current_opt
  }else{
    opt_ests = (alpha_current_opt*outer(Psi_u_current_opt, Psi_v_current_opt, getting_ests_opt)) + beta_current_opt
  }
  
  
  
  current_likelihood = -sum(u2[,1]*log(as.vector(opt_ests)), na.rm = T) -sum((1-u2[,1])*log(as.vector(1-opt_ests)), na.rm = T)
  best_likelihood = current_likelihood
  #start an iterator and "error", that is difference in likelihood, which we will use to stop the algorithm 
  iterator = currentiter
  
  update_eps()
  
  #error_from_truth tracks the likelihood of our estimated matrix
  error_from_truth=rep(0, maxiter)
  
  
  while(iterator<maxiter){ 
    iterator= iterator+1
    print(paste("On iteration", iterator, "of",  maxiter, "for", ifelse(func_form=="normal", "symmetric", ifelse(func_form=="exponential", "concave", ifelse(func_form=="negexponential", "convex", ifelse(func_form=="linear", "linear", ifelse(func_form=="uniform", "LS1", "LS2"))))))) 
    #First we assume our parameters are right, and get the most likely epsilon values for each edge based on the observed network, and update our estimated probabilities based on thes estimated epsilon values
    #eps_estimate = rnorm(unlist(lapply(1:length(eps_estimate), est_eps_lapply), use.names = F))
    
    
    update_eps()
    
    update_estimates3()
    
    
    #randomize whether we update whole network parameters, Psi_u, or Psi_v first. We need special functions for when the matrix is symmetric, as we constrain rho to be equal, and we require that Psi_u=Psi_v
    opt_order = sample(1:4, 1)
    if(symmetry_check==FALSE){
      if(opt_order==1){
        update_a_b_rho_sig()
        update_Psi_u()
        update_Psi_v()
      }
      if(opt_order==2){
        update_a_b_rho_sig()
        update_Psi_v()
        update_Psi_u()
      }
      
      if(opt_order==3){
        update_Psi_u()
        update_Psi_v()
        update_a_b_rho_sig()
      }
      
      if(opt_order==4){
        update_Psi_v()
        update_Psi_u()
        update_a_b_rho_sig()
      }
    }else{
      if(opt_order%%2==0){
        update_a_b_rho_sig_sym()
        update_Psi_u_sym()
      }
      if(opt_order%%2==1){
        update_Psi_u_sym()  
        update_a_b_rho_sig_sym()
      }
      
      
    }
    #Update Psi_u and Psi_V
    
    
    
    #Update our estimated probabilities, get the likelihood based on these estimates, ensure the likelihood is improving, 
    if(func_form %in% c("uniform", "neguniform")){
      opt_ests = (alpha_current_opt*outer(u1_u[1:num_nodes2], u2_u[seq(1, length(u2_u), by = num_nodes)], getting_ests_opt) ) + beta_current_opt
    }else{
      opt_ests = (alpha_current_opt*outer(Psi_u_current_opt, Psi_v_current_opt, getting_ests_opt)) + beta_current_opt
    }
    
    current_likelihood = -sum(u2[,1]*log(as.vector(opt_ests)), na.rm = T) -sum((1-u2[,1])*log(as.vector(1-opt_ests)), na.rm = T)
    best_likelihood = min(current_likelihood, best_likelihood)
    error_from_truth[iterator] = current_likelihood#mean((opt_ests-no_error_probs_matrix)^2)
    
    
    
    #in case the likelihood doesn't improve, we want to hold onto the best parameters
    if(current_likelihood ==best_likelihood){
      Psi_u_best = Psi_u_current_opt
      Psi_v_best = Psi_v_current_opt
      alpha_best = alpha_current_opt
      beta_best = beta_current_opt
      rho_best = rho_current_opt
      sigma_best = sigma_current_opt
    }
    
  } 
  
  
  
  if(currentiter>0){
    
    Psi_u_current_opt = Psi_u_best
    Psi_v_current_opt = Psi_v_best
    alpha_current_opt = alpha_best
    beta_current_opt = beta_best
    rho_current_opt = rho_best 
    sigma_current_opt = sigma_best
    eps_estimate = rep(0, length(local_matrix))
    
    update_estimates3()
    #eps_estimate = unlist(lapply(1:length(eps_estimate), est_eps_lapply), use.names = F)
    # eps_estimate=   (u2[,1]*(alpha_current_opt*dnorm(qnorm(entries))*abs(sigma_current_opt/sqrt(1+sigma_current_opt^2)))/(alpha_current_opt*entries+beta_current_opt)) - ((1-u2[,1])*(alpha_current_opt*dnorm(qnorm(entries)))*abs(sigma_current_opt/sqrt(1+sigma_current_opt^2))/(1-(alpha_current_opt*entries+beta_current_opt)))
    eps_estimate <<- unlist(lapply(1:length(eps_estimate), FUN = epsFinder_lapply))
    #eps_estimate =  qnorm(order(eps_estimate)/(length(eps_estimate)+1))
    update_estimates3()
    sigma_current_opt = 0 #CHANGESIGMA #optim(par = 1, fn = sigma_analytical_func3, method = "Brent", lower = .1, upper= 10)$par
    
  }
  
  
  #opt_ests is an estimate of the expected value of probs_matrix - it ignores epsilons but divides by sqrt(1+2*sigma^2) 
  if(func_form %in% c("uniform", "neguniform")){
    opt_ests = (alpha_current_opt*outer(u1_u[1:num_nodes2], u2_u[seq(1, length(u2_u), by = num_nodes)], getting_ests_opt) ) + beta_current_opt
  }else{
    opt_ests = (alpha_current_opt*outer(Psi_u_current_opt, Psi_v_current_opt, getting_ests_opt)) + beta_current_opt
  }
  
  #opt_ests_with_epsilons is the estimate of the exact value of probs_matrix, that is, including estimated epsilons
  opt_ests_with_epsilons = alpha_current_opt*(pnorm((q2 + (sigma_current_opt* eps_estimate))/sqrt(1+ sigma_current_opt^2)))+beta_current_opt
  
  #no_error_opt_ests is the estimate of no_error_probs_matrix - it doesn't include sigma or epsilon at all 
  no_error_opt_ests = (alpha_current_opt*pnorm(q2))+beta_current_opt
  
  return(list(alpha_current_opt, beta_current_opt, sigma_current_opt, rho_current_opt, Psi_u_current_opt, Psi_v_current_opt, func_form, best_likelihood, opt_ests, opt_ests_with_epsilons, no_error_opt_ests))
}




full_subnetwork_estimation_sigma_zero = function(subnetwork){  
local_start_time = proc.time()
preiter = 20
 try2 <- function(code, silent = FALSE) {
tryCatch(code, error = function(c) {
    if (!silent) {NA}
    else{code}})}
    
  normal5 = try2(estimate_local_pars_sigma_zero(local_matrix = subnetwork, func_form = "normal", maxiter = preiter))
  exp5 = try2(estimate_local_pars_sigma_zero(local_matrix = subnetwork, func_form = "exponential", maxiter = preiter))
  negexp5 = try2(estimate_local_pars_sigma_zero(local_matrix = subnetwork, func_form = "negexponential", maxiter = preiter))
  linear5 = try2(estimate_local_pars_sigma_zero(local_matrix = subnetwork, func_form = "linear", maxiter = preiter))
  
  
  if(max(is.na(normal5))){normal5 = (list(NA, NA, NA, NA, NA, NA, NA, Inf, NA, NA, NA))}
  if(max(is.na(exp5))){exp5 = (list(NA, NA, NA, NA, NA, NA, NA, Inf, NA, NA, NA))}
  if(max(is.na(negexp5))){negexp5 =(list(NA, NA, NA, NA, NA, NA, NA, Inf, NA, NA, NA))}
  if(max(is.na(linear5))){linear5 =(list(NA, NA, NA, NA, NA, NA, NA, Inf, NA, NA, NA))}
  
  if(normal5[[8]] == min(normal5[[8]], exp5[[8]], negexp5[[8]], linear5[[8]])){
    result = estimate_local_pars_sigma_zero(local_matrix = subnetwork, func_form = "normal", maxiter = 100, alpha=normal5[[1]], beta=normal5[[2]], sigma=normal5[[3]], rho=normal5[[4]], Psi_u = normal5[[5]], Psi_v = normal5[[6]], currentiter = preiter)
  }else{
    if(exp5[[8]] == min(normal5[[8]], exp5[[8]], negexp5[[8]], linear5[[8]])){
      result = estimate_local_pars_sigma_zero(local_matrix = subnetwork, func_form = "exponential", maxiter = 100, alpha=exp5[[1]], beta=exp5[[2]], sigma=exp5[[3]], rho=exp5[[4]], Psi_u = exp5[[5]], Psi_v = exp5[[6]], currentiter = preiter)
    }else{
      if(negexp5[[8]] == min(normal5[[8]], exp5[[8]], negexp5[[8]], linear5[[8]])){
        result = estimate_local_pars_sigma_zero(local_matrix = subnetwork, func_form = "negexponential", maxiter = 100, alpha=negexp5[[1]], beta=negexp5[[2]], sigma=negexp5[[3]], rho=negexp5[[4]], Psi_u = negexp5[[5]], Psi_v = negexp5[[6]], currentiter = preiter)
      }else{
        result = estimate_local_pars_sigma_zero(local_matrix = subnetwork, func_form = "linear", maxiter = 100, alpha=linear5[[1]], beta=linear5[[2]], sigma=linear5[[3]], rho=linear5[[4]], Psi_u = linear5[[5]], Psi_v = linear5[[6]], currentiter = preiter)
        
      }
      
    }
  }
result  = append(result, (proc.time() - local_start_time)[[3]])
}

  


load("Redo_DOROTHY_Acronym200.RData")



start_time =  Sys.time() 


library(foreach)
library(doParallel)
#cl <- makeCluster(6, type = "SOCK")
num_cores <- detectCores() - 1  # Use one less than the total number of cores
print(paste("num_cores", num_cores))
cl <- makeCluster(num_cores)
registerDoParallel(cl)


parallelOutput_sigma_zero <- foreach(k = 1:total_runs, .errorhandling="pass")%dopar%{
   if(mean(all_adjacency_matrix[[k]], na.rm=T)<.5){ 
   full_subnetwork_estimation_sigma_zero(all_adjacency_matrix[[k]])
  }else{
   full_subnetwork_estimation_sigma_zero(1-all_adjacency_matrix[[k]])
   }
  
}




total_runtime =  Sys.time() - start_time
print(total_runtime)


save.image(file = "Redo_DOROTHY_Acronym200Sigma0.RData")



parallelOutput_sigma_zero2 <- foreach(k = 1:total_runs, .errorhandling="pass")%dopar%{
if(length(parallelOutput_sigma_zero[[k]]) == 12){
parallelOutput_sigma_zero[[k]]
}else{ 
    if(mean(all_adjacency_matrix[[k]], na.rm=T)<.5){ 
   full_subnetwork_estimation_sigma_zero(all_adjacency_matrix[[k]])
  }else{
   full_subnetwork_estimation_sigma_zero(1-all_adjacency_matrix[[k]])
   }
  }
}
total_runtime =  Sys.time() - start_time
print(total_runtime)


save.image(file = "Redo_DOROTHY_Acronym200Sigma0.RData")













while(min(unlist(lapply(1:240, function(x){length(parallelOutput_sigma_zero2[[x]])})))!=12){
parallelOutput_sigma_zero = parallelOutput_sigma_zero2


parallelOutput_sigma_zero2 <- foreach(k = 1:total_runs, .errorhandling="pass")%dopar%{
if(length(parallelOutput_sigma_zero[[k]]) == 12){
parallelOutput_sigma_zero[[k]]
}else{ 
    if(mean(all_adjacency_matrix[[k]], na.rm=T)<.5){ 
   full_subnetwork_estimation_sigma_zero(all_adjacency_matrix[[k]])
  }else{
   full_subnetwork_estimation_sigma_zero(1-all_adjacency_matrix[[k]])
   }
  }
}
total_runtime =  Sys.time() - start_time
print(total_runtime)


save.image(file = "Redo_DOROTHY_Acronym200Sigma0.RData")


}




stopCluster(cl)
