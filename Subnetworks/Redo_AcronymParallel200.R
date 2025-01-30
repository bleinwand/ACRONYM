set.seed(2024)
library(foreach)
library(doParallel)



source("BackgroundFunctions.R")

total_runs = 240
all_sig_within =  vector(mode='list', length = total_runs)
all_sig_between =  vector(mode='list', length = total_runs)
all_node_groups = vector(mode='list', length = total_runs)
all_node_sociabilities =  vector(mode='list', length = total_runs)
all_node_types =  vector(mode='list', length = total_runs)
all_network_params =  vector(mode='list', length = total_runs)
all_probs_matrix =  vector(mode='list', length = total_runs)
all_no_error_probs_matrix =  vector(mode='list', length = total_runs)
all_adjacency_matrix =  vector(mode='list', length = total_runs)
all_estimated_node_groups =  vector(mode='list', length = total_runs)
all_parallel_output =  vector(mode='list', length = total_runs)
all_run_times = vector(mode='list', length = total_runs)
all_comm_pairs = vector(mode='list', length = total_runs)


HForms = rep(c('normal', 'linear', 'exponential', 'negexponential'), each= 60)
alphaForms = rep(rep(c(1, .5, .2), each=20), 4) 
betaForms = rep(rep(c(0, .1, .2), each=20), 4)


for(run_no in 1:total_runs){ 
start_time = Sys.time()
##MATRIX CREATION PART


sig_within = .5 #We have 4 different network types, so we will include each of 6 sigma levels for each type
sig_between = .5 #this isn't used in this demo, but can be used later when we have multiple communities


nodes= 400
node_groups= rep(1:2, each=200) 
node_types = rep("uniform",2) 
node_sociabilities = c(sort(runif(200)), sort(runif(200)))

#sociability_generator = function(node_type, n){
#if(node_type == "uniform"){
#return( sort(runif(n=n)))
#}
#if(node_type == "gap"){
#
#gap_min = runif(1, .25, .7)
#gap_max = runif(1, gap_min+.1, min(gap_min +.3, .8) )
#values =  runif(n=n)
#for(i in intersect(which(values> gap_min), which(values<gap_max))){
#  if(runif(1)>(gap_min/(1-(gap_max - gap_min)))){
#    values[i]=runif(1, 0, gap_min)
#  }else{values[i]=runif(1, gap_max, 1)}
#
#
#}

#return(sort(values))
#}
#if(node_type == "beta"){
#a = runif(1, .6, 5/3)
#b = runif(1, .5, 2) * a
#return(sort(rbeta(n = n, shape1=a, shape2 =b)))
#}


#}



#node_sociabilities[which(node_groups==1)] = sociability_generator(node_type = node_types[1], n = sum(node_groups==1))
#node_sociabilities[which(node_groups==2)] = sociability_generator(node_type = node_types[2], n = sum(node_groups==2))
#node_sociabilities[which(node_groups==3)] = sociability_generator(node_type = node_types[3], n = sum(node_groups==3))
#node_sociabilities[which(node_groups==4)] = sociability_generator(node_type = node_types[4], n = sum(node_groups==4))
#node_sociabilities[which(node_groups==5)] = sociability_generator(node_type = node_types[5], n = sum(node_groups==5))
#node_sociabilities[which(node_groups==6)] = sociability_generator(node_type = node_types[6], n = sum(node_groups==6))

u_values = node_sociabilities #this is redundant/vestigial







#network_params = vector(mode='list', length = 1) 
#for(i in 1:length(network_params)){
#network_params[[i]]= vector(mode='list', length = 6)
#for(j in i:length(network_params)){
#network_params[[i]][[j]][[1]] =ifelse(runif(n = 1) < .1, 0, rbeta(1, 2,2) ) #alpha
#network_params[[i]][[j]][[2]] = runif(1, 0, 1-network_params[[i]][[j]][[1]]) #beta
#if(i == j){
#network_params[[i]][[j]][[3]] = 1
#network_params[[i]][[j]][[4]] = sample(c('p', 'n'), 1) #assoc
#} else{
#balcheck = runif(1)
#network_params[[i]][[j]][[3]] = ifelse(balcheck < .4, 1+rexp(1, rate= 5/3), ifelse(balcheck>.6, (1/(1+rexp(1, rate= 5/3))), 1))#balance
#network_params[[i]][[j]][[4]] = sample(c('p', 'n', 'pn', 'np'), 1) #assoc
#}
#network_params[[i]][[j]][[5]] = sample(c('normal', 'linear', 'exponential', 'negexponential'), 1) #function



#}
#}






probs_matrix = matrix(0, nrow=200, ncol=200)
no_error_probs_matrix = probs_matrix
adjacency_matrix = probs_matrix
for(i in 1:200){
  print(i)
  for(j in 1:200){
    
    probs_matrix[i, j] = edge_creator(u1 = node_sociabilities[i], u2 = node_sociabilities[200+j], func_form = HForms[run_no], assoc = "p",  intercept = betaForms[run_no], multiplier = alphaForms[run_no], balance = 1, sig= .5)
    
    no_error_probs_matrix[i, j] = zero_epsilon_edge_creator(u1 = node_sociabilities[i], u2 = node_sociabilities[200+j], func_form = HForms[run_no], assoc = "p",  intercept = betaForms[run_no], multiplier = alphaForms[run_no], balance = 1, sig= .5)
  
  adjacency_matrix[i,j] =rbinom(n=1,size=1, p =probs_matrix[i,j]) 
  }
}

#probs_matrix = probs_matrix +t(probs_matrix)
#no_error_probs_matrix = no_error_probs_matrix +t(no_error_probs_matrix)


#reordering = cbind(1:nodes, node_groups, node_sociabilities)

#reordering = reordering[order(reordering[,2], reordering[,3]),]


#filled.contour(no_error_probs_matrix[reordering[,1], reordering[,1]], color.palette = colorRampPalette(c("blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05))

#################################




#adjacency_matrix = adjacency_matrix + t(adjacency_matrix)


estimated_node_groups = node_groups #we're not doing community estimation here, so this is fictitious


#comm_pairs = matrix(0, nrow= (num_groups*(num_groups+1)/2), ncol=2)
#comm_pairs_iterator = 0
#for(i in 1:num_groups){
#  for(j in i:num_groups){
#    comm_pairs_iterator = comm_pairs_iterator + 1
#    comm_pairs[comm_pairs_iterator,] = c(i, j)
#    }
#}


all_sig_within[[run_no]] =  sig_within
all_sig_between[[run_no]] =  sig_between
all_node_groups[[run_no]] = node_groups
all_node_sociabilities[[run_no]] =  node_sociabilities
all_node_types[[run_no]] =  node_types
#all_network_params[[run_no]] =  network_params
all_probs_matrix[[run_no]] =  probs_matrix
all_no_error_probs_matrix[[run_no]] =  no_error_probs_matrix
all_adjacency_matrix[[run_no]] =  adjacency_matrix
all_estimated_node_groups[[run_no]] =  estimated_node_groups
#all_comm_pairs[[run_no]] = comm_pairs
}

save.image(file = "Redo_DOROTHY_Acronym200.RData")



##### THE ACTUAL ESTIMATION
num_cores <- detectCores() - 1  # Use one less than the total number of cores
print(paste("num_cores", num_cores))
cl <- makeCluster(num_cores)
registerDoParallel(cl)



parallelOutput <- foreach(k = 1:total_runs, .errorhandling="pass")%dopar%{
   if(mean(all_adjacency_matrix[[k]], na.rm=T)<.5){ 
   full_subnetwork_estimation(all_adjacency_matrix[[k]])
  }else{
   full_subnetwork_estimation(1-all_adjacency_matrix[[k]])
   }
  
}

total_runtime =  Sys.time() - start_time
print(total_runtime)


save.image(file = "Redo_DOROTHY_Acronym200.RData")

parallelOutput2 <- foreach(k = 1:total_runs, .errorhandling="pass")%dopar%{
if(length(parallelOutput[[k]]) == 12){
parallelOutput[[k]]
}else{ 
   if(mean(all_adjacency_matrix[[k]], na.rm=T)<.5){ 
   full_subnetwork_estimation(all_adjacency_matrix[[k]])
  }else{
   full_subnetwork_estimation(1-all_adjacency_matrix[[k]])
   }
  }
}
total_runtime =  Sys.time() - start_time
print(total_runtime)

save.image(file = "Redo_DOROTHY_Acronym200.RData")

while(min(unlist(lapply(1:240, function(x){length(parallelOutput2[[x]])})))==2){
parallelOutput = parallelOutput2

parallelOutput2 <- foreach(k = 1:total_runs, .errorhandling="pass")%dopar%{
if(length(parallelOutput[[k]]) == 12){
parallelOutput[[k]]
}else{ 
   if(mean(all_adjacency_matrix[[k]], na.rm=T)<.5){ 
   full_subnetwork_estimation(all_adjacency_matrix[[k]])
  }else{
   full_subnetwork_estimation(1-all_adjacency_matrix[[k]])
   }
  }
}
total_runtime =  Sys.time() - start_time
print(total_runtime)

save.image(file = "Redo_DOROTHY_Acronym200.RData")
}











parallelOutputNORMAL <- foreach(k = 1:total_runs, .errorhandling="pass")%dopar%{
if(parallelOutput[[k]][[7]] == "normal"){
parallelOutput[[k]]
}else{
if(mean(all_adjacency_matrix[[k]], na.rm=T)<.5){ 
   estimate_local_pars(local_matrix = all_adjacency_matrix[[k]], func_form = "normal", maxiter = 100)
  }else{
    estimate_local_pars(local_matrix = 1-all_adjacency_matrix[[k]], func_form = "normal", maxiter = 100)
   }
  } 
}

save.image(file = "Redo_DOROTHY_Acronym200.RData")

parallelOutputNORMAL2 <- foreach(k = 1:total_runs, .errorhandling="pass")%dopar%{
if(length(parallelOutputNORMAL[[k]]) == 12){
parallelOutputNORMAL[[k]]
}else{
if(mean(all_adjacency_matrix[[k]], na.rm=T)<.5){ 
   estimate_local_pars(local_matrix = all_adjacency_matrix[[k]], func_form = "normal", maxiter = 100)
  }else{
    estimate_local_pars(local_matrix = 1-all_adjacency_matrix[[k]], func_form = "normal", maxiter = 100)
   }
  } 
}
save.image(file = "Redo_DOROTHY_Acronym200.RData")



while(min(unlist(lapply(1:240, function(x){length(parallelOutputNORMAL2[[x]])})))==2){
parallelOutputNORMAL = parallelOutputNORMAL2

parallelOutputNORMAL2 <- foreach(k = 1:total_runs, .errorhandling="pass")%dopar%{
if(length(parallelOutputNORMAL[[k]]) == 12){
parallelOutputNORMAL[[k]]
}else{
if(mean(all_adjacency_matrix[[k]], na.rm=T)<.5){ 
   estimate_local_pars(local_matrix = all_adjacency_matrix[[k]], func_form = "normal", maxiter = 100)
  }else{
    estimate_local_pars(local_matrix = 1-all_adjacency_matrix[[k]], func_form = "normal", maxiter = 100)
   }
  } 
}

total_runtime =  Sys.time() - start_time
print(total_runtime)

save.image(file = "Redo_DOROTHY_Acronym200.RData")

}




stopCluster(cl)





