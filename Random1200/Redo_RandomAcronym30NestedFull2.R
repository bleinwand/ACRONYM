set.seed(2024)
library(foreach)
library(doParallel)

start_time = Sys.time()

source("BackgroundFunctions.R")


############################
##Defining network parameters



total_runs = 30
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




############################
##MATRIX CREATION PART

for(run_no in 1:total_runs){ 

sig_within = rexp(1, rate=3) #We have 4 different network types, so we will include each of 6 sigma levels for each type
sig_between = rexp(1,rate=3) #this isn't used in this demo, but can be used later when we have multiple communities


nodes= 1200
node_groups= sort(sample(1:6, nodes, replace=T))
node_types = sample(c("uniform", "gap", "beta"), 6, prob = c(.4, .3, .3), replace=T) 
node_sociabilities = runif(nodes)

sociability_generator = function(node_type, n){
if(node_type == "uniform"){
return( sort(runif(n=n)))
}
if(node_type == "gap"){

gap_min = runif(1, .25, .7)
gap_max = runif(1, gap_min+.1, min(gap_min +.3, .8) )
values =  runif(n=n)
for(i in intersect(which(values> gap_min), which(values<gap_max))){
  if(runif(1)>(gap_min/(1-(gap_max - gap_min)))){
    values[i]=runif(1, 0, gap_min)
  }else{values[i]=runif(1, gap_max, 1)}


}

return(sort(values))
}
if(node_type == "beta"){
a = runif(1, .6, 5/3)
b = runif(1, .5, 2) * a
return(sort(rbeta(n = n, shape1=a, shape2 =b)))
}


}



node_sociabilities[which(node_groups==1)] = sociability_generator(node_type = node_types[1], n = sum(node_groups==1))
node_sociabilities[which(node_groups==2)] = sociability_generator(node_type = node_types[2], n = sum(node_groups==2))
node_sociabilities[which(node_groups==3)] = sociability_generator(node_type = node_types[3], n = sum(node_groups==3))
node_sociabilities[which(node_groups==4)] = sociability_generator(node_type = node_types[4], n = sum(node_groups==4))
node_sociabilities[which(node_groups==5)] = sociability_generator(node_type = node_types[5], n = sum(node_groups==5))
node_sociabilities[which(node_groups==6)] = sociability_generator(node_type = node_types[6], n = sum(node_groups==6))

u_values = node_sociabilities #this is redundant/vestigial







network_params = vector(mode='list', length = 6) 
for(i in 1:length(network_params)){
network_params[[i]]= vector(mode='list', length = 6)
for(j in i:length(network_params)){
network_params[[i]][[j]][[1]] =ifelse(runif(n = 1) < .1, 0, rbeta(1, 2,2) ) #alpha
network_params[[i]][[j]][[2]] = runif(1, 0, 1-network_params[[i]][[j]][[1]]) #beta
if(i == j){
network_params[[i]][[j]][[3]] = 1
network_params[[i]][[j]][[4]] = sample(c('p', 'n'), 1) #assoc
} else{
balcheck = runif(1)
network_params[[i]][[j]][[3]] = ifelse(balcheck < .4, 1+rexp(1, rate= 5/3), ifelse(balcheck>.6, (1/(1+rexp(1, rate= 5/3))), 1))#balance
network_params[[i]][[j]][[4]] = sample(c('p', 'n', 'pn', 'np'), 1) #assoc
}
network_params[[i]][[j]][[5]] = sample(c('normal', 'linear', 'exponential', 'negexponential'), 1) #function



}
}






probs_matrix = matrix(0, nrow=nodes, ncol=nodes)
no_error_probs_matrix = probs_matrix

for(i in 1:(nodes-1)){
  print(i)
  for(j in (i+1):nodes){
    
    probs_matrix[i, j] = edge_creator(u1 = node_sociabilities[i], u2 = node_sociabilities[j], func_form = network_params[[node_groups[i]]][[node_groups[j]]][[5]], assoc = network_params[[node_groups[i]]][[node_groups[j]]][[4]],  intercept = network_params[[node_groups[i]]][[node_groups[j]]][[2]], multiplier = network_params[[node_groups[i]]][[node_groups[j]]][[1]], balance = network_params[[node_groups[i]]][[node_groups[j]]][[3]], sig= ifelse(node_groups[i]==node_groups[j], sig_within, sig_between))
    no_error_probs_matrix[i, j] = zero_epsilon_edge_creator(u1 = node_sociabilities[i], u2 = node_sociabilities[j], func_form = network_params[[node_groups[i]]][[node_groups[j]]][[5]], assoc = network_params[[node_groups[i]]][[node_groups[j]]][[4]],  intercept = network_params[[node_groups[i]]][[node_groups[j]]][[2]], multiplier = network_params[[node_groups[i]]][[node_groups[j]]][[1]], balance = network_params[[node_groups[i]]][[node_groups[j]]][[3]], sig= ifelse(node_groups[i]==node_groups[j], sig_within, sig_between))
    
  }
}

probs_matrix = probs_matrix +t(probs_matrix)
no_error_probs_matrix = no_error_probs_matrix +t(no_error_probs_matrix)


#reordering = cbind(1:nodes, node_groups, node_sociabilities)

#reordering = reordering[order(reordering[,2], reordering[,3]),]


#filled.contour(no_error_probs_matrix[reordering[,1], reordering[,1]], color.palette = colorRampPalette(c("blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05))

#################################



adjacency_matrix = matrix(0, nrow=nodes, ncol=nodes)
for(i in 1:(nodes-1)){
  for(j in (i+1):nodes){
    adjacency_matrix[i,j] =rbinom(n=1,size=1, p =probs_matrix[i,j])
  }
}

adjacency_matrix = adjacency_matrix + t(adjacency_matrix)








############################
#Community Detection Part - note: we use 8 as the number of vectors to keep (k_l), but in practice, the user will need to look at the eigenvalues of eig_norm and look for a gap. 


normalized_matrix = apply(adjacency_matrix, 2, function(x){(x- mean(x))/sd(x)})
eig_norm = eigen(normalized_matrix)
#plot(abs(eig_norm$values))


k_l = 8 #Need to check after the fact that 8 is correct. 
sph_mat = eig_norm$vectors
for(i in 1:nrow(sph_mat)){
  scaler = sqrt(sum(sph_mat[i,1:k_l]^2))
  sph_mat[i,1:k_l] = sph_mat[i,1:k_l]/scaler
  
}


cos_sph = cosine(t(sph_mat[,1:k_l]))
hdist = hclust(as.dist(1- cos_sph), method = "single")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
estimated_node_groups = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)
num_groups = max(estimated_node_groups)
print(num_groups)








comm_pairs = matrix(0, nrow= (num_groups*(num_groups+1)/2), ncol=2)
comm_pairs_iterator = 0
for(i in 1:num_groups){
  for(j in i:num_groups){
    comm_pairs_iterator = comm_pairs_iterator + 1
    comm_pairs[comm_pairs_iterator,] = c(i, j)
    }
}



############################
#Recording all the information about our generating processes, matrices, and community detection results

all_sig_within[[run_no]] =  sig_within
all_sig_between[[run_no]] =  sig_between
all_node_groups[[run_no]] = node_groups
all_node_sociabilities[[run_no]] =  node_sociabilities
all_node_types[[run_no]] =  node_types
all_network_params[[run_no]] =  network_params
all_probs_matrix[[run_no]] =  probs_matrix
all_no_error_probs_matrix[[run_no]] =  no_error_probs_matrix
all_adjacency_matrix[[run_no]] =  adjacency_matrix
all_estimated_node_groups[[run_no]] =  estimated_node_groups
all_comm_pairs[[run_no]] = comm_pairs
}



save.image(file = "Redo_DOROTHY_RandomAcronym30NestedFULL2.RData")
############################
##### THE ACTUAL ESTIMATION Parallelizing (within parallelizing) everything

library(foreach)
library(doParallel)
#cl <- makeCluster(6, type = "SOCK")
num_cores <- detectCores() - 1  # Use one less than the total number of cores
print(paste("num_cores", num_cores))
cl <- makeCluster(num_cores)
registerDoParallel(cl)



#full_network_predictions = vector(mode="list", length = num_groups)
#for(node_group in 1:num_groups){full_network_predictions[[node_group]] = vector(mode="list", length = num_groups)
#}

#estimated_full_network_no_error_probs = matrix(0, nrow= nrow(adjacency_matrix), ncol(adjacency_matrix))
#foreach(first_group = 1:num_groups, .errorhandling="remove")%dopar% {
  #for(second_group in first_group:num_groups){
    #full_network_predictions[[first_group]][[second_group]] = full_subnetwork_estimation(adjacency_matrix[which(estimated_node_groups==first_group), which(estimated_node_groups==second_group)])
    #estimated_full_network_no_error_probs[which(estimated_node_groups==first_group), which(estimated_node_groups==second_group)] =  full_network_predictions[[first_group]][[second_group]][[11]]
    
    #if(second_group>first_group){
      #estimated_full_network_no_error_probs[which(estimated_node_groups==second_group), which(estimated_node_groups==first_group)] = t(matrix(full_network_predictions[[first_group]][[second_group]][[11]], nrow=length(which(estimated_node_groups==first_group))))
      
    #}
    
  #}
#}


parallelOutput <- foreach(run = 1:total_runs)%:%
foreach(k = 1:nrow(all_comm_pairs[[run]]), .errorhandling="pass")%dopar%{
   if(mean(all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])], na.rm=T)<.5){ 
   full_subnetwork_estimation(all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
  }else{
   full_subnetwork_estimation(1-all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
   }
  
}



total_runtime =  Sys.time() - start_time
print(total_runtime)

save.image(file = "Redo_DOROTHY_RandomAcronym30NestedFULL2.RData")

############################
#we do a second estimation try just to clean up any issues that may have come up, but only fix those with errors

parallelOutput2 <- foreach(run = 1:total_runs)%:%
foreach(k = 1:nrow(all_comm_pairs[[run]]), .errorhandling="pass")%dopar%{
if(length(parallelOutput[[run]][[k]])==12){
parallelOutput[[run]][[k]]
}else{
   if(mean(all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])], na.rm=T)<.5){ 
   full_subnetwork_estimation(all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
  }else{
   full_subnetwork_estimation(1-all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
   }
  
}
}

total_runtime =  Sys.time() - start_time
print(total_runtime)





while(min(unlist(lapply(1:length(unlist(parallelOutput2, recursive=F)), function(x){length(unlist(parallelOutput2, recursive=F)[[x]])})))==2){

parallelOutput = parallelOutput2 

parallelOutput2 <- foreach(run = 1:total_runs)%:%
foreach(k = 1:nrow(all_comm_pairs[[run]]), .errorhandling="pass")%dopar%{
if(length(parallelOutput[[run]][[k]])==12){
parallelOutput[[run]][[k]]
}else{
   if(mean(all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])], na.rm=T)<.5){ 
   full_subnetwork_estimation(all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
  }else{
   full_subnetwork_estimation(1-all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
   }
  
}
}


total_runtime =  Sys.time() - start_time
print(total_runtime)

save.image(file = "Redo_DOROTHY_RandomAcronym30NestedFULL2.RData")
}




stopCluster(cl)

############################
#The above doesn't give you our estimates of the no_error_prob_matrix values, so here, we collect everything into the easy to compare list of matrices reconstructedEstimates

reconstructedEstimates = vector(mode='list', length = total_runs)
for(run in 1:total_runs){
print(run)
reconstructedEstimates[[run]] = 0*all_adjacency_matrix[[run]]
for(k in 1:length(parallelOutput2[[run]])){
if(mean(all_adjacency_matrix[[run]][which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 1]), which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 2])], na.rm=T)<.5){ 
reconstructedEstimates[[run]][which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 1]), which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 2])] = parallelOutput2[[run]][[k]][[9]] 
if(all_comm_pairs[[run]][k, 2] != all_comm_pairs[[run]][k, 1]){
reconstructedEstimates[[run]][which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 2]), which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 1])] = t(parallelOutput2[[run]][[k]][[9]]) 
}
	}else{
	reconstructedEstimates[[run]][which(all_estimated_node_groups[[run]] ==  all_comm_pairs[[run]][k, 1]), which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 2])] = 1 - parallelOutput2[[run]][[k]][[9]] 
if(all_comm_pairs[[run]][k, 2] != all_comm_pairs[[run]][k, 1]){
reconstructedEstimates[[run]][which(all_estimated_node_groups[[run]] ==  all_comm_pairs[[run]][k, 2]), which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 1])] = 1- t(parallelOutput2[[run]][[k]][[9]]) 
	
	
	}
	}
	
	
	}
	}
	
	
	
all_matrix_MSEs = unlist(lapply(1:30, function(x){mean((reconstructedEstimates[[x]][lower.tri(reconstructedEstimates[[x]])] - all_no_error_probs_matrix[[x]][lower.tri(all_no_error_probs_matrix[[x]])])^2)}))




save.image(file = "Redo_DOROTHY_RandomAcronym30NestedFULL2.RData")

