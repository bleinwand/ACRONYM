load("Redo_DOROTHY_RandomAcronym30NestedFULL2.RData")

all_comm_pairs_SVD =  vector(mode='list', length = total_runs)


start_time = Sys.time()



############################
#Community Detection Part - note: we use 8 as the number of vectors to keep (k_l), but in practice, the user will need to look at the eigenvalues of eig_norm and look for a gap. 



library(foreach)
library(doParallel)
library(mclust)



#cl <- makeCluster(6, type = "SOCK")
num_cores <- detectCores() - 1  # Use one less than the total number of cores
print(paste("num_cores", num_cores))
cl <- makeCluster(num_cores)
registerDoParallel(cl)



all_estimated_node_groups_SVD <- foreach(run_no = 1:total_runs)%dopar%{
normalized_matrix = apply(all_adjacency_matrix[[run_no]], 2, function(x){(x- mean(x))/sd(x)})
svd_attempt = svd(normalized_matrix)
cos_sph2 = cosine(t(svd_attempt$v[,1:8])) #Need to check after the fact that 8 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "single")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)
}

 





for(run_no in 1:total_runs){ 
num_groups = max(all_estimated_node_groups_SVD[[run_no]])
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

all_comm_pairs_SVD[[run_no]] = comm_pairs
}

############################
##### THE ACTUAL ESTIMATION Parallelizing (within parallelizing) everything



#Don't re-estimate when the SVD approach gives the same communities as the eigenvector approach
parallelOutput_SVD <- foreach(run = 1:total_runs)%:%
foreach(k = 1:nrow(all_comm_pairs_SVD[[run]]), .errorhandling="pass")%dopar%{
if(all.equal(all_estimated_node_groups_SVD[[run]], all_estimated_node_groups[[run]])==1){
vector(mode='list', length = 12)
}else{ 
   if(mean(all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,2])], na.rm=T)<.5){ 
   full_subnetwork_estimation(all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,2])])
  }else{
   full_subnetwork_estimation(1-all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,2])])
   }
  
}
}


total_runtime =  Sys.time() - start_time
print(total_runtime)

save.image(file = "Redo_SingularValueTestV.RData")

############################
#we do a second estimation try just to clean up any issues that may have come up, but only fix those with errors

parallelOutput_SVD2 <- foreach(run = 1:total_runs)%:%
foreach(k = 1:nrow(all_comm_pairs_SVD[[run]]), .errorhandling="pass")%dopar%{
if(length(parallelOutput_SVD[[run]][[k]])==12){
parallelOutput_SVD[[run]][[k]]
}else{
   if(mean(all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,2])], na.rm=T)<.5){ 
   full_subnetwork_estimation(all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,2])])
  }else{
   full_subnetwork_estimation(1-all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,2])])
   }
  
}
}

total_runtime =  Sys.time() - start_time
print(total_runtime)
save.image(file = "Redo_SingularValueTestV.RData")



while(min(unlist(lapply(1:length(unlist(parallelOutput_SVD2, recursive=F)), function(x){length(unlist(parallelOutput_SVD2, recursive=F)[[x]])})))==2){

parallelOutput_SVD = parallelOutput_SVD2

parallelOutput_SVD2 <- foreach(run = 1:total_runs)%:%
foreach(k = 1:nrow(all_comm_pairs_SVD[[run]]), .errorhandling="pass")%dopar%{
if(length(parallelOutput_SVD[[run]][[k]])==12){
parallelOutput_SVD[[run]][[k]]
}else{
   if(mean(all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,2])], na.rm=T)<.5){ 
   full_subnetwork_estimation(all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,2])])
  }else{
   full_subnetwork_estimation(1-all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs_SVD[[run]][k,2])])
   }
  
}
}

total_runtime =  Sys.time() - start_time
print(total_runtime)
save.image(file = "Redo_SingularValueTestV.RData")

}









stopCluster(cl)

############################
#The above doesn't give you our estimates of the no_error_prob_matrix values, so here, we collect everything into the easy to compare list of matrices reconstructedEstimates, except those we already have from the eigenvector method

reconstructedEstimates_SVD = vector(mode='list', length = total_runs)
for(run in 1:total_runs){
print(run)
if(all.equal(all_estimated_node_groups_SVD[[run]], all_estimated_node_groups[[run]])!=1){
reconstructedEstimates_SVD[[run]] = 0*all_adjacency_matrix[[run]]
for(k in 1:length(parallelOutput_SVD2[[run]])){
if(mean(all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]] == all_comm_pairs_SVD[[run]][k, 1]), which(all_estimated_node_groups_SVD[[run]] == all_comm_pairs_SVD[[run]][k, 2])], na.rm=T)<.5){ 
reconstructedEstimates_SVD[[run]][which(all_estimated_node_groups_SVD[[run]] == all_comm_pairs_SVD[[run]][k, 1]), which(all_estimated_node_groups_SVD[[run]] == all_comm_pairs_SVD[[run]][k, 2])] = parallelOutput_SVD2[[run]][[k]][[9]] 
if(all_comm_pairs_SVD[[run]][k, 2] != all_comm_pairs_SVD[[run]][k, 1]){
reconstructedEstimates_SVD[[run]][which(all_estimated_node_groups_SVD[[run]] == all_comm_pairs_SVD[[run]][k, 2]), which(all_estimated_node_groups_SVD[[run]] == all_comm_pairs_SVD[[run]][k, 1])] = t(parallelOutput_SVD2[[run]][[k]][[9]]) 
}
	}else{
	reconstructedEstimates_SVD[[run]][which(all_estimated_node_groups_SVD[[run]] ==  all_comm_pairs_SVD[[run]][k, 1]), which(all_estimated_node_groups_SVD[[run]] == all_comm_pairs_SVD[[run]][k, 2])] = 1 - parallelOutput_SVD2[[run]][[k]][[9]] 
if(all_comm_pairs_SVD[[run]][k, 2] != all_comm_pairs_SVD[[run]][k, 1]){
	reconstructedEstimates_SVD[[run]][which(all_estimated_node_groups_SVD[[run]] == all_comm_pairs_SVD[[run]][k, 2]), which(all_estimated_node_groups_SVD[[run]] == all_comm_pairs_SVD[[run]][k, 1])] = 1- t(parallelOutput_SVD2[[run]][[k]][[9]]) 

	
	}
	}
	
	
	}
	}
	}
	
	
all_matrix_MSEs_SVD = unlist(lapply(1:30, function(x){mean((reconstructedEstimates_SVD[[x]][lower.tri(all_no_error_probs_matrix[[x]])] - all_no_error_probs_matrix[[x]][lower.tri(all_no_error_probs_matrix[[x]])])^2)}))




save.image(file = "Redo_SingularValueTestV.RData")

##all_parallel_output[[run_no]] =  parallelOutput
##all_run_times[[run_no]] = total_runtime

#for(run in 1:total_runs){
#for(k in 1:nrow(all_comm_pairs[[run]])){
#if(mean(all_adjacency_matrix[[run]][which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups_SVD[[run]]==all_comm_pairs[[run]][k,2])], na.rm=T)>=.5){ 
#parallelOutput[[run]][[k]][[2]] = 1 - parallelOutput[[run]][[k]][[2]] - parallelOutput[[run]][[k]][[1]] 
#parallelOutput[[run]][[k]] = 1 - parallelOutput[[run]][[k]][[9]] 
#parallelOutput[[run]][[k]][[10]] = 1 - parallelOutput[[run]][[k]][[10]] 
#parallelOutput[[run]][[k]][[11]] = 1 - parallelOutput[[run]][[k]][[11]] 
#}
#}
#}

#save.image(file = "DOROTHY_RandomAcronym30NestedFULL.RData")

#for(k in 1:nrow(comm_pairs)){
#full_network_predictions[[comm_pairs[k,1]]][[comm_pairs[k,2]]] = parallelOutput[[k]]
#}









