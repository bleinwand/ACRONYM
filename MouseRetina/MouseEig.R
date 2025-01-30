set.seed(2024)
library(rjson)
library(RSpectra)
library(caret)
library(foreach)
library(doParallel)



source("BackgroundFunctions.R")

load("MouseCommDet.RData")

all_estimated_node_groups = all_estimated_node_groupsEig
all_comm_pairs = all_comm_pairsEig 

#cl <- makeCluster(6, type = "SOCK")
num_cores <- detectCores() - 1  # Use one less than the total number of cores
print(paste("num_cores", num_cores))
cl <- makeCluster(num_cores)
registerDoParallel(cl)

start_time = Sys.time() 


parallelOutput <- foreach(run = 1:total_runs)%:%
foreach(k = 1:nrow(all_comm_pairs[[run]]), .errorhandling="pass")%dopar%{
if(EigRepeats[run]>0){
vector(mode='list', length = 12)
}else{

    if(run == 1){
    new_adjacency_matrix = adjacency_matrix
    }else{
    new_adjacency_matrix= adjacency_matrix
    new_adjacency_matrix[lower.tri(new_adjacency_matrix)]=0
    new_adjacency_matrix[upper.tri(new_adjacency_matrix)][fold_contents[[run]]]=NA
    new_adjacency_matrix= new_adjacency_matrix+t(new_adjacency_matrix)
    }
   if(mean(new_adjacency_matrix[which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])], na.rm=T)<.5){ 
   full_subnetwork_estimation(new_adjacency_matrix[which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
  }else{
   full_subnetwork_estimation(1-new_adjacency_matrix[which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
   }
 } 
}


total_runtime =  Sys.time() - start_time
print(total_runtime)

save.image(file = "MouseFitEig.RData")


parallelOutput2 <- foreach(run = 1:total_runs)%:%
foreach(k = 1:nrow(all_comm_pairs[[run]]), .errorhandling="pass")%dopar%{
if(length(parallelOutput[[run]][[k]])==12){
parallelOutput[[run]][[k]]
}else{

if(run == 1){
    new_adjacency_matrix = adjacency_matrix
    }else{
    new_adjacency_matrix= adjacency_matrix
    new_adjacency_matrix[lower.tri(new_adjacency_matrix)]=0
    new_adjacency_matrix[upper.tri(new_adjacency_matrix)][fold_contents[[run]]]=NA
    new_adjacency_matrix= new_adjacency_matrix+t(new_adjacency_matrix)
    }

   if(mean(new_adjacency_matrix[which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])], na.rm=T)<.5){ 
   full_subnetwork_estimation(new_adjacency_matrix[which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
  }else{
   full_subnetwork_estimation(1-new_adjacency_matrix[which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
   }
  
}
}

total_runtime =  Sys.time() - start_time
print(total_runtime)


save.image(file = "MouseFitEig.RData")








#Now just keep updating as necessary
while(min(unlist(lapply(1:length(unlist(parallelOutput2, recursive=F)), function(x){length(unlist(parallelOutput2, recursive=F)[[x]])})))==2){

parallelOutput = parallelOutput2


parallelOutput2 <- foreach(run = 1:total_runs)%:%
foreach(k = 1:nrow(all_comm_pairs[[run]]), .errorhandling="pass")%dopar%{
if(length(parallelOutput[[run]][[k]])==12){
parallelOutput[[run]][[k]]
}else{

if(run == 1){
    new_adjacency_matrix = adjacency_matrix
    }else{
    new_adjacency_matrix= adjacency_matrix
    new_adjacency_matrix[lower.tri(new_adjacency_matrix)]=0
    new_adjacency_matrix[upper.tri(new_adjacency_matrix)][fold_contents[[run]]]=NA
    new_adjacency_matrix= new_adjacency_matrix+t(new_adjacency_matrix)
    }

   if(mean(new_adjacency_matrix[which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])], na.rm=T)<.5){ 
   full_subnetwork_estimation(new_adjacency_matrix[which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
  }else{
   full_subnetwork_estimation(1-new_adjacency_matrix[which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,1]), which(all_estimated_node_groups[[run]]==all_comm_pairs[[run]][k,2])])
   }
  
}
}

total_runtime =  Sys.time() - start_time
print(total_runtime)


save.image(file = "MouseFitEig.RData")
}







stopCluster(cl)










reconstructedEstimates = vector(mode='list', length = total_runs)
for(run in 1:total_runs){
if(EigRepeats[run]>0){
reconstructedEstimates[[run]] = NA
}else{
if(run == 1){
    new_adjacency_matrix = adjacency_matrix
    }else{
    new_adjacency_matrix= adjacency_matrix
    new_adjacency_matrix[lower.tri(new_adjacency_matrix)]=0
    new_adjacency_matrix[upper.tri(new_adjacency_matrix)][fold_contents[[run]]]=NA
    new_adjacency_matrix= new_adjacency_matrix+t(new_adjacency_matrix)
    }
print(run)
reconstructedEstimates[[run]] = 0*adjacency_matrix
for(k in 1:length(parallelOutput2[[run]])){
if(mean(new_adjacency_matrix[which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 1]), which(all_estimated_node_groups[[run]] == all_comm_pairs[[run]][k, 2])], na.rm=T)<.5){ 
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
 
 }
 
	save.image(file = "MouseFitEig.RData")

 
 
 
 
#for(run in 1:total_runs){
#if(EigRepeats[run]==1){
#reconstructedEstimates[[run]] = reconstructedEstimatesU[[run]] 
#} 
#if(EigRepeats[run]==2){
#reconstructedEstimates[[run]] = reconstructedEstimatesV[[run]] 
#} 


# }
 
 
#BlindreconstructedEstimates = vector(mode='list', length = 10)
#for(blinder in 1:10){
#BlindreconstructedEstimates[[blinder]] =  0*adjacency_matrix
#for(k in (2+(10*(blinder-1))):((10*blinder)+1)){
#BlindreconstructedEstimates[[blinder]][upper.tri(BlindreconstructedEstimates[[blinder]])][fold_contents[[k]]] = reconstructedEstimates[[k]][upper.tri(reconstructedEstimates[[k]])][fold_contents[[k]]]
#}
# BlindreconstructedEstimates[[blinder]] = BlindreconstructedEstimates[[blinder]] + t(BlindreconstructedEstimates[[blinder]])
#}
# save.image(file = "MouseFitEig.RData")