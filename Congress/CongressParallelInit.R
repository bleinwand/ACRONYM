set.seed(2024)
library(rjson)
library(RSpectra)
library(caret)
library(foreach)
library(doParallel)


source("BackgroundFunctions.R")

#######
##MATRIX CREATION PART


myData <- fromJSON(file="congress_network_data.json" )


names(myData[[1]])


inlist = matrix(0, nrow= length(myData[[1]][[5]]), ncol=length(myData[[1]][[5]]))
outlist = matrix(0, nrow= length(myData[[1]][[5]]), ncol= length(myData[[1]][[5]]))
for(user in 1:nrow(inlist)){
  if(length(myData[[1]][[1]][[user]])>0){
    inlist[user, myData[[1]][[1]][[user]]+1] = 1}
  if(length(myData[[1]][[3]][[user]])>0){
    outlist[user, myData[[1]][[3]][[user]]+1] = 1
    
  }
}


total_list = inlist+outlist
total_list[total_list==2] = 1


adjacency_matrix = total_list




total_runs = 101
#all_sig_within =  vector(mode='list', length = total_runs)
#all_sig_between =  vector(mode='list', length = total_runs)
#all_node_groups = vector(mode='list', length = total_runs)
#all_node_sociabilities =  vector(mode='list', length = total_runs)
#all_node_types =  vector(mode='list', length = total_runs)
#all_network_params =  vector(mode='list', length = total_runs)
#all_probs_matrix =  vector(mode='list', length = total_runs)
#all_no_error_probs_matrix =  vector(mode='list', length = total_runs)
#all_adjacency_matrix =  vector(mode='list', length = total_runs)
all_estimated_node_groupsEig =  vector(mode='list', length = total_runs)
all_estimated_node_groupsU =  vector(mode='list', length = total_runs)
all_estimated_node_groupsV =  vector(mode='list', length = total_runs)

all_comm_pairsEig = vector(mode='list', length = total_runs)
all_comm_pairsU = vector(mode='list', length = total_runs)
all_comm_pairsV = vector(mode='list', length = total_runs)


VRepeats = rep(0, total_runs)
EigRepeats = rep(0, total_runs)

all_parallel_output =  vector(mode='list', length = total_runs)
all_run_times = vector(mode='list', length = total_runs)
#library(sinkR)


fold_contents = append(list(NA), replicate(n=10, createFolds(1:length(adjacency_matrix[upper.tri(adjacency_matrix)]), k=10)))

start_time = Sys.time()
for(run_no in 1:total_runs){ 
if(run_no == 1){
normalized_matrix = apply(adjacency_matrix, 2, function(x){(x- mean(x))/sd(x)})
svd_attempt = svd(normalized_matrix)
cos_sph2 = cosine(t(svd_attempt$u[,1:3])) #Need to check after the fact that 8 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "ward.D")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
all_estimated_node_groupsU[[run_no]] = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)

#Now also do for V!
cos_sph2 = cosine(t(svd_attempt$v[,1:3])) #Need to check after the fact that 8 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "ward.D")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
all_estimated_node_groupsV[[run_no]] = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)



eig_norm = eigen(normalized_matrix)
k_l =2
  sph_mat = Re(eig_norm$vectors)
  for(i in 1:nrow(sph_mat)){
    scaler = sqrt(sum(sph_mat[i,1:k_l]^2))
    sph_mat[i,1:k_l] = sph_mat[i,1:k_l]/scaler
  }  
  
  
  cos_sph = cosine(t(Re(sph_mat[,1:k_l])))
  hdist = hclust(as.dist(1- cos_sph), method = "ward.D")
  heightmin = which.max(diff(hdist$height))
  
  all_estimated_node_groupsEig[[run_no]]= cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)
 

}else{  
new_adjacency_matrix= adjacency_matrix
new_adjacency_matrix[lower.tri(new_adjacency_matrix)]=0
new_adjacency_matrix[upper.tri(new_adjacency_matrix)][fold_contents[[run_no]]]=NA
new_adjacency_matrix= new_adjacency_matrix+t(new_adjacency_matrix)

normalized_matrix = apply(new_adjacency_matrix, 2, function(x){(x- mean(x, na.rm=T))/sd(x, na.rm=T)})
Dineof = dineof(normalized_matrix)
svd_attempt = svd(Dineof$Xa)
cos_sph2 = cosine(t(svd_attempt$u[,1:3])) #Need to check after the fact that 3 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "ward.D")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
all_estimated_node_groupsU[[run_no]] = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)


#Now also do for V!
cos_sph2 = cosine(t(svd_attempt$v[,1:3])) #Need to check after the fact that 8 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "ward.D")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
all_estimated_node_groupsV[[run_no]] = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)

#Eigen Version
eig_norm = eigen(Dineof$Xa)
k_l =2
  sph_mat = Re(eig_norm$vectors)
  for(i in 1:nrow(sph_mat)){
    scaler = sqrt(sum(sph_mat[i,1:k_l]^2))
    sph_mat[i,1:k_l] = sph_mat[i,1:k_l]/scaler
  }  
  
  cos_sph = cosine(t(Re(sph_mat[,1:k_l])))
  hdist = hclust(as.dist(1- cos_sph), method = "ward.D")
  heightmin = which.max(diff(hdist$height))
   all_estimated_node_groupsEig[[run_no]]= cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)
}

num_groupsU = max(all_estimated_node_groupsU[[run_no]])
num_groupsV = max(all_estimated_node_groupsV[[run_no]])
num_groupsEig = max(all_estimated_node_groupsEig[[run_no]])


comm_pairsU = matrix(0, nrow= (num_groupsU*(num_groupsU+1)/2), ncol=2)
comm_pairs_iteratorU = 0
for(i in 1:num_groupsU){
  for(j in i:num_groupsU){
    comm_pairs_iteratorU = comm_pairs_iteratorU + 1
    comm_pairsU[comm_pairs_iteratorU,] = c(i, j)
    }
}

all_comm_pairsU[[run_no]] = comm_pairsU



comm_pairsV = matrix(0, nrow= (num_groupsV*(num_groupsV+1)/2), ncol=2)
comm_pairs_iteratorV = 0
for(i in 1:num_groupsV){
  for(j in i:num_groupsV){
    comm_pairs_iteratorV = comm_pairs_iteratorV + 1
    comm_pairsV[comm_pairs_iteratorV,] = c(i, j)
    }
}

all_comm_pairsV[[run_no]] = comm_pairsV



comm_pairsEig = matrix(0, nrow= (num_groupsEig*(num_groupsEig+1)/2), ncol=2)
comm_pairs_iteratorEig = 0
for(i in 1:num_groupsEig){
  for(j in i:num_groupsEig){
    comm_pairs_iteratorEig = comm_pairs_iteratorEig + 1
    comm_pairsEig[comm_pairs_iteratorEig,] = c(i, j)
    }
}

all_comm_pairsEig[[run_no]] = comm_pairsEig



if(all.equal(all_estimated_node_groupsV[[run_no]], all_estimated_node_groupsU[[run_no]])==1){
VRepeats[run_no] = 1
}


if(all.equal(all_estimated_node_groupsEig[[run_no]], all_estimated_node_groupsV[[run_no]])==1){
EigRepeats[run_no] = 2
}


if(all.equal(all_estimated_node_groupsEig[[run_no]], all_estimated_node_groupsU[[run_no]])==1){
EigRepeats[run_no] = 1
}






}




save.image(file = "CongressCommDet.RData")

