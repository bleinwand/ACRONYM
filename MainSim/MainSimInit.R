set.seed(2024)

library(RSpectra)
library(caret)
library(foreach)
library(doParallel)


source("BackgroundFunctions.R")



#######
##MATRIX CREATION PART

node_sociabilities = runif(1200)
nodes= 1200

node_sociabilities[1:200] = sort(rbeta(200, shape1 = 1.5, shape2=.9)) #used to be sort(rbeta(num_nodes, shape1 = 2.1, shape2=1.2))
node_sociabilities[201:400] = sort(rbeta(200, shape1 = .8, shape2=1.4))# used to be sort(rbeta(num_nodes2, shape1 = 1.6, shape2=2.8))
for(i in intersect(intersect(which(node_sociabilities>.4), which(node_sociabilities<.6)), 401:600)){
  if(runif(1)>.5){
    node_sociabilities[i]=runif(1, 0, .4)
  }else{node_sociabilities[i]=runif(1, .6, 1)}
  
}
node_sociabilities[401:600]= sort(node_sociabilities[401:600])

for(i in intersect(intersect(which(node_sociabilities>.55), which(node_sociabilities<.8)), 601:800)){
  if(runif(1)>(.55/.75)){
    node_sociabilities[i]=runif(1, 0, .55)
  }else{node_sociabilities[i]=runif(1, .8, 1)}
  
}
node_sociabilities[601:800]= sort(node_sociabilities[601:800])
node_sociabilities[801:1000] = sort(node_sociabilities[801:1000])
node_sociabilities[1001:1200] = sort(node_sociabilities[1001:1200])


u_values = node_sociabilities #this is redundant/vestigial
node_groups = rep(1:6, each=200)





sig_within = .3 #We have 4 different network types, so we will include each of 6 sigma levels for each type
sig_between = .4 #this isn't used in this demo, but can be used later when we have multiple communities

network_params = vector(mode='list', length = 6) 
for(i in 1:length(network_params)){network_params[[i]]= vector(mode='list', length = 6)}
network_params[[1]][[1]][[1]] =.6 #alpha
network_params[[1]][[1]][[2]] = .2 #beta
network_params[[1]][[1]][[3]] = 1#balance
network_params[[1]][[1]][[4]] =  "p" #assoc
network_params[[1]][[1]][[5]] = 'normal' #function

network_params[[1]][[2]][[1]] = 1#alpha
network_params[[1]][[2]][[2]] = 0#beta
network_params[[1]][[2]][[3]] = 1.3 #balance
network_params[[1]][[2]][[4]] = 'n' #assoc
network_params[[1]][[2]][[5]] = 'negexponential' #function

network_params[[1]][[3]][[1]] =.8 #alpha
network_params[[1]][[3]][[2]] = .1#beta
network_params[[1]][[3]][[3]] = .5#balance
network_params[[1]][[3]][[4]] = 'np' #assoc
network_params[[1]][[3]][[5]] = 'linear' #function


network_params[[1]][[4]][[1]] = .7#alpha
network_params[[1]][[4]][[2]] = .3#beta
network_params[[1]][[4]][[3]] = 2.5 #balance
network_params[[1]][[4]][[4]] = 'n' #assoc
network_params[[1]][[4]][[5]] = 'normal'#function


network_params[[1]][[5]][[1]] = .8 #alpha
network_params[[1]][[5]][[2]] = .2#beta
network_params[[1]][[5]][[3]] = .8#balance
network_params[[1]][[5]][[4]] = 'p'#assoc
network_params[[1]][[5]][[5]] = 'exponential' #function

network_params[[1]][[6]][[1]] = 0 #alpha
network_params[[1]][[6]][[2]] = .6 #beta
network_params[[1]][[6]][[3]] = 1#balance
network_params[[1]][[6]][[4]] = 'p'#assoc
network_params[[1]][[6]][[5]] = 'normal'#function


network_params[[2]][[2]][[1]] = .6#alpha
network_params[[2]][[2]][[2]] = .4#beta
network_params[[2]][[2]][[3]] = 1#balance
network_params[[2]][[2]][[4]] = 'p' #assoc
network_params[[2]][[2]][[5]] = 'linear'#function

network_params[[2]][[3]][[1]] = .9#alpha
network_params[[2]][[3]][[2]] = .05#beta
network_params[[2]][[3]][[3]] = .4#balance
network_params[[2]][[3]][[4]] = 'pn'#assoc
network_params[[2]][[3]][[5]] = 'exponential'#function


network_params[[2]][[4]][[1]] = .5#alpha
network_params[[2]][[4]][[2]] = .25#beta
network_params[[2]][[4]][[3]] = 1.1#balance
network_params[[2]][[4]][[4]] = 'p'#assoc
network_params[[2]][[4]][[5]] = 'linear' #function


network_params[[2]][[5]][[1]] = .6#alpha
network_params[[2]][[5]][[2]] = .3#beta
network_params[[2]][[5]][[3]] = 1.3#balance
network_params[[2]][[5]][[4]] = 'pn'#assoc
network_params[[2]][[5]][[5]] = 'negexponential' #function

network_params[[2]][[6]][[1]] = .6#alpha
network_params[[2]][[6]][[2]] = 0#beta
network_params[[2]][[6]][[3]] = 3#balance
network_params[[2]][[6]][[4]] = 'np'#assoc
network_params[[2]][[6]][[5]] = 'negexponential'#function


network_params[[3]][[3]][[1]] = .85#alpha
network_params[[3]][[3]][[2]] = .1#beta
network_params[[3]][[3]][[3]] = 1 #balance
network_params[[3]][[3]][[4]] = 'p'#assoc
network_params[[3]][[3]][[5]] = 'negexponential'#function

network_params[[3]][[4]][[1]] = 0#alpha
network_params[[3]][[4]][[2]] = .3#beta
network_params[[3]][[4]][[3]] = 1#balance
network_params[[3]][[4]][[4]] = 'n'#assoc
network_params[[3]][[4]][[5]] = 'exponential'#function


network_params[[3]][[5]][[1]] = 1#alpha
network_params[[3]][[5]][[2]] = 0#beta
network_params[[3]][[5]][[3]] = 2.2#balance
network_params[[3]][[5]][[4]] = 'np'#assoc
network_params[[3]][[5]][[5]] = 'negexponential'#function


network_params[[3]][[6]][[1]] = .5#alpha
network_params[[3]][[6]][[2]] = 0#beta
network_params[[3]][[6]][[3]] = .75#balance
network_params[[3]][[6]][[4]] = 'pn'#assoc
network_params[[3]][[6]][[5]] = 'exponential'#function






network_params[[4]][[4]][[1]] =1 #alpha
network_params[[4]][[4]][[2]] = 0#beta
network_params[[4]][[4]][[3]] = 1#balance
network_params[[4]][[4]][[4]] = 'n'#assoc
network_params[[4]][[4]][[5]] = 'normal'#function


network_params[[4]][[5]][[1]] = .9#alpha
network_params[[4]][[5]][[2]] = 0#beta
network_params[[4]][[5]][[3]] = 1.4#balance
network_params[[4]][[5]][[4]] = 'p'#assoc
network_params[[4]][[5]][[5]] = 'normal'#function

network_params[[4]][[6]][[1]] = .4#alpha
network_params[[4]][[6]][[2]] = .3#beta
network_params[[4]][[6]][[3]] = 1.6#balance
network_params[[4]][[6]][[4]] = 'pn'#assoc
network_params[[4]][[6]][[5]] = 'linear'#function





network_params[[5]][[5]][[1]] = 0#alpha
network_params[[5]][[5]][[2]] = .4#beta
network_params[[5]][[5]][[3]] = 1#balance
network_params[[5]][[5]][[4]] = 'p'#assoc
network_params[[5]][[5]][[5]] = 'linear'#function

network_params[[5]][[6]][[1]] = 1#alpha
network_params[[5]][[6]][[2]] = 0#beta
network_params[[5]][[6]][[3]] = .3#balance
network_params[[5]][[6]][[4]] = 'n'#assoc
network_params[[5]][[6]][[5]] = 'linear'#function




network_params[[6]][[6]][[1]] =.6 #alpha
network_params[[6]][[6]][[2]] = .3#beta
network_params[[6]][[6]][[3]] = 1#balance
network_params[[6]][[6]][[4]] = 'p'#assoc
network_params[[6]][[6]][[5]] = 'exponential'#function









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







total_runs = 31
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




fold_contents = append(list(NA), replicate( n=10, createFolds(1:length(adjacency_matrix[upper.tri(adjacency_matrix)]), k=3)) )

start_time = Sys.time()
for(run_no in 1:total_runs){ 
if(run_no == 1){
normalized_matrix = apply(adjacency_matrix, 2, function(x){(x- mean(x))/sd(x)})
svd_attempt = svd(normalized_matrix)
cos_sph2 = cosine(t(svd_attempt$u[,1:8])) #Need to check after the fact that 8 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "single")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
all_estimated_node_groupsU[[run_no]] = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)

#Now also do for V!
cos_sph2 = cosine(t(svd_attempt$v[,1:8])) #Need to check after the fact that 8 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "single")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
all_estimated_node_groupsV[[run_no]] = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)



eig_norm = eigen(normalized_matrix)
k_l =8
  sph_mat = Re(eig_norm$vectors)
  for(i in 1:nrow(sph_mat)){
    scaler = sqrt(sum(sph_mat[i,1:k_l]^2))
    sph_mat[i,1:k_l] = sph_mat[i,1:k_l]/scaler
  }  
  
  
  cos_sph = cosine(t(Re(sph_mat[,1:k_l])))
  hdist = hclust(as.dist(1- cos_sph), method = "single")
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
cos_sph2 = cosine(t(svd_attempt$u[,1:8])) #Need to check after the fact that 3 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "single")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
all_estimated_node_groupsU[[run_no]] = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)


#Now also do for V!
cos_sph2 = cosine(t(svd_attempt$v[,1:8])) #Need to check after the fact that 8 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "single")
#plot(hdist$height)
#choose the h that shows a big gap
heightmin = which.max(diff(hdist$height))
all_estimated_node_groupsV[[run_no]] = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)

#Eigen Version
eig_norm = eigen(Dineof$Xa)
k_l =8
  sph_mat = Re(eig_norm$vectors)
  for(i in 1:nrow(sph_mat)){
    scaler = sqrt(sum(sph_mat[i,1:k_l]^2))
    sph_mat[i,1:k_l] = sph_mat[i,1:k_l]/scaler
  }  
  
  cos_sph = cosine(t(Re(sph_mat[,1:k_l])))
  hdist = hclust(as.dist(1- cos_sph), method = "single")
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


save.image(file = "SimCommDet.RData")

