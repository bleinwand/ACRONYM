load("CongressFitU.RData")



UFullLikelihood =  -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])


BlindreconstructedEstimates = vector(mode='list', length = 10)
for(blinder in 1:10){
  BlindreconstructedEstimates[[blinder]] =  0*adjacency_matrix
  for(k in (2+(10*(blinder-1))):((10*blinder)+1)){
    BlindreconstructedEstimates[[blinder]][upper.tri(BlindreconstructedEstimates[[blinder]])][fold_contents[[k]]] = reconstructedEstimates[[k]][upper.tri(reconstructedEstimates[[k]])][fold_contents[[k]]]
  }
  BlindreconstructedEstimates[[blinder]] = BlindreconstructedEstimates[[blinder]] + t(BlindreconstructedEstimates[[blinder]])
}



#BlindMSEs = rep(0, 10)
BlindLikelihoods = rep(0, 10)

for(i in 1:10){
  #BlindMSEs[i] = mean((BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])] - no_error_probs_matrix[lower.tri(no_error_probs_matrix)])^2)
  BlindLikelihoods[i] = -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])
}



reconstructedEstimatesU = reconstructedEstimates
UBlindreconstructedEstimates = BlindreconstructedEstimates
#UBlindMSEs = BlindMSEs
UBlindLikelihoods = BlindLikelihoods 

Utotal_runtime = total_runtime
UFullyObservedEstTime = sum(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))
UFullyObservedRealisticTime = max(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))

U100Runtimes = unlist(lapply(2:101, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
U100RealisticTime = unlist(lapply(2:101, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))



load("CongressFitV.RData")




for(run in 1:total_runs){
  if(VRepeats[run]>0){
    reconstructedEstimates[[run]] = reconstructedEstimatesU[[run]] 
  } 
}


VFullLikelihood =  -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])



BlindreconstructedEstimates = vector(mode='list', length = 10)
for(blinder in 1:10){
  BlindreconstructedEstimates[[blinder]] =  0*adjacency_matrix
  for(k in (2+(10*(blinder-1))):((10*blinder)+1)){
    BlindreconstructedEstimates[[blinder]][upper.tri(BlindreconstructedEstimates[[blinder]])][fold_contents[[k]]] = reconstructedEstimates[[k]][upper.tri(reconstructedEstimates[[k]])][fold_contents[[k]]]
  }
  BlindreconstructedEstimates[[blinder]] = BlindreconstructedEstimates[[blinder]] + t(BlindreconstructedEstimates[[blinder]])
}





#BlindMSEs = rep(0, 10)
BlindLikelihoods = rep(0, 10)

for(i in 1:10){
  #BlindMSEs[i] = mean((BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])] - no_error_probs_matrix[lower.tri(no_error_probs_matrix)])^2)
  BlindLikelihoods[i] = -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])
}



reconstructedEstimatesV = reconstructedEstimates
VBlindreconstructedEstimates = BlindreconstructedEstimates
#VBlindMSEs = BlindMSEs
VBlindLikelihoods = BlindLikelihoods 

Vtotal_runtime = total_runtime
VFullyObservedEstTime = sum(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))
VFullyObservedRealisticTime = max(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))


V100Runtimes = unlist(lapply(2:101, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
V100RealisticTime = unlist(lapply(2:101, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))





load("CongressFitEig.RData")

for(run in 1:total_runs){
  if(EigRepeats[run]==1){
    reconstructedEstimates[[run]] = reconstructedEstimatesU[[run]] 
  } 
  if(EigRepeats[run]==2){
    reconstructedEstimates[[run]] = reconstructedEstimatesV[[run]] 
  } 
  
  
}





EigFullLikelihood =  -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])



BlindreconstructedEstimates = vector(mode='list', length = 10)
for(blinder in 1:10){
  BlindreconstructedEstimates[[blinder]] =  0*adjacency_matrix
  for(k in (2+(10*(blinder-1))):((10*blinder)+1)){
    BlindreconstructedEstimates[[blinder]][upper.tri(BlindreconstructedEstimates[[blinder]])][fold_contents[[k]]] = reconstructedEstimates[[k]][upper.tri(reconstructedEstimates[[k]])][fold_contents[[k]]]
  }
  BlindreconstructedEstimates[[blinder]] = BlindreconstructedEstimates[[blinder]] + t(BlindreconstructedEstimates[[blinder]])
}





#BlindMSEs = rep(0, 10)
BlindLikelihoods = rep(0, 10)

for(i in 1:10){
  #BlindMSEs[i] = mean((BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])] - no_error_probs_matrix[lower.tri(no_error_probs_matrix)])^2)
  BlindLikelihoods[i] = -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])
}



reconstructedEstimatesEig = reconstructedEstimates
EigBlindreconstructedEstimates = BlindreconstructedEstimates
#EigBlindMSEs = BlindMSEs
EigBlindLikelihoods = BlindLikelihoods 


Eigtotal_runtime = total_runtime
EigFullyObservedEstTime = sum(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))
EigFullyObservedRealisticTime = max(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))

Eig100Runtimes = unlist(lapply(2:101, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
Eig100RealisticTime = unlist(lapply(2:101, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))





UFullLikelihood
VFullLikelihood
EigFullLikelihood
mean(UBlindLikelihoods)
mean(VBlindLikelihoods)
mean(EigBlindLikelihoods)


mean(unlist(lapply(1:length(all_estimated_node_groupsU), function(x){max(all_estimated_node_groupsU[[x]])})))
mean(unlist(lapply(1:length(all_estimated_node_groupsV), function(x){max(all_estimated_node_groupsV[[x]])})))
mean(unlist(lapply(1:length(all_estimated_node_groupsEig), function(x){max(all_estimated_node_groupsEig[[x]])})))


mean(unlist(lapply(2:101, function(x){adjustedRandIndex(all_estimated_node_groupsU[[x]], all_estimated_node_groupsU[[1]])})))
 mean(unlist(lapply(2:101, function(x){adjustedRandIndex(all_estimated_node_groupsV[[x]], all_estimated_node_groupsV[[1]])})))
 mean(unlist(lapply(2:101, function(x){adjustedRandIndex(all_estimated_node_groupsEig[[x]], all_estimated_node_groupsEig[[1]])})))


################################################################################



length(which(unlist(reconstructedEstimatesV) < 0))
length(which(unlist(reconstructedEstimatesV) > 1))
length(reconstructedEstimatesV[[1]][upper.tri(reconstructedEstimatesV[[1]])])








Adjust_within_community_estimates  = function(Estimates, g){
  K <- length(unique(g))
  #B <- matrix(0, K, K)
  P.hat <- Estimates
  for (i in 1:K) {
    N.i <- which(g == i)
    A_prime = Estimates[N.i, N.i]
    B  <- sum(A_prime, na.rm = TRUE) + 
      0.001
    Psi <- rowSums(A_prime, na.rm = TRUE)
    Psi <- Psi/sum(Psi)
    
    P.hat[N.i, N.i] = B*Psi%*%t(Psi) #/ +  # (1-outer(Psi, Psi, "+"))#((1-Psi)%*%t(1-Psi2))
    toADD = t(outer(Psi, 1-Psi, "/"))*diag(P.hat[N.i, N.i])
    P.hat[N.i, N.i] = P.hat[N.i, N.i] + ((toADD + t(toADD))/2) 
    diag(P.hat[N.i, N.i]) = 0
  }
  return(P.hat)
  
}



library(nett)

bethe_hessian_select(adjacency_matrix, 10)

zh = spec_clust(adjacency_matrix, K=7)

sparse_adj <- as(adjacency_matrix, "sparseMatrix") 

nett_DCBM_est = estim_dcsbm(sparse_adj, zh)
library(randnet)


randnet_DCBM_est_using_nett_labels = DCSBM.estimate(sparse_adj, zh)
range(randnet_DCBM_est_using_nett_labels$Phat)
randnet_DCBM_est_using_nett_labels$Phat  = Adjust_within_community_estimates(randnet_DCBM_est_using_nett_labels$Phat, zh)
valid_nett_7_entries = intersect(which(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)] <=1), which(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)] >=0))
paste("Negative log-likelihood of DCBM estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_7_entries]*log(as.vector(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_7_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_7_entries])*log(as.vector(1-randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_7_entries]))), 3))




randnet_num_clus = BHMC.estimate(sparse_adj, K.max = 10)$K


randnet_clus =  reg.SP(sparse_adj,K=7,lap=TRUE)
randnet_clus2 =  reg.SSP(sparse_adj,K=7)

randnet_DCBM_est1 =  DCSBM.estimate(sparse_adj, randnet_clus$cluster)
randnet_DCBM_est2 = DCSBM.estimate(sparse_adj, randnet_clus2$cluster)
randnet_DCBM_est1$Phat  = Adjust_within_community_estimates(randnet_DCBM_est1$Phat, randnet_clus$cluster)
randnet_DCBM_est2$Phat  = Adjust_within_community_estimates(randnet_DCBM_est2$Phat, randnet_clus2$cluster)


range(randnet_DCBM_est1$Phat)
range(randnet_DCBM_est2$Phat)

valid_randnet_SP_7_entries = intersect(which(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)] <=1), which(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)] >=0))
valid_randnet_SSP_7_entries = intersect(which(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)] <=1), which(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)] >=0))
paste("Negative log-likelihood of randnet_DCBM_est1 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_7_entries]*log(as.vector(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_7_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_7_entries])*log(as.vector(1-randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_7_entries]))), 3))


paste("Negative log-likelihood of randnet_DCBM_est2 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_7_entries]*log(as.vector(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_7_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_7_entries])*log(as.vector(1-randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_7_entries]))), 3))

length(valid_nett_7_entries)
length(valid_randnet_SP_7_entries)
length(valid_randnet_SSP_7_entries)






randnet_true_comms =  DCSBM.estimate(sparse_adj, all_estimated_node_groupsV[[1]])
randnet_true_comms$Phat  = Adjust_within_community_estimates(randnet_true_comms$Phat, all_estimated_node_groupsV[[1]])

valid_randnet_true_comms = intersect(which(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)] <=1), which(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)] >=0))
#paste("Negative log-likelihood of DCBM with true community estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)]))), 3))

paste("Negative log-likelihood of DCBM with our community estimate is",  round(-sum(log(as.vector(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)][valid_randnet_true_comms][which(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_true_comms]==1)]))) -sum(log(as.vector(1-randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)][valid_randnet_true_comms][which(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_true_comms]==0)]))), 3))

length(valid_randnet_true_comms)





PABM_estimate = function(A, g)
{
  K <- length(unique(g))
  #B <- matrix(0, K, K)
  P.hat <- 0*A
  for (i in 1:K) {
    for (j in 1:K) {
      N.i <- which(g == i)
      N.j <- which(g == j)
      A_prime = A[N.i, N.j]
      B  <- sum(A_prime, na.rm = TRUE) + 
        0.001
      Psi <- rowSums(A_prime, na.rm = TRUE)
      Psi <- Psi/sum(Psi)
      Psi2 <- colSums(A_prime, na.rm = TRUE)
      Psi2 <- Psi2/sum(Psi2)
      P.hat[which(g == i), which(g == j)] = B*Psi%*%t(Psi2)
      
      if(i == j){
        P.hat[N.i, N.j] = B*Psi%*%t(Psi2) #/ +  # (1-outer(Psi, Psi, "+"))#((1-Psi)%*%t(1-Psi2))
        toADD = t(outer(Psi, 1-Psi, "/"))*diag(P.hat[N.i, N.j])
        P.hat[N.i, N.j] = P.hat[N.i, N.j] + ((toADD + t(toADD))/2) 
        diag(P.hat[N.i, N.j]) = 0
      }
    }
    
  }
  
  
  return(P.hat)
}

#check_PABM = PABM_estimate(total_list, all_estimated_node_groupsV[[1]])


#valid_PABM_entries = intersect(which(check_PABM[upper.tri(check_PABM)] <=1), which(check_PABM[upper.tri(check_PABM)] >=0))

library(T4cluster)
eigforSSC = eigen(adjacency_matrix)
numComms = max(all_estimated_node_groupsV[[1]])
numVecs = numComms*(numComms-1)/2
vecsToKeep = c(1:numVecs,  (ncol(eigforSSC$vectors)-numVecs+1):ncol(eigforSSC$vectors))
newSSC = eigforSSC$vectors[,vecsToKeep]
newSSC = sqrt(nrow(newSSC))*newSSC
SSCStartTime = Sys.time()
SSCclust = SSC(newSSC, k=numComms)
SSCRuntime = Sys.time() - SSCStartTime
check_PABM =  PABM_estimate(adjacency_matrix, SSCclust$cluster)
range(check_PABM[upper.tri(check_PABM)])


paste("Negative log-likelihood of PABM estimate using SSC communities is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(check_PABM[upper.tri(check_PABM)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-check_PABM[upper.tri(check_PABM)]))), 3))

#valid_PABM_entries = intersect(which(check_PABM[upper.tri(check_PABM)] <=1), which(check_PABM[upper.tri(check_PABM)] >=0))

valid_PABM_entries = intersect(which(check_PABM[upper.tri(check_PABM)] <1), which(check_PABM[upper.tri(check_PABM)] >0))


paste("Negative log-likelihood of PABM estimate using SSC communities is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries]*log(as.vector(check_PABM[upper.tri(check_PABM)][valid_PABM_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries])*log(as.vector(1-check_PABM[upper.tri(check_PABM)][valid_PABM_entries]))), 3))



check_PABM_truncated = check_PABM
check_PABM_truncated[which(check_PABM_truncated>1)] = .999
paste("Negative log-likelihood of truncated PABM estimate using our communities is",  round(-sum(log(as.vector(check_PABM_truncated[upper.tri(check_PABM_truncated)][which(adjacency_matrix[upper.tri(adjacency_matrix)]==1)]))) -sum(log(as.vector(1-check_PABM_truncated[upper.tri(check_PABM_truncated)][which(adjacency_matrix[upper.tri(adjacency_matrix)]==0)]))), 3))

length(intersect(which(check_PABM_truncated[upper.tri(check_PABM_truncated)] <=1), which(check_PABM_truncated[upper.tri(check_PABM_truncated)] >=0)))






Utotal_runtime*60
Vtotal_runtime*60
Eigtotal_runtime*60

UFullyObservedEstTime/60
VFullyObservedEstTime/60
EigFullyObservedEstTime/60

101- sum(VRepeats)
101- length(which(EigRepeats>0))


UFullyObservedRealisticTime/60
VFullyObservedRealisticTime/60
EigFullyObservedRealisticTime/60



sum(U100Runtimes)/60
mean(U100RealisticTime)/60
sum(V100Runtimes)/60
mean(V100RealisticTime)/60
sum(Eig100Runtimes)/60
mean(Eig100RealisticTime)/60









library(igraph)




gamma <- seq(0.25,2,0.025)
nc <- vector("numeric",length(gamma))
for (i in 1:length(gamma)){
  gc <- cluster_leiden(graph_from_adjacency_matrix(total_list, mode = "undirected"), objective_function = "modularity",
                       n_iterations = 3, resolution_parameter = gamma[i])
  nc[i] <- length(gc)
}
plot(gamma,nc,xlab="gamma",ylab="# Communities",main="Congress")
check_leiden <- cluster_leiden(graph_from_adjacency_matrix(total_list, mode = "undirected"), objective_function = "modularity",
                               n_iterations = 3, resolution_parameter = .75)
plot(check_leiden,graph_from_adjacency_matrix(total_list, mode = "undirected"))


library(mclust)

tocomp = cbind(all_estimated_node_groupsU[[1]], all_estimated_node_groupsV[[1]], 
               all_estimated_node_groupsEig[[1]], SSCclust$cluster,  zh, 
               randnet_clus$cluster, randnet_clus2$cluster,  check_leiden$membership
)



ARI_mat = matrix(0, nrow=ncol(tocomp), ncol=ncol(tocomp))
#NMI_mat = matrix(0, nrow=6, ncol=6)
for(i in 1:ncol(tocomp)){
  for(j in 1:ncol(tocomp)){
    
    ARI_mat[i,j] = adjustedRandIndex(tocomp[,i], tocomp[,j])
    #NMI_mat[i,j] = NMI(tocomp[,i], tocomp[,j])
  }
}

ARI_mat





normalized_matrix = apply(adjacency_matrix, 2, function(x){(x- mean(x))/sd(x)})
#svd_attempt = svd(normalized_matrix)
eig_norm = eigen(normalized_matrix)
k_l =2
sph_mat = Re(eig_norm$vectors)
for(i in 1:nrow(sph_mat)){
  scaler = sqrt(sum(sph_mat[i,1:k_l]^2))
  sph_mat[i,1:k_l] = sph_mat[i,1:k_l]/scaler
}  
plot(sph_mat, col = all_estimated_node_groupsU[[1]])
plot(sph_mat, col = all_estimated_node_groupsV[[1]])
plot(sph_mat, col = all_estimated_node_groupsEig[[1]])


k_l =4
sph_mat = Re(eig_norm$vectors)
for(i in 1:nrow(sph_mat)){
  scaler = sqrt(sum(sph_mat[i,1:k_l]^2))
  sph_mat[i,1:k_l] = sph_mat[i,1:k_l]/scaler
}  

dim(sph_mat)
plot(sph_mat, col = all_estimated_node_groupsU[[1]])
plot(sph_mat, col = all_estimated_node_groupsV[[1]])
plot(sph_mat, col = all_estimated_node_groupsEig[[1]])


svd_attempt = svd(normalized_matrix)



plot(svd_attempt$u, col =all_estimated_node_groupsU[[1]])
plot(svd_attempt$u, col =all_estimated_node_groupsV[[1]])
plot(svd_attempt$u, col =all_estimated_node_groupsEig[[1]])

plot(svd_attempt$v, col =all_estimated_node_groupsU[[1]])
plot(svd_attempt$v, col =all_estimated_node_groupsV[[1]])
plot(svd_attempt$v, col =all_estimated_node_groupsEig[[1]])




#########################################
###Figure 1
##########################################
within_group_degree = rep(0, 475)
for(congmember in 1:length(within_group_degree)){
  within_group_degree[congmember] = sum(adjacency_matrix[congmember, which(all_estimated_node_groupsV[[1]] == all_estimated_node_groupsV[[1]][congmember])] )
}

reordering = cbind(1:475, all_estimated_node_groupsV[[1]], within_group_degree)

reordering = reordering[order(reordering[,2], reordering[,3]),]




image(adjacency_matrix[reordering[,1], reordering[,1]], col=c("white", "black"))
image(adjacency_matrix, col=c("white", "black"))


filled.contour(VBlindreconstructedEstimates[[1]][reordering[,1], reordering[,1]], color.palette = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05))
filled.contour(VBlindreconstructedEstimates[[1]], color.palette = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05))

filled.contour(check_PABM[reordering[,1], reordering[,1]], color.palette = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05))
filled.contour(check_PABM, color.palette = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05))

filled.contour(randnet_DCBM_est_using_nett_labels$Phat[reordering[,1], reordering[,1]], color.palette = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05))
filled.contour(randnet_DCBM_est_using_nett_labels$Phat, color.palette = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05))

