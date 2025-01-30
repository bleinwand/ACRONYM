load("MouseFitU.RData")



UFullLikelihood =  -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])


BlindreconstructedEstimates = vector(mode='list', length = 10)
for(blinder in 1:10){
  BlindreconstructedEstimates[[blinder]] =  0*adjacency_matrix
  for(k in (2+(3*(blinder-1))):((3*blinder)+1)){
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


U30Runtimes = unlist(lapply(2:31, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
U30RealisticTime = unlist(lapply(2:31, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))



load("MouseFitV.RData")




for(run in 1:total_runs){
  if(VRepeats[run]>0){
    reconstructedEstimates[[run]] = reconstructedEstimatesU[[run]] 
  } 
}


VFullLikelihood =  -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])



BlindreconstructedEstimates = vector(mode='list', length = 10)
for(blinder in 1:10){
  BlindreconstructedEstimates[[blinder]] =  0*adjacency_matrix
  for(k in (2+(3*(blinder-1))):((3*blinder)+1)){
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

V30Runtimes = unlist(lapply(2:31, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
V30RealisticTime = unlist(lapply(2:31, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))






load("MouseFitEig.RData")

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
  for(k in (2+(3*(blinder-1))):((3*blinder)+1)){
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

Eig30Runtimes = unlist(lapply(2:31, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
Eig30RealisticTime = unlist(lapply(2:31, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))




UFullLikelihood
VFullLikelihood
EigFullLikelihood
mean(UBlindLikelihoods)
mean(VBlindLikelihoods)
mean(EigBlindLikelihoods)


mean(unlist(lapply(1:length(all_estimated_node_groupsU), function(x){max(all_estimated_node_groupsU[[x]])})))
mean(unlist(lapply(1:length(all_estimated_node_groupsV), function(x){max(all_estimated_node_groupsV[[x]])})))
mean(unlist(lapply(1:length(all_estimated_node_groupsEig), function(x){max(all_estimated_node_groupsEig[[x]])})))

mean(unlist(lapply(2:31, function(x){adjustedRandIndex(all_estimated_node_groupsU[[x]], all_estimated_node_groupsU[[1]])})))
mean(unlist(lapply(2:31, function(x){adjustedRandIndex(all_estimated_node_groupsV[[x]], all_estimated_node_groupsV[[1]])})))
mean(unlist(lapply(2:31, function(x){adjustedRandIndex(all_estimated_node_groupsEig[[x]], all_estimated_node_groupsEig[[1]])})))

31- sum(VRepeats)
31- length(which(EigRepeats>0))


UFullyObservedEstTime/60
VFullyObservedEstTime/60
EigFullyObservedEstTime/60

UFullyObservedRealisticTime/60
VFullyObservedRealisticTime/60
EigFullyObservedRealisticTime/60

Utotal_runtime*60
Vtotal_runtime*60
Eigtotal_runtime*60

sum(U30Runtimes)/60
sum(V30Runtimes)/60
sum(Eig30Runtimes)/60

firstcol =c(mean(U30RealisticTime)/60,
mean(V30RealisticTime)/60,
mean(Eig30RealisticTime)/60)
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





library(randnet)
sparse_adj <- as(adjacency_matrix, "sparseMatrix") 

randnet_num_clus = BHMC.estimate(sparse_adj, K.max = 10)$K

randnet_clus =  reg.SP(sparse_adj,K=randnet_num_clus,lap=TRUE)
randnet_clus2 =  reg.SSP(sparse_adj,K=randnet_num_clus)

randnet_DCBM_est1 =  DCSBM.estimate(sparse_adj, randnet_clus$cluster)
randnet_DCBM_est2 = DCSBM.estimate(sparse_adj, randnet_clus2$cluster)
randnet_DCBM_est1$Phat  = Adjust_within_community_estimates(randnet_DCBM_est1$Phat, randnet_clus$cluster)
randnet_DCBM_est2$Phat  = Adjust_within_community_estimates(randnet_DCBM_est2$Phat, randnet_clus2$cluster)


range(randnet_DCBM_est1$Phat)
range(randnet_DCBM_est2$Phat)

valid_randnet_SP_9_entries = intersect(which(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)] <=1), which(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)] >=0))
valid_randnet_SSP_9_entries = intersect(which(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)] <=1), which(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)] >=0))
paste("Negative log-likelihood of randnet_DCBM_est1 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_9_entries]*log(as.vector(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_9_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_9_entries])*log(as.vector(1-randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_9_entries]))), 3))


paste("Negative log-likelihood of randnet_DCBM_est2 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_9_entries]*log(as.vector(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_9_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_9_entries])*log(as.vector(1-randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_9_entries]))), 3))


length(valid_randnet_SP_9_entries)
length(valid_randnet_SSP_9_entries)


library(nett)

bethe_hessian_select(adjacency_matrix, 10)

zh = spec_clust(adjacency_matrix, K=9)


#nett_DCBM_est = estim_dcsbm(sparse_adj, zh)
#library(randnet)


randnet_DCBM_est_using_nett_labels = DCSBM.estimate(sparse_adj, zh)
range(randnet_DCBM_est_using_nett_labels$Phat)
randnet_DCBM_est_using_nett_labels$Phat  = Adjust_within_community_estimates(randnet_DCBM_est_using_nett_labels$Phat, zh)
valid_nett_9_entries = intersect(which(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)] <=1), which(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)] >=0))
paste("Negative log-likelihood of DCBM estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_9_entries]*log(as.vector(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_9_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_9_entries])*log(as.vector(1-randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_9_entries]))), 3))


length(valid_nett_9_entries)











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

#check_PABM = PABM_estimate(adjacency_matrix, all_estimated_node_groupsV[[1]])


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

valid_PABM_entries = intersect(which(check_PABM[upper.tri(check_PABM)] <=1), which(check_PABM[upper.tri(check_PABM)] >=0))

paste("Negative log-likelihood of PABM estimate using our communities is",  round(-sum(log(as.vector(check_PABM[upper.tri(check_PABM)][valid_PABM_entries][which(adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries]==1)]))) -sum(log(as.vector(1-check_PABM[upper.tri(check_PABM)][valid_PABM_entries][which(adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries]==0)]))), 3))




check_PABM_truncated = check_PABM
check_PABM_truncated[which(check_PABM_truncated>1)] = .999
paste("Negative log-likelihood of truncated PABM estimate using our communities is",  round(-sum(log(as.vector(check_PABM_truncated[upper.tri(check_PABM_truncated)][which(adjacency_matrix[upper.tri(adjacency_matrix)]==1)]))) -sum(log(as.vector(1-check_PABM_truncated[upper.tri(check_PABM_truncated)][which(adjacency_matrix[upper.tri(adjacency_matrix)]==0)]))), 3))

length(intersect(which(check_PABM_truncated[upper.tri(check_PABM_truncated)] <=1), which(check_PABM_truncated[upper.tri(check_PABM_truncated)] >=0)))











randnet_true_comms =  DCSBM.estimate(sparse_adj, all_estimated_node_groupsV[[1]])
randnet_true_comms$Phat  = Adjust_within_community_estimates(randnet_true_comms$Phat, all_estimated_node_groupsV[[1]])

valid_randnet_true_comms = intersect(which(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)] <=1), which(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)] >=0))
#paste("Negative log-likelihood of DCBM with true community estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)]))), 3))

paste("Negative log-likelihood of DCBM with our community estimate is",  round(-sum(log(as.vector(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)][valid_randnet_true_comms][which(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_true_comms]==1)]))) -sum(log(as.vector(1-randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)][valid_randnet_true_comms][which(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_true_comms]==0)]))), 3))

length(valid_randnet_true_comms)





library(RColorBrewer)
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

par(mfrow=c(4, 3))

par(mar=c(1.2,0 ,1.1, 1.2))
par(pty="s")

plot(x = MRsoma$x[-no_contacts], y = MRsoma$y[-no_contacts], col=c25[all_estimated_node_groupsV[[1]]], pch=20, cex=.9, ylab=NA, xaxt= 'n', yaxt='n')
title(ylab="y position", line=.5, cex.lab=1.2)
title(xlab = "x position", line=0, cex.lab=1.2)

plot(x = MRsoma$x[-no_contacts], y = MRsoma$z[-no_contacts], col=c25[all_estimated_node_groupsV[[1]]], pch=20, cex=.9, ylab=NA, xaxt= 'n', yaxt='n', main ="ACRONYM")
title(ylab="z position", line=.5, cex.lab=1.2)
title(xlab = "x position", line=0, cex.lab=1.2)

plot(x = MRsoma$y[-no_contacts], y = MRsoma$z[-no_contacts], col=c25[all_estimated_node_groupsV[[1]]], pch=20, cex=.9,ylab=NA, xaxt= 'n', yaxt='n')
title(ylab="z position", line=.5, cex.lab=1.2)
title(xlab = "y position", line=0, cex.lab=1.2)



plot(x = MRsoma$x[-no_contacts], y = MRsoma$y[-no_contacts], col=c25[as.numeric(randnet_clus$cluster)], pch=20, cex=.9, ylab=NA, xaxt= 'n', yaxt='n')
title(ylab="y position", line=.5, cex.lab=1.2)
title(xlab = "x position", line=0, cex.lab=1.2)

plot(x = MRsoma$x[-no_contacts], y = MRsoma$z[-no_contacts], col=c25[as.numeric(randnet_clus$cluster)], pch=20, cex=.9, ylab=NA, xaxt= 'n', yaxt='n', main ="randnet Spectral Clustering")
title(ylab="z position", line=.5, cex.lab=1.2)
title(xlab = "x position", line=0, cex.lab=1.2)

plot(x = MRsoma$y[-no_contacts], y = MRsoma$z[-no_contacts], col=c25[as.numeric(randnet_clus$cluster)], pch=20, cex=.9,ylab=NA, xaxt= 'n', yaxt='n')
title(ylab="z position", line=.5, cex.lab=1.2)
title(xlab = "y position", line=0, cex.lab=1.2)


plot(x = MRsoma$x[-no_contacts], y = MRsoma$y[-no_contacts], col=c25[as.numeric(randnet_clus2$cluster)], pch=20, cex=.9, ylab=NA, xaxt= 'n', yaxt='n')
title(ylab="y position", line=.5, cex.lab=1.2)
title(xlab = "x position", line=0, cex.lab=1.2)

plot(x = MRsoma$x[-no_contacts], y = MRsoma$z[-no_contacts], col=c25[as.numeric(randnet_clus2$cluster)], pch=20, cex=.9, ylab=NA, xaxt= 'n', yaxt='n', main ="randnet Spherical Spectral Clustering")
title(ylab="z position", line=.5, cex.lab=1.2)
title(xlab = "x position", line=0, cex.lab=1.2)

plot(x = MRsoma$y[-no_contacts], y = MRsoma$z[-no_contacts], col=c25[as.numeric(randnet_clus2$cluster)], pch=20, cex=.9,ylab=NA, xaxt= 'n', yaxt='n')
title(ylab="z position", line=.5, cex.lab=1.2)
title(xlab = "y position", line=0, cex.lab=1.2)


plot(x = MRsoma$x[-no_contacts], y = MRsoma$y[-no_contacts], col=c25[as.numeric(zh)], pch=20, cex=.9, ylab=NA, xaxt= 'n', yaxt='n')
title(ylab="y position", line=.5, cex.lab=1.2)
title(xlab = "x position", line=0, cex.lab=1.2)

plot(x = MRsoma$x[-no_contacts], y = MRsoma$z[-no_contacts], col=c25[as.numeric(zh)], pch=20, cex=.9, ylab=NA, xaxt= 'n', yaxt='n', main ="nett Spectral Clustering")
title(ylab="z position", line=.5, cex.lab=1.2)
title(xlab = "x position", line=0, cex.lab=1.2)

plot(x = MRsoma$y[-no_contacts], y = MRsoma$z[-no_contacts], col=c25[as.numeric(zh)], pch=20, cex=.9,ylab=NA, xaxt= 'n', yaxt='n')
title(ylab="z position", line=.5, cex.lab=1.2)
title(xlab = "y position", line=0, cex.lab=1.2)






within_group_degree = rep(0, nrow(MRgraph2))
for(cell in 1:length(within_group_degree)){
  within_group_degree[cell] = sum(MRgraph2[cell, which(all_estimated_node_groupsV[[1]] == all_estimated_node_groupsV[[1]][cell])] )
}


reordering = cbind(1:nrow(MRgraph2), all_estimated_node_groupsV[[1]], within_group_degree)

reordering = reordering[order(reordering[,2], reordering[,3]),]

dev.off()
image(MRgraph2[reordering[,1], reordering[,1]], col=c("white", "black"))

check_PABM_for_plot = check_PABM
check_PABM_for_plot[which(check_PABM_for_plot>1)] = NA

DCBM_for_plot = randnet_DCBM_est_using_nett_labels$Phat
DCBM_for_plot[which(DCBM_for_plot>1)] = NA


image(check_PABM_for_plot[reordering[,1], reordering[,1]], col = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(20))
image(DCBM_for_plot[reordering[,1], reordering[,1]], col = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(20))
image(reconstructedEstimatesV[[1]][reordering[,1], reordering[,1]], col = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(20))
image(VBlindreconstructedEstimates[[1]][reordering[,1], reordering[,1]], col = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(20))

fullacronym = reconstructedEstimatesV[[1]]
diag(fullacronym) = 0
image(MRgraph2[reordering[which(reordering[,2]==4),1][192:241], reordering[which(reordering[,2]==4),1][192:241]], col=c("white", "black"))
filled.contour(randnet_DCBM_est_using_nett_labels$Phat[reordering[which(reordering[,2]==4),1][192:241], reordering[which(reordering[,2]==4),1][192:241]], color.palette = colorRampPalette(c("purple","blue", "green", "yellow", "orange", "red")), levels = seq(0, 7, by = .5), plot.axes = 0)
filled.contour(fullacronym[reordering[which(reordering[,2]==4),1][192:241], reordering[which(reordering[,2]==4),1][192:241]], color.palette = colorRampPalette(c("purple","blue", "green", "yellow", "orange", "red")), levels = seq(0, 7, by = .5), plot.axes = 0)
filled.contour(VBlindreconstructedEstimates[[1]][reordering[which(reordering[,2]==4),1][192:241], reordering[which(reordering[,2]==4),1][192:241]], color.palette = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red")), levels = seq(0, 7, by =.5 ), plot.axes = 0 )
filled.contour(check_PABM[reordering[which(reordering[,2]==4),1][192:241], reordering[which(reordering[,2]==4),1][192:241]], color.palette = colorRampPalette(c("purple","blue", "green", "yellow", "orange", "red")), levels = seq(0, 7, by = .5), plot.axes = 0)






library(igraph)
gamma <- seq(0.25,1,0.01)
nc <- vector("numeric",length(gamma))
for (i in 1:length(gamma)){
  gc <- cluster_leiden(graph_from_adjacency_matrix(MRgraph2, mode = "undirected"), objective_function = "modularity",
                       n_iterations = 3, resolution_parameter = gamma[i])
  nc[i] <- length(gc)
}
plot(gamma,nc,xlab="gamma",ylab="# Communities",main="Mouse Retina")

check_leiden <- cluster_leiden(graph_from_adjacency_matrix(MRgraph2, mode = "undirected"), objective_function = "modularity",
                               n_iterations = 50, resolution_parameter = .75)




tocomp = cbind(all_estimated_node_groupsU[[1]], all_estimated_node_groupsV[[1]], all_estimated_node_groupsEig[[1]], SSCclust$cluster,  zh, randnet_clus$cluster, randnet_clus2$cluster, check_leiden$membership, 
               factor(MRtype[MRcells[-no_contacts, 2], 5]),factor(MRtype[MRcells[-no_contacts, 2], 5] )
)

tocomp[which(is.na(tocomp[,9])),10] = 6

ARI_mat = matrix(0, nrow=10, ncol=10)
#NMI_mat = matrix(0, nrow=6, ncol=6)
for(i in 1:10){
  for(j in 1:10){
    
    ARI_mat[i,j] = adjustedRandIndex(tocomp[,i], tocomp[,j])
    #NMI_mat[i,j] = NMI(tocomp[,i], tocomp[,j])
  }
}

ARI_mat
