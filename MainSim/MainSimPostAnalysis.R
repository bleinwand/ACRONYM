
################################
####REDOING FIGURE 4 - 5
################################
#First part will be done using the output from our main estimate.
#However, since we didn't want to re-estimate things when the same communities are recovered,
#later, we will use the output from SimFitU as all 3 methods gave the same communities

load("SimFitV.RData")

normalized_matrix = apply(adjacency_matrix, 2, function(x){(x- mean(x))/sd(x)})
svd_attempt = svd(normalized_matrix)
cos_sph2 = cosine(t(svd_attempt$v[,1:8])) #Need to check after the fact that 8 is correct.
hdist = hclust(as.dist(1- cos_sph2), method = "single")
heightmin = which.max(diff(hdist$height))

fullyobservedcomms = cutree(hdist, h = sum(hdist$height[heightmin:(heightmin+1)] )/2)


plot(svd_attempt$d)
plot(hdist$height)
abline(h = sum(hdist$height[heightmin:(heightmin+1)] )/2, col='red')

plot(svd_attempt$v, col=fullyobservedcomms)
plot(svd_attempt$v[,3:4], col=fullyobservedcomms)

library(viridis)

plot(svd_attempt$v, col=viridis(100)[as.numeric(cut(node_sociabilities,breaks = 100))])

plot(svd_attempt$v[,3:4], col=viridis(100)[as.numeric(cut(node_sociabilities,breaks = 100))])





plot(1:200, node_sociabilities[1:200], col=1, pch=1, cex = .75)
points(1:200, node_sociabilities[201:400], col=2, pch=2, cex = .75)
points(1:200, node_sociabilities[401:600], col=3, pch=3, cex = .75)
points(1:200, node_sociabilities[601:800], col=4, pch=4, cex = .75)
points(1:200, node_sociabilities[801:1000], col=5, pch=5, cex = .75)
points(1:200, node_sociabilities[1001:1200], col=6, pch=6, cex = .75)


legend("topleft",  title="Community",
       c("1",  "2", "3", "4", "5", "6"), fill= 1:6, cex=0.8)





rm(list = ls())



################################
####REDOING FIGURE 5
################################

load("SimFitU.RData")



par(mfrow = c(6, 6))
par(mar = c(1, 1, 1, 1))
num_groups = 6 
for(groupno in num_groups:1){
  # Identify rows with the relevant group for plotting
local_comm_pairs = which(apply(all_comm_pairsU[[1]], 1, function(row) any(row == groupno)))


for(k in local_comm_pairs){
i = all_comm_pairsU[[1]][k, 1]
j = all_comm_pairsU[[1]][k, 2]
if(groupno == i){
      if(network_params[[i]][[j]][[5]] %in% c("exponential", "negexponential" )){
        relcol = ifelse(network_params[[i]][[j]][[3]]<=.5, "darkgreen", ifelse(network_params[[i]][[j]][[3]]<1,'green',  ifelse(network_params[[i]][[j]][[3]]>=2, 'red', ifelse(network_params[[i]][[j]][[3]]>1, 'lightpink', 4))))
      }else{
        
        relcol = ifelse(network_params[[i]][[j]][[3]]<=.5, "red", ifelse(network_params[[i]][[j]][[3]]<1,'lightpink',  ifelse(network_params[[i]][[j]][[3]]>=2, 'darkgreen', ifelse(network_params[[i]][[j]][[3]]>1, 'green', 4))))
      }
      if(network_params[[i]][[j]][[4]] %in% c('p', 'pn')){
        if(mean(adjacency_matrix[which(node_groups==i), which(node_groups==j)])<.5){
        plot(parallelOutput2[[1]][[k]][[5]], xlim= c(1,200), ylim = c(0,1), ylab = NA, xlab=NA, xaxt ='n', yaxt ='n', axes =F)
        }else{
          plot(1-parallelOutput2[[1]][[k]][[5]], xlim= c(1,200), ylim = c(0,1), ylab = NA, xlab=NA, xaxt ='n', yaxt ='n', axes =F)
        }
      }else{
        
        if(mean(adjacency_matrix[which(node_groups==i), which(node_groups==j)])<.5){
        plot(1-parallelOutput2[[1]][[k]][[5]], xlim= c(1,200), ylim = c(0,1), ylab = NA, xlab=NA, xaxt ='n', yaxt ='n', axes =F)
        }else{
          plot(parallelOutput2[[1]][[k]][[5]], xlim= c(1,200), ylim = c(0,1), ylab = NA, xlab=NA, xaxt ='n', yaxt ='n', axes =F)
        }
        
      }
      points(1:200, node_sociabilities[(((i-1)*200)+1):(((i-1)*200)+200)], col=relcol)
      }else{
        if(network_params[[i]][[j]][[5]] %in% c("exponential", "negexponential" )){
          relcol = ifelse(network_params[[i]][[j]][[3]]<=.5, "red", ifelse(network_params[[i]][[j]][[3]]<1,'lightpink',  ifelse(network_params[[i]][[j]][[3]]>=2, 'darkgreen', ifelse(network_params[[i]][[j]][[3]]>1, 'green', 4))))}else{
            relcol = ifelse(network_params[[i]][[j]][[3]]<=.5, "darkgreen", ifelse(network_params[[i]][[j]][[3]]<1,'green',  ifelse(network_params[[i]][[j]][[3]]>=2, 'red', ifelse(network_params[[i]][[j]][[3]]>1, 'lightpink', 4))))
          }
        
        
        if(network_params[[i]][[j]][[4]] %in% c('p', 'np')){
          if(mean(adjacency_matrix[which(node_groups==i), which(node_groups==j)])<.5){
          plot(parallelOutput2[[1]][[k]][[6]], xlim= c(1,200), ylim = c(0,1), ylab = NA, xlab=NA, xaxt ='n', yaxt ='n', axes =F)
          }else{
            
            plot(1-parallelOutput2[[1]][[k]][[6]], xlim= c(1,200), ylim = c(0,1), ylab = NA, xlab=NA, xaxt ='n', yaxt ='n', axes =F)
            
          }
            
            
            }else{
              if(mean(adjacency_matrix[which(node_groups==i), which(node_groups==j)])<.5){
                
            plot(1-parallelOutput2[[1]][[k]][[6]], xlim= c(1,200), ylim = c(0,1), ylab = NA, xlab=NA, xaxt ='n', yaxt ='n', axes =F)
              }else{
                plot(parallelOutput2[[1]][[k]][[6]], xlim= c(1,200), ylim = c(0,1), ylab = NA, xlab=NA, xaxt ='n', yaxt ='n', axes =F)
                
                }
              
              
            }
        
        points(1:200, node_sociabilities[(((j-1)*200)+1):(((j-1)*200)+200)], col=relcol)
        
      }
  }

}












################################
####REDOING FIGURE 6
################################


par(mfrow = c(1, 1))



#filled.contour(probs_matrix, color.palette = colorRampPalette(c("purple","blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05))


image(no_error_probs_matrix, col = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(20))
image(probs_matrix, col = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(20))
image(reconstructedEstimates[[1]], col = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(20))


image(adjacency_matrix, col = c("white", "black"))










################################
####REDOING FIGURE 4 and table 1 - calculating other estimatees
################################






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
bethe_hessian_select(adjacency_matrix, 24)
zh = spec_clust(adjacency_matrix, K=5)
#compute_mutual_info(node_groups, zh)


sparse_adj <- as(adjacency_matrix, "sparseMatrix") 

nett_DCBM_est = estim_dcsbm(sparse_adj, zh)



library(randnet)


randnet_DCBM_est_using_nett_labels = DCSBM.estimate(sparse_adj, zh)
range(randnet_DCBM_est_using_nett_labels$Phat)
randnet_DCBM_est_using_nett_labels$Phat  = Adjust_within_community_estimates(randnet_DCBM_est_using_nett_labels$Phat, zh)

valid_nett_5_entries = intersect(which(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)] <=1), which(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)] >=0))

paste("Negative log-likelihood of DCBM estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_5_entries]*log(as.vector(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_5_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_5_entries])*log(as.vector(1-randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_5_entries]))), 3))

round(mean((randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_5_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_nett_5_entries])^2), 5)


randnet_num_clus = BHMC.estimate(sparse_adj, K.max = 24)$K
randnet_clus =  reg.SP(sparse_adj,K=5,lap=TRUE)
randnet_clus2 =  reg.SSP(sparse_adj,K=5)

randnet_DCBM_est1 =  DCSBM.estimate(sparse_adj, randnet_clus$cluster)
randnet_DCBM_est2 = DCSBM.estimate(sparse_adj, randnet_clus2$cluster)
randnet_DCBM_est1$Phat  = Adjust_within_community_estimates(randnet_DCBM_est1$Phat, randnet_clus$cluster)
randnet_DCBM_est2$Phat  = Adjust_within_community_estimates(randnet_DCBM_est2$Phat, randnet_clus2$cluster)


range(randnet_DCBM_est1$Phat)
range(randnet_DCBM_est2$Phat)

valid_randnet_SP_5_entries = intersect(which(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)] <=1), which(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)] >=0))
valid_randnet_SSP_5_entries = intersect(which(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)] <=1), which(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)] >=0))
paste("Negative log-likelihood of randnet_DCBM_est1 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_5_entries]*log(as.vector(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_5_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_5_entries])*log(as.vector(1-randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_5_entries]))), 3))
round(mean((randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_5_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_SP_5_entries])^2), 5)


paste("Negative log-likelihood of randnet_DCBM_est2 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_5_entries]*log(as.vector(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_5_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_5_entries])*log(as.vector(1-randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_5_entries]))), 3))

round(mean((randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_5_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_SSP_5_entries])^2), 5)


randnet_LR_num_clus = LRBIC(sparse_adj, 24)
randnet_LR_clus =  reg.SP(sparse_adj,K=randnet_LR_num_clus$DCSBM.K,lap=TRUE)
randnet_LR_clus2 =  reg.SSP(sparse_adj,K=randnet_LR_num_clus$DCSBM.K)

randnet_LR_DCBM_est1 =  DCSBM.estimate(sparse_adj, randnet_LR_clus$cluster)
randnet_LR_DCBM_est2 = DCSBM.estimate(sparse_adj, randnet_LR_clus2$cluster)

randnet_LR_DCBM_est1$Phat  = Adjust_within_community_estimates(randnet_LR_DCBM_est1$Phat, randnet_LR_clus$cluster)
randnet_LR_DCBM_est2$Phat  = Adjust_within_community_estimates(randnet_LR_DCBM_est2$Phat, randnet_LR_clus2$cluster)

range(randnet_LR_DCBM_est1$Phat)
range(randnet_LR_DCBM_est2$Phat)

valid_randnet_SP_16_entries = intersect(which(randnet_LR_DCBM_est1$Phat[upper.tri(randnet_LR_DCBM_est1$Phat)] <=1), which(randnet_LR_DCBM_est1$Phat[upper.tri(randnet_LR_DCBM_est1$Phat)] >=0))
valid_randnet_SSP_16_entries = intersect(which(randnet_LR_DCBM_est2$Phat[upper.tri(randnet_LR_DCBM_est2$Phat)] <=1), which(randnet_LR_DCBM_est2$Phat[upper.tri(randnet_LR_DCBM_est2$Phat)] >=0))


paste("Negative log-likelihood of randnet_LR_DCBM_est1 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_16_entries]*log(as.vector(randnet_LR_DCBM_est1$Phat[upper.tri(randnet_LR_DCBM_est1$Phat)][valid_randnet_SP_16_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_16_entries])*log(as.vector(1-randnet_LR_DCBM_est1$Phat[upper.tri(randnet_LR_DCBM_est1$Phat)][valid_randnet_SP_16_entries]))), 3))
round(mean((randnet_LR_DCBM_est1$Phat[upper.tri(randnet_LR_DCBM_est1$Phat)][valid_randnet_SP_16_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_SP_16_entries])^2), 5)

paste("Negative log-likelihood of randnet_LR_DCBM_est2 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_16_entries]*log(as.vector(randnet_LR_DCBM_est2$Phat[upper.tri(randnet_LR_DCBM_est2$Phat)][valid_randnet_SSP_16_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_16_entries])*log(as.vector(1-randnet_LR_DCBM_est2$Phat[upper.tri(randnet_LR_DCBM_est2$Phat)][valid_randnet_SSP_16_entries]))), 3))
round(mean((randnet_LR_DCBM_est2$Phat[upper.tri(randnet_LR_DCBM_est2$Phat)][valid_randnet_SSP_16_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_SSP_16_entries])^2), 5)


zh2 = spec_clust(adjacency_matrix, K=randnet_LR_num_clus$DCSBM.K)
randnet_LR_DCBM_est_using_nett_labels = DCSBM.estimate(sparse_adj, zh2)
randnet_LR_DCBM_est_using_nett_labels$Phat  = Adjust_within_community_estimates(randnet_LR_DCBM_est_using_nett_labels$Phat, zh2)





range(randnet_LR_DCBM_est_using_nett_labels$Phat)


valid_nett_16_entries = intersect(which(randnet_LR_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_LR_DCBM_est_using_nett_labels$Phat)] <=1), which(randnet_LR_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_LR_DCBM_est_using_nett_labels$Phat)] >=0))
paste("Negative log-likelihood of nett_16 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_16_entries]*log(as.vector(randnet_LR_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_LR_DCBM_est_using_nett_labels$Phat)][valid_nett_16_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_16_entries])*log(as.vector(1-randnet_LR_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_LR_DCBM_est_using_nett_labels$Phat)][valid_nett_16_entries]))), 3))
round(mean((randnet_LR_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_LR_DCBM_est_using_nett_labels$Phat)][valid_nett_16_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_nett_16_entries])^2), 5)



#randnet_LR_DCBM_est_using_nett_labels = DCSBM.estimate(sparse_adj, zh2)

randnet_true_comms =  DCSBM.estimate(sparse_adj, node_groups)
randnet_true_comms$Phat  = Adjust_within_community_estimates(randnet_true_comms$Phat, node_groups)

valid_randnet_true_comms = intersect(which(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)] <=1), which(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)] >=0))
paste("Negative log-likelihood of DCBM with true community estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)]))), 3))

round(mean((randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)][valid_randnet_true_comms]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_true_comms])^2), 5)


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




#check_PABM = PABM_estimate(adjacency_matrix, node_groups)
#check_PABM2 = PABM_estimate(adjacency_matrix, zh2)
##PABM Clustering using SSC
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
check_PABM = check_PABM = PABM_estimate(adjacency_matrix, SSCclust$cluster)
range(check_PABM[upper.tri(check_PABM)])


paste("Negative log-likelihood of PABM estimate using SSC communities is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(check_PABM[upper.tri(check_PABM)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-check_PABM[upper.tri(check_PABM)]))), 3))

valid_PABM_entries = intersect(which(check_PABM[upper.tri(check_PABM)] <=1), which(check_PABM[upper.tri(check_PABM)] >=0))



paste("Negative log-likelihood of PABM estimate using SSC communities is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries]*log(as.vector(check_PABM[upper.tri(check_PABM)][valid_PABM_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries])*log(as.vector(1-check_PABM[upper.tri(check_PABM)][valid_PABM_entries]))), 3))

paste("Negative log-likelihood of our estimate using true communities using only those valid values of PABM is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries]*log(as.vector(reconstructedEstimates[[1]][upper.tri(reconstructedEstimates[[1]])][valid_PABM_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries])*log(as.vector(1-reconstructedEstimates[[1]][upper.tri(reconstructedEstimates[[1]])][valid_PABM_entries]))), 3))
round(mean((check_PABM[upper.tri(check_PABM)][valid_PABM_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_PABM_entries])^2), 5)

paste("Negative log-likelihood of no_error_probs_matrix using true communities using only those valid values of PABM is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries]*log(as.vector(no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_PABM_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries])*log(as.vector(1-no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_PABM_entries]))), 3))


check_PABM_truncated = check_PABM
check_PABM_truncated[which(check_PABM_truncated>=1)] = .999
check_PABM_truncated[which(check_PABM_truncated<=0)] = .001
length(intersect(which(check_PABM_truncated[upper.tri(check_PABM_truncated)] <=1), which(check_PABM_truncated[upper.tri(check_PABM_truncated)] >=0)))
paste("Negative log-likelihood of truncated PABM estimate using SSC is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(check_PABM_truncated[upper.tri(check_PABM_truncated)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-check_PABM_truncated[upper.tri(check_PABM_truncated)]))), 3))
round(mean((check_PABM_truncated[upper.tri(check_PABM_truncated)]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2), 5)



subnetwork_degree_equality_checker = function(A, B, g)
{
  K <- length(unique(g))
  unequal_degrees=0
  #B <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      N.i <- which(g == i)
      N.j <- which(g == j)
      A_prime = sum(A[N.i, N.j])
      B_prime  <- sum(B[N.i, N.j])
      if(abs(A_prime-B_prime)>.1){unequal_degrees=unequal_degrees+1}
    }
  }
  
  return(unequal_degrees)
}



subnetwork_degree_equality_checker(adjacency_matrix, check_PABM, node_groups)



randnetDCBMforplot  = randnet_DCBM_est_using_nett_labels$Phat
randnetDCBMforplot[which(randnetDCBMforplot>1)] =NA 

PABMforplot  = check_PABM
PABMforplot[which(PABMforplot>1)] =NA 


image(randnetDCBMforplot, col = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(20))
image(PABMforplot, col = colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))(20))




#table 1
#P* is no_error_probs_matrix

length(intersect(which(no_error_probs_matrix[upper.tri(no_error_probs_matrix)]>=0), which(no_error_probs_matrix[upper.tri(no_error_probs_matrix)]<=1)))
round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(no_error_probs_matrix[upper.tri(no_error_probs_matrix)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-no_error_probs_matrix[upper.tri(no_error_probs_matrix)]))), 3)

#P is probs_matrix
length(intersect(which(probs_matrix[upper.tri(probs_matrix)]>=0), which(probs_matrix[upper.tri(probs_matrix)]<=1)))
round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(probs_matrix[upper.tri(probs_matrix)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-probs_matrix[upper.tri(probs_matrix)]))), 3)
round(mean((probs_matrix[upper.tri(probs_matrix)]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2), 5)


#P-tilde is reconstructedEstimates[[1]]

length(intersect(which(reconstructedEstimates[[1]][upper.tri(reconstructedEstimates[[1]])]>=0), which(reconstructedEstimates[[1]][upper.tri(reconstructedEstimates[[1]])]<=1)))
round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(reconstructedEstimates[[1]][upper.tri(reconstructedEstimates[[1]])]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-reconstructedEstimates[[1]][upper.tri(reconstructedEstimates[[1]])]))), 3)
round(mean((reconstructedEstimates[[1]][upper.tri(reconstructedEstimates[[1]])]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2), 5)




#DCBM true communities is below
paste("Negative log-likelihood of DCBM with true community estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)]))), 3))
round(mean((randnet_true_comms$Phat[upper.tri(randnet_true_comms$Phat)][valid_randnet_true_comms]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_true_comms])^2), 5)

#PABM SSC communities is below
length(valid_PABM_entries)
paste("Negative log-likelihood of PABM estimate using SSC communities is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries]*log(as.vector(check_PABM[upper.tri(check_PABM)][valid_PABM_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries])*log(as.vector(1-check_PABM[upper.tri(check_PABM)][valid_PABM_entries]))), 3))

paste("Negative log-likelihood of our estimate using true communities using only those valid values of PABM is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries]*log(as.vector(reconstructedEstimates[[1]][upper.tri(reconstructedEstimates[[1]])][valid_PABM_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries])*log(as.vector(1-reconstructedEstimates[[1]][upper.tri(reconstructedEstimates[[1]])][valid_PABM_entries]))), 3))
round(mean((check_PABM[upper.tri(check_PABM)][valid_PABM_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_PABM_entries])^2), 5)


#P* at valid values of PABM
paste("Negative log-likelihood of no_error_probs_matrix using true communities using only those valid values of PABM is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries]*log(as.vector(no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_PABM_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_PABM_entries])*log(as.vector(1-no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_PABM_entries]))), 3))

#PABM true communities truncated is below
paste("Negative log-likelihood of truncated PABM estimate using true communities is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)]*log(as.vector(check_PABM_truncated[upper.tri(check_PABM_truncated)]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)])*log(as.vector(1-check_PABM_truncated[upper.tri(check_PABM_truncated)]))), 3))
round(mean((check_PABM_truncated[upper.tri(check_PABM_truncated)]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2), 5)

#DCBM with spectral clustering with regularization(nett)
length(valid_nett_5_entries)
paste("Negative log-likelihood of DCBM estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_5_entries]*log(as.vector(randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_5_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_5_entries])*log(as.vector(1-randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_5_entries]))), 3))
round(mean((randnet_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_DCBM_est_using_nett_labels$Phat)][valid_nett_5_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_nett_5_entries])^2), 5)



#DCBM with regularized spectral clustering(randnet) 

length(valid_randnet_SP_5_entries)


paste("Negative log-likelihood of randnet_DCBM_est1 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_5_entries]*log(as.vector(randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_5_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_5_entries])*log(as.vector(1-randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_5_entries]))), 3))
round(mean((randnet_DCBM_est1$Phat[upper.tri(randnet_DCBM_est1$Phat)][valid_randnet_SP_5_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_SP_5_entries])^2), 5)

#DCBM with regularized spherical spectral clustering(randnet) 

length(valid_randnet_SSP_5_entries)
paste("Negative log-likelihood of randnet_DCBM_est2 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_5_entries]*log(as.vector(randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_5_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_5_entries])*log(as.vector(1-randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_5_entries]))), 3))

round(mean((randnet_DCBM_est2$Phat[upper.tri(randnet_DCBM_est2$Phat)][valid_randnet_SSP_5_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_SSP_5_entries])^2), 5)


#DCBM with spectral clustering with regularization(nett)
length(valid_nett_16_entries)
paste("Negative log-likelihood of nett_16 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_16_entries]*log(as.vector(randnet_LR_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_LR_DCBM_est_using_nett_labels$Phat)][valid_nett_16_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_nett_16_entries])*log(as.vector(1-randnet_LR_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_LR_DCBM_est_using_nett_labels$Phat)][valid_nett_16_entries]))), 3))
round(mean((randnet_LR_DCBM_est_using_nett_labels$Phat[upper.tri(randnet_LR_DCBM_est_using_nett_labels$Phat)][valid_nett_16_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_nett_16_entries])^2), 5)


#DCBM with regularized spectral clustering(randnet)
length(valid_randnet_SP_16_entries)
paste("Negative log-likelihood of randnet_LR_DCBM_est1 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_16_entries]*log(as.vector(randnet_LR_DCBM_est1$Phat[upper.tri(randnet_LR_DCBM_est1$Phat)][valid_randnet_SP_16_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SP_16_entries])*log(as.vector(1-randnet_LR_DCBM_est1$Phat[upper.tri(randnet_LR_DCBM_est1$Phat)][valid_randnet_SP_16_entries]))), 3))
round(mean((randnet_LR_DCBM_est1$Phat[upper.tri(randnet_LR_DCBM_est1$Phat)][valid_randnet_SP_16_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_SP_16_entries])^2), 5)



#DCBM with regularized spherical spectral clustering(randnet)
length(valid_randnet_SSP_16_entries)
paste("Negative log-likelihood of randnet_LR_DCBM_est2 estimate is",  round(-sum(adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_16_entries]*log(as.vector(randnet_LR_DCBM_est2$Phat[upper.tri(randnet_LR_DCBM_est2$Phat)][valid_randnet_SSP_16_entries]))) -sum((1-adjacency_matrix[upper.tri(adjacency_matrix)][valid_randnet_SSP_16_entries])*log(as.vector(1-randnet_LR_DCBM_est2$Phat[upper.tri(randnet_LR_DCBM_est2$Phat)][valid_randnet_SSP_16_entries]))), 3))
round(mean((randnet_LR_DCBM_est2$Phat[upper.tri(randnet_LR_DCBM_est2$Phat)][valid_randnet_SSP_16_entries]- no_error_probs_matrix[upper.tri(no_error_probs_matrix)][valid_randnet_SSP_16_entries])^2), 5)




###timing

sum(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))/60
max(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))/60

































##################################
##OTHER ANALYSIS - comparing community detection methods
##################################




UFullLikelihood =  -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-reconstructedEstimates[[1]][lower.tri(reconstructedEstimates[[1]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])


BlindMSEs = rep(0, 10)
BlindLikelihoods = rep(0, 10)

for(i in 1:10){
  BlindMSEs[i] = mean((BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])] - no_error_probs_matrix[lower.tri(no_error_probs_matrix)])^2)
  BlindLikelihoods[i] = -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])
}



reconstructedEstimatesU = reconstructedEstimates
UBlindreconstructedEstimates = BlindreconstructedEstimates
UBlindMSEs = BlindMSEs
UBlindLikelihoods = BlindLikelihoods 

Utotal_runtime = total_runtime
UFullyObservedEstTime = sum(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))
UFullyObservedRealisticTime = max(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))

U31Runtimes = unlist(lapply(1:31, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
U31RealisticTime = unlist(lapply(1:31, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))








load("SimFitV.RData")

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





BlindMSEs = rep(0, 10)
BlindLikelihoods = rep(0, 10)

for(i in 1:10){
  BlindMSEs[i] = mean((BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])] - no_error_probs_matrix[lower.tri(no_error_probs_matrix)])^2)
  BlindLikelihoods[i] = -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])
}



reconstructedEstimatesV = reconstructedEstimates
VBlindreconstructedEstimates = BlindreconstructedEstimates
VBlindMSEs = BlindMSEs
VBlindLikelihoods = BlindLikelihoods 

Vtotal_runtime = total_runtime
if(VRepeats[1]==1){
  VFullyObservedEstTime = UFullyObservedEstTime
  VFullyObservedRealisticTime = UFullyObservedRealisticTime
}else{
VFullyObservedEstTime =  sum(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))
VFullyObservedRealisticTime = max(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))
}


V31Runtimes = unlist(lapply(1:31, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
V31RealisticTime = unlist(lapply(1:31, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
for(k in 1:31){
  if(VRepeats[k] == 1){
    V31Runtimes[k] = U31Runtimes[k]
    V31RealisticTime[k] = U31RealisticTime[k]
  }
}






load("SimFitEig.RData")

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





BlindMSEs = rep(0, 10)
BlindLikelihoods = rep(0, 10)

for(i in 1:10){
  BlindMSEs[i] = mean((BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])] - no_error_probs_matrix[lower.tri(no_error_probs_matrix)])^2)
  BlindLikelihoods[i] = -sum((as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])* log(as.vector(BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))]) -sum(((1-as.vector(adjacency_matrix[lower.tri(adjacency_matrix)]))*log(as.vector(1-BlindreconstructedEstimates[[i]][lower.tri(BlindreconstructedEstimates[[i]])])))[!(is.na(as.vector(adjacency_matrix[lower.tri(adjacency_matrix)])))])
}



reconstructedEstimatesEig = reconstructedEstimates
EigBlindreconstructedEstimates = BlindreconstructedEstimates
EigBlindMSEs = BlindMSEs
EigBlindLikelihoods = BlindLikelihoods 


Eigtotal_runtime = total_runtime

if(EigRepeats[1]==1){
  EigFullyObservedEstTime = UFullyObservedEstTime
  EigFullyObservedRealisticTime = UFullyObservedRealisticTime
}else{if(EigRepeats[1]==2){
  EigFullyObservedEstTime = VFullyObservedEstTime
  EigFullyObservedRealisticTime = VFullyObservedRealisticTime
}else{
  EigFullyObservedEstTime =  sum(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))
  EigFullyObservedRealisticTime = max(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))
}
}


Eig31Runtimes = unlist(lapply(1:31, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
Eig31RealisticTime = unlist(lapply(1:31, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
for(k in 1:31){
  if(EigRepeats[k] == 1){
    Eig31Runtimes[k] = U31Runtimes[k]
    Eig31RealisticTime[k] = U31RealisticTime[k]
  }
  if(EigRepeats[k] == 2){
    Eig31Runtimes[k] = V31Runtimes[k]
    Eig31RealisticTime[k] = V31RealisticTime[k]
  }
}



################################
####Plotting community detection results in different ways 
################################


plot(jitter(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsEig[[x]])})), amount = .2),
     sqrt(unlist(lapply(2:31, function(x){
       mean((reconstructedEstimatesEig[[x]][upper.tri(no_error_probs_matrix)] 
             -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)}
     ))),
     col='black', pch=16, ylim= c(0.03, .25), xlim = c(0, 8), , xlab = "Number of Estimated Communities", ylab = "RMSE")

points(jitter(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsU[[x]])})), amount = .2),
     sqrt(unlist(lapply(2:31, function(x){
       mean((reconstructedEstimatesU[[x]][upper.tri(no_error_probs_matrix)] 
               -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)}
       ))),
     col='red', pch=25)


points(jitter(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsV[[x]])})), amount = .2),
       sqrt(unlist(lapply(2:31, function(x){
         mean((reconstructedEstimatesV[[x]][upper.tri(no_error_probs_matrix)] 
               -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)}
       ))),
       col='green', pch=4)

legend("bottomleft",  title="Method",
       c("Eigenvector",  "U", "V"), col= c("black","red", "green"), cex=0.8, pch = c(16, 25, 4))



#plot(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsEig[[x]])})), unlist(lapply(2:31, function(x){mean((reconstructedEstimatesEig[[x]][upper.tri(no_error_probs_matrix)] -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)})), pch=16)
#points(1:30, unlist(lapply(2:31, function(x){max(all_estimated_node_groupsU[[x]])})), col='red', pch=25)
#points(1:30, unlist(lapply(2:31, function(x){max(all_estimated_node_groupsV[[x]])})), col='green', pch=4)

##########Adjusted rand indices

library(mclust)
EigARI = unlist(lapply(2:31, function(x){adjustedRandIndex(all_estimated_node_groupsEig[[x]], node_groups)}))
UARI = unlist(lapply(2:31, function(x){adjustedRandIndex(all_estimated_node_groupsU[[x]], node_groups)}))
VARI = unlist(lapply(2:31, function(x){adjustedRandIndex(all_estimated_node_groupsV[[x]], node_groups)}))


plot(EigARI,  sqrt(unlist(lapply(2:31, function(x){
  mean((reconstructedEstimatesEig[[x]][upper.tri(no_error_probs_matrix)] 
        -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)}
  ))),
pch = as.character(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsEig[[x]])}))), ylim=c(.04, .22), xlab = "Adjusted Rand Index with True Communities", ylab= "RMSE", cex = .6 )


points(UARI,  sqrt(unlist(lapply(2:31, function(x){
  mean((reconstructedEstimatesU[[x]][upper.tri(no_error_probs_matrix)] 
        -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)}
))),
pch = as.character(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsU[[x]])}))), col='red', cex= .6)

points(VARI,  sqrt(unlist(lapply(2:31, function(x){
  mean((reconstructedEstimatesV[[x]][upper.tri(no_error_probs_matrix)] 
        -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)}
))),
pch = as.character(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsV[[x]])}))), col='green', cex=.6)

legend("bottomleft",  title="Method",
       c("Eigenvector",  "U", "V"), fill= c("black","red", "green"), cex=0.8)





plot(     sqrt(unlist(lapply(2:31, function(x){
       mean((reconstructedEstimatesEig[[x]][upper.tri(no_error_probs_matrix)] 
             -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)}
     ))),
     col='black', cex= .7, pch=as.character(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsEig[[x]])}))), ylim= c(0, .3),   xlab = "Replication Number", ylab = "RMSE")

points(1:30,
       sqrt(unlist(lapply(2:31, function(x){
         mean((reconstructedEstimatesU[[x]][upper.tri(no_error_probs_matrix)] 
               -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)}
       ))),
       col='red',cex= .7, pch=as.character(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsU[[x]])}))))


points(1:30,
       sqrt(unlist(lapply(2:31, function(x){
         mean((reconstructedEstimatesV[[x]][upper.tri(no_error_probs_matrix)] 
               -no_error_probs_matrix[upper.tri(no_error_probs_matrix)])^2)}
       ))),
       col='green',cex= .7, pch=as.character(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsV[[x]])}))))

legend("topright",
       c("Eigenvector",  "U", "V"), col= c("black","red", "green"), cex=0.6, pch = c(16, 25, 4))




################################
####Table for community detection results in different ways 
################################


UFullLikelihood
VFullLikelihood
EigFullLikelihood
mean(UBlindLikelihoods)
mean(VBlindLikelihoods)
mean(EigBlindLikelihoods)
mean(UBlindMSEs)
mean(VBlindMSEs)
mean(EigBlindMSEs)

mean(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsU[[x]])})))
mean(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsV[[x]])})))
mean(unlist(lapply(2:31, function(x){max(all_estimated_node_groupsEig[[x]])})))

mean(UARI)
mean(VARI)
mean(EigARI)
sum(UARI ==1)
sum(VARI ==1)
sum(EigARI ==1)



################################
####Making Ajusted Rand index table
################################




tocomp = cbind(node_groups, all_estimated_node_groupsV[[1]], SSCclust$cluster, 
               zh, randnet_clus$cluster, randnet_clus2$cluster,  zh2, 
               randnet_LR_clus$cluster, randnet_LR_clus2$cluster
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



################################
####More runtimes from simulations
################################

Utotal_runtime*60
sum(U31Runtimes)
mean(U31RealisticTime)


Vtotal_runtime*60
sum(V31Runtimes)
mean(V31RealisticTime)

Eigtotal_runtime*60
sum(Eig31Runtimes)
mean(Eig31RealisticTime)
