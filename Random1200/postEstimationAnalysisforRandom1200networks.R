load("Redo_DOROTHY_RandomAcronym30NestedFULL2.RData")
library(mclust)


all_estimated_node_groupsEig = all_estimated_node_groups
EigMSEs = all_matrix_MSEs
EigARI = unlist(lapply(1:30, function(x){adjustedRandIndex(all_estimated_node_groups[[x]], all_node_groups[[x]])}))
Eigtotal_runtime = total_runtime #14.27036*3600
Eig31Runtimes = unlist(lapply(1:30, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))
Eig31RealisticTime = unlist(lapply(1:30, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput2[[y]]), function(x){parallelOutput2[[y]][[x]][[12]]}))))}))


load("Redo_SingularValueTest.RData")
all_estimated_node_groupsU = all_estimated_node_groups_SVD
UMSEs = all_matrix_MSEs_SVD
UARI = unlist(lapply(1:30, function(x){adjustedRandIndex(all_estimated_node_groups_SVD[[x]], all_node_groups[[x]])}))
Utotal_runtime = total_runtime
U31Runtimes = unlist(lapply(1:30, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput_SVD2[[y]]), function(x){parallelOutput_SVD2[[y]][[x]][[12]]}))))}))
U31RealisticTime = unlist(lapply(1:30, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput_SVD2[[y]]), function(x){parallelOutput_SVD2[[y]][[x]][[12]]}))))}))


for(k in 1:30){
  if(all.equal(all_estimated_node_groups_SVD[[k]], all_estimated_node_groups[[k]])==1){
    U31Runtimes[k]= Eig31Runtimes[k]
    U31RealisticTime[k]= Eig31RealisticTime[k]
    
  
}
}


load("Redo_SingularValueTestV.RData")
all_estimated_node_groupsV = all_estimated_node_groups_SVD
VMSEs = all_matrix_MSEs_SVD
VARI = unlist(lapply(1:30, function(x){adjustedRandIndex(all_estimated_node_groups_SVD[[x]], all_node_groups[[x]])}))
Vtotal_runtime = total_runtime
V31Runtimes = unlist(lapply(1:30, function(y){sum(na.omit(unlist(lapply(1:length(parallelOutput_SVD2[[y]]), function(x){parallelOutput_SVD2[[y]][[x]][[12]]}))))}))
V31RealisticTime = unlist(lapply(1:30, function(y){max(na.omit(unlist(lapply(1:length(parallelOutput_SVD2[[y]]), function(x){parallelOutput_SVD2[[y]][[x]][[12]]}))))}))


for(k in 1:30){
  if(all.equal(all_estimated_node_groups_SVD[[k]], all_estimated_node_groups[[k]])==1){
    V31Runtimes[k]= Eig31Runtimes[k]
    V31RealisticTime[k]= Eig31RealisticTime[k]
    
    
  }
}


for(k in 1:30){
  if(all.equal(all_estimated_node_groupsV[[k]], all_estimated_node_groupsEig[[k]])==1){
    VMSEs[k]= EigMSEs[k]
  }
  if(all.equal(all_estimated_node_groupsU[[k]], all_estimated_node_groupsEig[[k]])==1){
    UMSEs[k]= EigMSEs[k]
  }
}



mean(EigMSEs)
mean(VMSEs)
mean(UMSEs)

mean(unlist(lapply(1:30, function(x){max(all_estimated_node_groupsU[[x]])})))
mean(unlist(lapply(1:30, function(x){max(all_estimated_node_groupsV[[x]])})))
mean(unlist(lapply(1:30, function(x){max(all_estimated_node_groupsEig[[x]])})))

mean(UARI)
mean(VARI)
mean(EigARI)
sum(UARI ==1)
sum(VARI ==1)
sum(EigARI ==1)


Utotal_runtime*60
sum(U31Runtimes)/60
mean(U31RealisticTime)/60

Vtotal_runtime*60
sum(V31Runtimes)/60
mean(V31RealisticTime)/60

Eigtotal_runtime*60
sum(Eig31Runtimes)/60
mean(Eig31RealisticTime)/60

EigFullyObservedEstTime = sum(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))
EigFullyObservedRealisticTime = max(unlist(lapply(1:length(parallelOutput2[[1]]), function(x){parallelOutput2[[1]][[x]][[12]]})))













plot( sqrt(EigMSEs),
col='black',  pch=as.character(unlist(lapply(1:30, function(x){max(all_estimated_node_groupsEig[[x]])}))), ylim= c(0, .25),   xlab = "Replication Number", ylab = "RMSE")

points(1:30,
       sqrt(UMSEs),
       col='red', pch=as.character(unlist(lapply(1:30, function(x){max(all_estimated_node_groupsU[[x]])}))
       ))


points(1:30,
       sqrt(VMSEs),
       col='green', pch=as.character(unlist(lapply(1:30, function(x){max(all_estimated_node_groupsV[[x]])}))))

legend("topright",
       c("Eigenvector",  "U", "V"), col= c("black","red", "green"), fill=c("black","red", "green"), cex=0.7)





plot(unlist(lapply(1:30, function(x){max(all_estimated_node_groups[[x]])})), lapply(1:30, function(x){adjustedRandIndex(all_estimated_node_groups[[x]], all_node_groups[[x]])}))


plot(unlist(lapply(1:30, function(x){max(all_estimated_node_groups[[x]])})), all_matrix_MSEs)

plot(unlist(lapply(1:30, function(x){adjustedRandIndex(all_estimated_node_groups[[x]], all_node_groups[[x]])})), all_matrix_MSEs)

