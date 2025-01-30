load("Redo_DOROTHY_Acronym200.RData")

check_result = function(x){
  if(mean(all_adjacency_matrix[[x]])<.5){sqrt(mean((parallelOutput2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
  else{sqrt(mean((1-parallelOutput2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
}


check_result_psi = function(x){
  if(mean(all_adjacency_matrix[[x]])<.5){sqrt(mean((c(parallelOutput2[[x]][[5]], parallelOutput2[[x]][[6]])- node_sociabilities)^2))}
  else{sqrt(mean((1-c(parallelOutput2[[x]][[5]], parallelOutput2[[x]][[6]])- node_sociabilities)^2))}
}


check_resultNORMAL = function(x){
  if(mean(all_adjacency_matrix[[x]])<.5){sqrt(mean((parallelOutputNORMAL2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
  else{sqrt(mean((1-parallelOutputNORMAL2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
}


check_resultNORMAL_psi = function(x){
  if(mean(all_adjacency_matrix[[x]])<.5){sqrt(mean((c(parallelOutputNORMAL2[[x]][[5]], parallelOutputNORMAL2[[x]][[6]])- node_sociabilities)^2))}
  else{sqrt(mean((1-c(parallelOutputNORMAL2[[x]][[5]], parallelOutputNORMAL2[[x]][[6]])- node_sociabilities)^2))}

}


mean(unlist(lapply(1:240, function(x){parallelOutput2[[x]][[7]]==HForms[[x]]})))
correct_ests = which(unlist(lapply(1:240, function(x){parallelOutput2[[x]][[7]] == HForms[x] })) )



HAccuracy = rep(0, 12)
MSEs = rep(0, 12)
for(setting in 1:12){
  relevant_runs = ((20*(setting-1))+1):(20*setting)
  HAccuracy[setting] = mean(unlist(lapply(relevant_runs, function(x){parallelOutput2[[x]][[7]]==HForms[[x]]})))
  MSEs[setting] = mean(unlist(lapply(relevant_runs , check_result)))
  
}


firstcolumn = c(mean(HAccuracy),
                mean(HAccuracy[c(1, 4, 7, 10)]),
                mean(HAccuracy[c(1, 4, 7, 10)+1]),
                mean(HAccuracy[c(1, 4, 7, 10)+2]),
                mean(HAccuracy[1:3]),
                mean(HAccuracy[4:6]),
                mean(HAccuracy[7:9]),
                mean(HAccuracy[10:12])
)


second_column = c(mean(MSEs),
                  mean(MSEs[c(1, 4, 7, 10)]),
                  mean(MSEs[c(1, 4, 7, 10)+1]),
                  mean(MSEs[c(1, 4, 7, 10)+2]),
                  mean(MSEs[1:3]),
                  mean(MSEs[4:6]),
                  mean(MSEs[7:9]),
                  mean(MSEs[10:12])
) 

third_column = c(mean(unlist(lapply( correct_ests , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+20) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+40) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, 1:60) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, 61:120) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, 121:180) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, 181:240) , check_result)))
)


fourth_column = c( mean(unlist(lapply(setdiff(1:240, correct_ests) , check_result))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)) , check_result))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+20) , check_result))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+40) , check_result))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 1:60) , check_result))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 61:120) , check_result))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 121:180) , check_result))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 181:240) , check_result)))
)


fifth_column = c( mean(unlist(lapply( 1:240 , check_result_psi))),
                  mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200) , check_result_psi))),
                  mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200)+20 , check_result_psi))),
                  mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200)+40 , check_result_psi))),
                  mean(unlist(lapply( 1:60 , check_result_psi))),
                  mean(unlist(lapply( 61:120 , check_result_psi))),
                  mean(unlist(lapply( 121:180 , check_result_psi))),
                  mean(unlist(lapply( 181:240 , check_result_psi)))
)


sixth_column = c( mean(unlist(lapply( correct_ests , check_result_psi))), 
                  mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)) , check_result_psi))),
                  mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+20) , check_result_psi))),
                  mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+40) , check_result_psi))),
                  mean(unlist(lapply( intersect(correct_ests, 1:60) , check_result_psi))),
                  mean(unlist(lapply( intersect(correct_ests, 61:120) , check_result_psi))),
                  mean(unlist(lapply( intersect(correct_ests, 121:180) , check_result_psi))),
                  mean(unlist(lapply( intersect(correct_ests, 181:240) , check_result_psi)))
)


seventh_column = c(mean(unlist(lapply(setdiff(1:240, correct_ests) , check_result_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)) , check_result_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+20) , check_result_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+40) , check_result_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 1:60) , check_result_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 61:120) , check_result_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 121:180) , check_result_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 181:240) , check_result_psi)))
)



#comparison to DCBM/PABM
DCBM_symmetric_subnetwork_estimate = function(A)
{
  #B <- matrix(0, K, K)
  P.hat <- 0*A
  i=1 
  j=2
  N.i <- nrow(A)
  N.j <- ncol(A)
  A_prime = A
  B  <- sum(A_prime, na.rm = TRUE) + 
    0.001
  Psi <- rowSums(A_prime, na.rm = TRUE)
  Psi <- Psi/sum(Psi)
  Psi2 <- colSums(A_prime, na.rm = TRUE)
  Psi2 <- Psi2/sum(Psi2)
  P.hat = B*Psi%*%t(Psi2)
}



MSE_DCBM = rep(0, length(all_adjacency_matrix))
for(netnum in 1:length(all_adjacency_matrix)){
  
  MSE_DCBM[[netnum]] = sqrt(mean((DCBM_symmetric_subnetwork_estimate(all_adjacency_matrix[[netnum]]) - all_no_error_probs_matrix[[netnum]])^2))
  
}


mean(MSE_DCBM)




DCBM_column = c(mean(MSE_DCBM),
                mean(MSE_DCBM[c(1:20, 61:80, 121:140, 181:200)]),
                mean(MSE_DCBM[c(1:20, 61:80, 121:140, 181:200)+20]),
                mean(MSE_DCBM[c(1:20, 61:80, 121:140, 181:200)+40]),
                mean(MSE_DCBM[1:60]),
                mean(MSE_DCBM[61:120]),
                mean(MSE_DCBM[121:180]),
                mean(MSE_DCBM[181:240])
) 
















correct_ests = which(unlist(lapply(1:240, function(x){parallelOutputNORMAL2[[x]][[7]] == HForms[x] })) )



HAccuracy = rep(0, 12)
MSEs = rep(0, 12)
for(setting in 1:12){
  relevant_runs = ((20*(setting-1))+1):(20*setting)
  HAccuracy[setting] = mean(unlist(lapply(relevant_runs, function(x){parallelOutputNORMAL2[[x]][[7]]==HForms[[x]]})))
  MSEs[setting] = mean(unlist(lapply(relevant_runs , check_resultNORMAL)))
  
}


firstcolumn = c(mean(HAccuracy),
                mean(HAccuracy[c(1, 4, 7, 10)]),
                mean(HAccuracy[c(1, 4, 7, 10)+1]),
                mean(HAccuracy[c(1, 4, 7, 10)+2]),
                mean(HAccuracy[1:3]),
                mean(HAccuracy[4:6]),
                mean(HAccuracy[7:9]),
                mean(HAccuracy[10:12])
)


second_column = c(mean(MSEs),
                  mean(MSEs[c(1, 4, 7, 10)]),
                  mean(MSEs[c(1, 4, 7, 10)+1]),
                  mean(MSEs[c(1, 4, 7, 10)+2]),
                  mean(MSEs[1:3]),
                  mean(MSEs[4:6]),
                  mean(MSEs[7:9]),
                  mean(MSEs[10:12])
) 

third_column = c(mean(unlist(lapply( correct_ests , check_resultNORMAL))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)) , check_resultNORMAL))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+20) , check_resultNORMAL))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+40) , check_resultNORMAL))),
                 mean(unlist(lapply( intersect(correct_ests, 1:60) , check_resultNORMAL))),
                 mean(unlist(lapply( intersect(correct_ests, 61:120) , check_resultNORMAL))),
                 mean(unlist(lapply( intersect(correct_ests, 121:180) , check_resultNORMAL))),
                 mean(unlist(lapply( intersect(correct_ests, 181:240) , check_resultNORMAL)))
)


fourth_column = c( mean(unlist(lapply(setdiff(1:240, correct_ests) , check_resultNORMAL))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)) , check_resultNORMAL))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+20) , check_resultNORMAL))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+40) , check_resultNORMAL))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 1:60) , check_resultNORMAL))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 61:120) , check_resultNORMAL))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 121:180) , check_resultNORMAL))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 181:240) , check_resultNORMAL)))
)


fifth_column = c( mean(unlist(lapply( 1:240 , check_resultNORMAL_psi))),
                  mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200) , check_resultNORMAL_psi))),
                  mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200)+20 , check_resultNORMAL_psi))),
                  mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200)+40 , check_resultNORMAL_psi))),
                  mean(unlist(lapply( 1:60 , check_resultNORMAL_psi))),
                  mean(unlist(lapply( 61:120 , check_resultNORMAL_psi))),
                  mean(unlist(lapply( 121:180 , check_resultNORMAL_psi))),
                  mean(unlist(lapply( 181:240 , check_resultNORMAL_psi)))
)


sixth_column = c( mean(unlist(lapply( correct_ests , check_resultNORMAL_psi))), 
                  mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)) , check_resultNORMAL_psi))),
                  mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+20) , check_resultNORMAL_psi))),
                  mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+40) , check_resultNORMAL_psi))),
                  mean(unlist(lapply( intersect(correct_ests, 1:60) , check_resultNORMAL_psi))),
                  mean(unlist(lapply( intersect(correct_ests, 61:120) , check_resultNORMAL_psi))),
                  mean(unlist(lapply( intersect(correct_ests, 121:180) , check_resultNORMAL_psi))),
                  mean(unlist(lapply( intersect(correct_ests, 181:240) , check_resultNORMAL_psi)))
)


seventh_column = c(mean(unlist(lapply(setdiff(1:240, correct_ests) , check_resultNORMAL_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)) , check_resultNORMAL_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+20) , check_resultNORMAL_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+40) , check_resultNORMAL_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 1:60) , check_resultNORMAL_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 61:120) , check_resultNORMAL_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 121:180) , check_resultNORMAL_psi))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 181:240) , check_resultNORMAL_psi)))
)









plot(unlist(lapply(1:240 , check_result)),  col = rep(1:4, each=60), pch = (4*unlist(lapply(1:240, function(x){parallelOutput2[[x]][[7]] == HForms[x] }))), cex= .7,  main =paste("Estimation performance for subnetworks"), ylab ="RMSE"  , xaxt = "n", xlab ="", xlim = c(5, 235))

axis(1, at = seq(30, 210, by = 60), labels = c("Normal", "Linear", "Concave", "Convex"), mgp=c(0,.2,0))
axis(1, at = seq(10, 230, by = 20), labels = c("\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3"), lwd.tick = 0, mgp =c(0, 2,0)  )
abline(v = c(seq(20, 220, by=20)), lty=2, lwd = .5)


legend("topright",  c("Wrong H",  "Correct H"), pch = c(0, 4), cex=.7)

plot(unlist(lapply(1:240 , check_result_psi)),  col = rep(1:4, each=60), pch = (4*unlist(lapply(1:240, function(x){parallelOutput2[[x]][[7]] == HForms[x] }))), cex= .7,  main =paste("Psi Estimation MSE for subnetworks"), ylab ="RMSE"  , xaxt = "n", xlab ="", xlim = c(5, 235))

axis(1, at = seq(30, 210, by = 60), labels = c("Normal", "Linear", "Concave", "Convex"), mgp=c(0,.2,0))
axis(1, at = seq(10, 230, by = 20), labels = c("\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3"), lwd.tick = 0, mgp =c(0, 2,0)  )
abline(v = c(seq(20, 220, by=20)), lty=2, lwd = .5)


legend("topleft",  c("Wrong H",  "Correct H"), pch = c(0, 4), cex=.5)





plot(unlist(lapply(1:240 , check_resultNORMAL)),  col = rep(1:4, each=60), pch = (4*unlist(lapply(1:240, function(x){parallelOutputNORMAL2[[x]][[7]] == HForms[x] }))), cex= .7,  main =paste("Estimation performance for subnetworks using only Normal H-functions"), ylab ="RMSE"  , xaxt = "n", xlab ="", xlim = c(5, 235))

axis(1, at = seq(30, 210, by = 60), labels = c("Normal", "Linear", "Concave", "Convex"), mgp=c(0,.2,0))
axis(1, at = seq(10, 230, by = 20), labels = c("\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3"), lwd.tick = 0, mgp =c(0, 2,0)  )
abline(v = c(seq(20, 220, by=20)), lty=2, lwd = .5)


legend("topleft",  c("Wrong H",  "Correct H"), pch = c(0, 4), cex=.7)

plot(unlist(lapply(1:240 , check_resultNORMAL_psi)),  col = rep(1:4, each=60), pch = (4*unlist(lapply(1:240, function(x){parallelOutputNORMAL2[[x]][[7]] == HForms[x] }))), cex= .7,  main =paste("Psi Estimation MSE for subnetworks  using only Normal H-functions"), ylab ="RMSE"  , xaxt = "n", xlab ="", xlim = c(5, 235))

axis(1, at = seq(30, 210, by = 60), labels = c("Normal", "Linear", "Concave", "Convex"), mgp=c(0,.2,0))
axis(1, at = seq(10, 230, by = 20), labels = c("\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3"), lwd.tick = 0, mgp =c(0, 2,0)  )
abline(v = c(seq(20, 220, by=20)), lty=2, lwd = .5)


legend("topleft",  c("Wrong H",  "Correct H"), pch = c(0, 4), cex=.5)















load("Redo_DOROTHY_Acronym200Sigma0.RData")

check_result_zero = function(x){
  if(mean(all_adjacency_matrix[[x]])<.5){sqrt(mean((parallelOutput_sigma_zero2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
  else{sqrt(mean((1-parallelOutput_sigma_zero2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
}


check_result_psi_zero = function(x){
  if(mean(all_adjacency_matrix[[x]])<.5){sqrt(mean((c(parallelOutput_sigma_zero2[[x]][[5]], parallelOutput_sigma_zero2[[x]][[6]])- node_sociabilities)^2))}
  else{sqrt(mean((1-c(parallelOutput_sigma_zero2[[x]][[5]], parallelOutput_sigma_zero2[[x]][[6]])- node_sociabilities)^2))}
}

correct_ests = which(unlist(lapply(1:240, function(x){parallelOutput_sigma_zero2[[x]][[7]] == HForms[x] })) )



HAccuracy = rep(0, 12)
MSEs = rep(0, 12)
for(setting in 1:12){
  relevant_runs = ((20*(setting-1))+1):(20*setting)
  HAccuracy[setting] = mean(unlist(lapply(relevant_runs, function(x){parallelOutput_sigma_zero2[[x]][[7]]==HForms[[x]]})))
  MSEs[setting] = mean(unlist(lapply(relevant_runs , check_result_zero)))
  
}


firstcolumn = c(mean(HAccuracy),
                mean(HAccuracy[c(1, 4, 7, 10)]),
                mean(HAccuracy[c(1, 4, 7, 10)+1]),
                mean(HAccuracy[c(1, 4, 7, 10)+2]),
                mean(HAccuracy[1:3]),
                mean(HAccuracy[4:6]),
                mean(HAccuracy[7:9]),
                mean(HAccuracy[10:12])
)


second_column = c(mean(MSEs),
                  mean(MSEs[c(1, 4, 7, 10)]),
                  mean(MSEs[c(1, 4, 7, 10)+1]),
                  mean(MSEs[c(1, 4, 7, 10)+2]),
                  mean(MSEs[1:3]),
                  mean(MSEs[4:6]),
                  mean(MSEs[7:9]),
                  mean(MSEs[10:12])
) 

third_column = c(mean(unlist(lapply( correct_ests , check_result_zero))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)) , check_result_zero))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+20) , check_result_zero))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+40) , check_result_zero))),
                 mean(unlist(lapply( intersect(correct_ests, 1:60) , check_result_zero))),
                 mean(unlist(lapply( intersect(correct_ests, 61:120) , check_result_zero))),
                 mean(unlist(lapply( intersect(correct_ests, 121:180) , check_result_zero))),
                 mean(unlist(lapply( intersect(correct_ests, 181:240) , check_result_zero)))
)


fourth_column = c( mean(unlist(lapply(setdiff(1:240, correct_ests) , check_result_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)) , check_result_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+20) , check_result_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+40) , check_result_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 1:60) , check_result_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 61:120) , check_result_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 121:180) , check_result_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 181:240) , check_result_zero)))
)


fifth_column = c( mean(unlist(lapply( 1:240 , check_result_psi_zero))),
                  mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200) , check_result_psi_zero))),
                  mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200)+20 , check_result_psi_zero))),
                  mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200)+40 , check_result_psi_zero))),
                  mean(unlist(lapply( 1:60 , check_result_psi_zero))),
                  mean(unlist(lapply( 61:120 , check_result_psi_zero))),
                  mean(unlist(lapply( 121:180 , check_result_psi_zero))),
                  mean(unlist(lapply( 181:240 , check_result_psi_zero)))
)


sixth_column = c( mean(unlist(lapply( correct_ests , check_result_psi_zero))), 
                  mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)) , check_result_psi_zero))),
                  mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+20) , check_result_psi_zero))),
                  mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+40) , check_result_psi_zero))),
                  mean(unlist(lapply( intersect(correct_ests, 1:60) , check_result_psi_zero))),
                  mean(unlist(lapply( intersect(correct_ests, 61:120) , check_result_psi_zero))),
                  mean(unlist(lapply( intersect(correct_ests, 121:180) , check_result_psi_zero))),
                  mean(unlist(lapply( intersect(correct_ests, 181:240) , check_result_psi_zero)))
)


seventh_column = c(mean(unlist(lapply(setdiff(1:240, correct_ests) , check_result_psi_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)) , check_result_psi_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+20) , check_result_psi_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+40) , check_result_psi_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 1:60) , check_result_psi_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 61:120) , check_result_psi_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 121:180) , check_result_psi_zero))),
                   mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 181:240) , check_result_psi_zero)))
)




plot(unlist(lapply(1:240 , check_result_zero)),  col = rep(1:4, each=60), pch = (4*unlist(lapply(1:240, function(x){parallelOutput_sigma_zero[[x]][[7]] == HForms[x] }))), cex= .7,  main =paste("Estimation performance for subnetworks \n when σ is estimated to be 0"), ylab ="RMSE"  , xaxt = "n", xlab ="", xlim = c(5, 235))

axis(1, at = seq(30, 210, by = 60), labels = c("Normal", "Linear", "Concave", "Convex"), mgp=c(0,.2,0))
axis(1, at = seq(10, 230, by = 20), labels = c("\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3"), lwd.tick = 0, mgp =c(0, 2,0)  )
abline(v = c(seq(20, 220, by=20)), lty=2, lwd = .5)


legend("topleft",  c("Wrong H",  "Correct H"), pch = c(0, 4), cex=.7)

plot(unlist(lapply(1:240 , check_result_psi_zero)),  col = rep(1:4, each=60), pch = (4*unlist(lapply(1:240, function(x){parallelOutput_sigma_zero2[[x]][[7]] == HForms[x] }))), cex= .7,  main =paste("Psi Estimation MSE for subnetworks \n when σ is estimated to be 0"), ylab ="RMSE"  , xaxt = "n", xlab ="", xlim = c(5, 235))

axis(1, at = seq(30, 210, by = 60), labels = c("Normal", "Linear", "Concave", "Convex"), mgp=c(0,.2,0))
axis(1, at = seq(10, 230, by = 20), labels = c("\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3"), lwd.tick = 0, mgp =c(0, 2,0)  )
abline(v = c(seq(20, 220, by=20)), lty=2, lwd = .5)


legend("topleft",  c("Wrong H",  "Correct H"), pch = c(0, 4), cex=.5)






rm(list=ls())


load("DOROTHY_Acronym200Beta.RData")

check_result = function(x){
  if(mean(all_adjacency_matrix[[x]])<.5){sqrt(mean((parallelOutput2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
  else{sqrt(mean((1-parallelOutput2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
}


check_result_psi = function(x){
  if(mean(all_adjacency_matrix[[x]])<.5){sqrt(mean((c(parallelOutput2[[x]][[5]], parallelOutput2[[x]][[6]])- node_sociabilities)^2))}
  else{sqrt(mean((1-c(parallelOutput2[[x]][[5]], parallelOutput2[[x]][[6]])- node_sociabilities)^2))}
}


mean(unlist(lapply(1:240, function(x){parallelOutput2[[x]][[7]]==HForms[[x]]})))



correct_ests = which(unlist(lapply(1:240, function(x){parallelOutput2[[x]][[7]] == HForms[x] })) )



 HAccuracy = rep(0, 12)
 MSEs = rep(0, 12)
 for(setting in 1:12){
   relevant_runs = ((20*(setting-1))+1):(20*setting)
   HAccuracy[setting] = mean(unlist(lapply(relevant_runs, function(x){parallelOutput2[[x]][[7]]==HForms[[x]]})))
   MSEs[setting] = mean(unlist(lapply(relevant_runs , check_result)))
 
 }

 
 firstcolumn = c(mean(HAccuracy),
 mean(HAccuracy[c(1, 4, 7, 10)]),
 mean(HAccuracy[c(1, 4, 7, 10)+1]),
 mean(HAccuracy[c(1, 4, 7, 10)+2]),
 mean(HAccuracy[1:3]),
 mean(HAccuracy[4:6]),
 mean(HAccuracy[7:9]),
 mean(HAccuracy[10:12])
 )
 
 
 second_column = c(mean(MSEs),
 mean(MSEs[c(1, 4, 7, 10)]),
 mean(MSEs[c(1, 4, 7, 10)+1]),
 mean(MSEs[c(1, 4, 7, 10)+2]),
 mean(MSEs[1:3]),
 mean(MSEs[4:6]),
 mean(MSEs[7:9]),
 mean(MSEs[10:12])
 ) 
 
third_column = c(mean(unlist(lapply( correct_ests , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+20) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+40) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, 1:60) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, 61:120) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, 121:180) , check_result))),
                 mean(unlist(lapply( intersect(correct_ests, 181:240) , check_result)))
)
                 
                 
 fourth_column = c( mean(unlist(lapply(setdiff(1:240, correct_ests) , check_result))),
 mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)) , check_result))),
 mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+20) , check_result))),
 mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+40) , check_result))),
 mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 1:60) , check_result))),
 mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 61:120) , check_result))),
 mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 121:180) , check_result))),
 mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 181:240) , check_result)))
 )
 
 
 fifth_column = c( mean(unlist(lapply( 1:240 , check_result_psi))),
 mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200) , check_result_psi))),
 mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200)+20 , check_result_psi))),
 mean(unlist(lapply( c(1:20, 61:80, 121:140, 181:200)+40 , check_result_psi))),
 mean(unlist(lapply( 1:60 , check_result_psi))),
 mean(unlist(lapply( 61:120 , check_result_psi))),
 mean(unlist(lapply( 121:180 , check_result_psi))),
 mean(unlist(lapply( 181:240 , check_result_psi)))
 )
 
 
 sixth_column = c( mean(unlist(lapply( correct_ests , check_result_psi))), 
 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)) , check_result_psi))),
 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+20) , check_result_psi))),
 mean(unlist(lapply( intersect(correct_ests, c(1:20, 61:80, 121:140, 181:200)+40) , check_result_psi))),
 mean(unlist(lapply( intersect(correct_ests, 1:60) , check_result_psi))),
 mean(unlist(lapply( intersect(correct_ests, 61:120) , check_result_psi))),
 mean(unlist(lapply( intersect(correct_ests, 121:180) , check_result_psi))),
 mean(unlist(lapply( intersect(correct_ests, 181:240) , check_result_psi)))
 )
 
 
 seventh_column = c(mean(unlist(lapply(setdiff(1:240, correct_ests) , check_result_psi))),
  mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)) , check_result_psi))),
  mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+20) , check_result_psi))),
  mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), c(1:20, 61:80, 121:140, 181:200)+40) , check_result_psi))),
  mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 1:60) , check_result_psi))),
  mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 61:120) , check_result_psi))),
 mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 121:180) , check_result_psi))),
 mean(unlist(lapply(intersect(setdiff(1:240, correct_ests), 181:240) , check_result_psi)))
 )
 
 
 
 
 #comparison to DCBM/PABM
 DCBM_symmetric_subnetwork_estimate = function(A)
 {
   #B <- matrix(0, K, K)
   P.hat <- 0*A
   i=1 
   j=2
   N.i <- nrow(A)
   N.j <- ncol(A)
   A_prime = A
   B  <- sum(A_prime, na.rm = TRUE) + 
     0.001
   Psi <- rowSums(A_prime, na.rm = TRUE)
   Psi <- Psi/sum(Psi)
   Psi2 <- colSums(A_prime, na.rm = TRUE)
   Psi2 <- Psi2/sum(Psi2)
   P.hat = B*Psi%*%t(Psi2)
 }
 
 
 
 MSE_DCBM = rep(0, length(all_adjacency_matrix))
 for(netnum in 1:length(all_adjacency_matrix)){
   
   MSE_DCBM[[netnum]] = sqrt(mean((DCBM_symmetric_subnetwork_estimate(all_adjacency_matrix[[netnum]]) - all_no_error_probs_matrix[[netnum]])^2))
   
 }
 
 
 mean(MSE_DCBM)
 
 
 
 
 DCBM_column = c(mean(MSE_DCBM),
                 mean(MSE_DCBM[c(1:20, 61:80, 121:140, 181:200)]),
                 mean(MSE_DCBM[c(1:20, 61:80, 121:140, 181:200)+20]),
                 mean(MSE_DCBM[c(1:20, 61:80, 121:140, 181:200)+40]),
                 mean(MSE_DCBM[1:60]),
                 mean(MSE_DCBM[61:120]),
                 mean(MSE_DCBM[121:180]),
                 mean(MSE_DCBM[181:240])
 ) 
 
 
 
 
 
 
 
 # plot(HAccuracy, MSEs, col = rep(1:4, each=3), pch = as.character(rep(1:3, each=4)))
# 
# plot(unlist(lapply(1:240 , check_result)), col = rep(1:4, each=60), pch = as.character(rep(rep(1:3, each=20),4)))
# 
# plot(unlist(lapply(1:240 , check_result)), col = rep(1:4, each=60), pch = (4*unlist(lapply(1:240, function(x){parallelOutput2[[x]][[7]] == HForms[x] }))), cex= .7, xlab= "Regime", ylab ="MSE"  )
# abline(v = c(seq(20, 220, by=20)), lty=2, lwd = .5)
# 
# 




# 
# mean(unlist(lapply(setdiff(1:240,badSeeds), function(x){parallelOutput_sigma_zero2[[x]][[7]]==HForms[[x]]})))
# 
# mean(unlist(lapply(setdiff(1:240,badSeeds), function(x){parallelOutput_sigma_zero2[[x]][[7]]=='normal'})))
# 

plot(unlist(lapply(1:240 , check_result)),  col = rep(1:4, each=60), pch = (4*unlist(lapply(1:240, function(x){parallelOutput2[[x]][[7]] == HForms[x] }))), cex= .7,  main =paste("Estimation performance for subnetworks where Ψ \n values are drawn from Beta distributions"), ylab ="RMSE"  , xaxt = "n", xlab ="", xlim = c(5, 235))

axis(1, at = seq(30, 210, by = 60), labels = c("Normal", "Linear", "Concave", "Convex"), mgp=c(0,.2,0))
axis(1, at = seq(10, 230, by = 20), labels = c("\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3"), lwd.tick = 0, mgp =c(0, 2,0)  )
abline(v = c(seq(20, 220, by=20)), lty=2, lwd = .5)


legend("topleft",  c("Wrong H",  "Correct H"), pch = c(0, 4), cex=.7)

plot(unlist(lapply(1:240 , check_result_psi)),  col = rep(1:4, each=60), pch = (4*unlist(lapply(1:240, function(x){parallelOutput2[[x]][[7]] == HForms[x] }))), cex= .7,  main =paste("Psi Estimation MSE for subnetworks where Ψ \n values are drawn from Beta distributions"), ylab ="RMSE"  , xaxt = "n", xlab ="", xlim = c(5, 235))

axis(1, at = seq(30, 210, by = 60), labels = c("Normal", "Linear", "Concave", "Convex"), mgp=c(0,.2,0))
axis(1, at = seq(10, 230, by = 20), labels = c("\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3", 
                                               "\n \n State \n 1", "\n \n State \n 2", "\n \n State \n 3"), lwd.tick = 0, mgp =c(0, 2,0)  )
abline(v = c(seq(20, 220, by=20)), lty=2, lwd = .5)


legend("bottomright",  c("Wrong H",  "Correct H"), pch = c(0, 4), cex=.4)












check_result_sigma_zero = function(x){
  if(mean(all_adjacency_matrix[[x]])<.5){sqrt(mean((parallelOutput_sigma_zero2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
  else{sqrt(mean((1-parallelOutput_sigma_zero2[[x]][[9]]- all_no_error_probs_matrix[[x]])^2))}
}



mean(unlist(lapply(setdiff(1:240, badSeeds) , check_result_sigma_zero)))
mean(unlist(lapply(setdiff(1:240, badSeeds) , check_result)))





#Likelihood
mean(unlist(lapply(1:240, function(x){parallelOutput2[[x]][[8]] })) )
mean(unlist(lapply(1:240, function(x){parallelOutputNORMAL2[[x]][[8]] })) )










#t.test(unlist(lapply(setdiff(1:240,correct_ests) , check_result)), unlist(lapply(correct_ests , check_result)))


#extractVals = function(x){unlist(lapply(1:length(parallelOutput2), function(z){parallelOutput2[[z]][[x]]}   ))}
#Reporting: true H, true alpha, true beta, true sigma (.5), true rho (1),  estimated alpha, estimated beta, estimated sigma, estimated rho, estimated H, negative-log-likelihood,  MSE,  time in minutes

#formattedOutoput = data.frame("true H" = HForms, "true alpha" = alphaForms, "true beta" = betaForms, "true sigma" = .5, "true rho" = 1, "estimated alpha" = extractVals(1), "estimated beta"= extractVals(2), "estimated sigma"=extractVals(3), "estimated rho" = extractVals(4), "estimated H" = extractVals(7), "negative log likelihood" = extractVals(8), "MSE" =unlist(lapply(1:240 , check_result)), "time in minutes" =extractVals(12)/60  )
#colnames(output) = c("true H", "true alpha", , "true sigma", "true rho",  "estimated alpha", "estimated beta", "estimated sigma", "estimated rho", "estimated H", "negative log likelihood",  "MSE",  "time in minutes" )






filled.contour(all_no_error_probs_matrix[[127]], color.palette = colorRampPalette(c("purple","blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05), plot.axes = 0)
filled.contour(parallelOutput2[[127]][[9]], color.palette = colorRampPalette(c("purple","blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05), plot.axes = 0)
filled.contour(DCBM_symmetric_subnetwork_estimate(all_adjacency_matrix[[127]]), color.palette = colorRampPalette(c("purple","blue", "green", "yellow", "orange", "red")), levels = seq(0, 1, by = .05), plot.axes = 0)
image(all_adjacency_matrix[[127]], col=c("white", "black"))

