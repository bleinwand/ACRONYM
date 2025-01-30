source("dependencies.R")

library(foreach)
library(doParallel)
library(mclust)
library(nett)
library(randnet)
library(T4cluster)
library(igraph)
library(viridis)


#Congress
library(rjson)
library(RSpectra)
library(caret)

###MOUSE
library(plyr)
library(RSQLite)
library(RColorBrewer)



current_dir = getwd()
MainSim_dir = paste0(current_dir,"/MainSim")
setwd(MainSim_dir)

source("MainSimInit.R")
rm(list=ls())
source("MainSimU.R")
rm(list=ls())
source("MainSimV.R")
rm(list=ls())
source("MainSimEig.R")
rm(list=ls())
source("MainSimPostAnalysis.R")
save.image(file = "SimPostAnalysis.RData")
rm(list=ls())


setwd("..")
current_dir = getwd()
Random1200_dir = paste0(current_dir,"/Random1200")
setwd(Random1200_dir) 


source("Redo_RandomAcronym30NestedFull2.R")
rm(list=ls())
source("Redo_SingularValuesTest.R")
rm(list=ls())
source("Redo_SingularValuesTestV.R")
rm(list=ls())
source("postEstimationAnalysisforRandom1200networks.R")
save.image(file = "Random1200PostAnalysis.RData")
rm(list=ls())



setwd("..")
current_dir = getwd()
Subnetworks_dir = paste0(current_dir,"/Subnetworks")
setwd(Subnetworks_dir) 

source("Redo_AcronymParallel200.R")
rm(list=ls())
source("Redo_AcronymParallel200sigma0.R")
rm(list=ls())
source("AcronymParallel200Beta.R")
rm(list=ls())
source("postEstimationAnalysisfor200Subnetwork.R")
save.image(file = "SubnetworksPostAnalysis.RData")
rm(list=ls())



setwd("..")
current_dir = getwd()
Congress_dir = paste0(current_dir,"/Congress")
setwd(Congress_dir) 


source("CongressParallelInit.R")
rm(list=ls())
source("CongressParallelU.R")
rm(list=ls())
source("CongressParallelV.R")
rm(list=ls())
source("CongressParallelEig.R")
rm(list=ls())
source("CongressPostAnalysis2.R")
save.image(file = "CongressPostAnalysis.RData")
rm(list=ls())


setwd("..")
current_dir = getwd()
Mouse_dir = paste0(current_dir,"/MouseRetina")
setwd(Mouse_dir) 

source("MouseInit.R")
rm(list=ls())
source("MouseU.R")
rm(list=ls())
source("MouseV.R")
rm(list=ls())
source("MouseEig.R")
rm(list=ls())
source("MousePostAnalysis2.R")
save.image(file = "MousePostAnalysis.RData")



