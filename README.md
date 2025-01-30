# ACRONYM
ACRONYM: Augmented degree corrected, Community Reticulated Organized Network Yielding Model


This is the code used to for the paper 
"ACRONYM: Augmented degree corrected, Community Reticulated Organized Network Yielding Model"

There are 2 data examples and 3 types of simulations.
Each of these are in their own folder. Datasets are not included
but information about how to retrieve these datasets are included below.


Installing R library dependencies are included in the file "dependencies.R"

"RunEverything.R" runs all the code, but sequentially, so it will take a long time (like a week).
As it is, without the datasets, the simulations will run, but the Congress and Mouse analyses won't. 
Each of the 5 analyses can be run in parallel and parts of each can be run in parallel as well.
If chained together more intentionally using the instructions below, and given enough clusters,
it can all be done in under a day. 

  

Data:

1) Congressional Twitter Dataset 

	Runs ACRONYM on the full dataset and creates 100 partially hidden networks
	and runs ACRONYM on them, too. 	

	Data availability: 
	The .zip file can be retrieved from here
	https://zenodo.org/records/8253486
 
	We use the file 
	congress_network_data.json

	This is read in and formatted at the beginning  
	/Congress/CongressParallelInit.R  

	Running our code:
		- Download congress_network_data.json to /Congress
		- Change workspace to /Congress
		- Run CongressParallelInit.R 
		  This creates common networks and does community detection
		- Once that is complete, you can run all of 
	  	  CongressParallelU.R
	  	  CongressParallelV.R
	  	  CongressParallelEig.R
	  	  In parallel, as they are estimating the same networks using different community detection approaches
		- When those 3 files have run, use CongressPostAnalysis2.R to view results



2) Mouse Retina Dataset 

	Runs ACRONYM on the full dataset and creates 30 partially hidden networks
	and runs ACRONYM on them, too. 	

	Data availability:
	The mouse retina.db file can be retrieved from here
	github.com/ericmjonas/circuitdata/tree/master/mouseretina


	We use the file 
	mouseretina.db


	This is read in and formatted at the beginning  
	/MouseRetina/MouseInit.R  


	Running our code:
		- Download mouseretina.db to /MouseRetina
		- Change workspace to /MouseRetina
		- Run MouseInit.R
		  This creates common networks and does community detection
		- Once that is complete, you can run all of 
		  MouseU.R
		  MouseV.R
		  MouseEig.R
		  In parallel, as they are estimating the same networks using different community detection approaches
		- When those 3 files have run, use MousePostAnalysis2.R to view results



Simulations:

1) Main Simulation. Section 6 in the paper, with more detail in the appendix

	Creates one 1200x1200 network, and 30 partially hidden networks from this network.
	Runs ACRONYM on all of them.  

	Running our code:
		- Change workspace to /MainSim
		- Run MainSimInit.R  
		- This creates the networks we will be modeling, and does community detection
		- Once that is complete, you can run all of 
		  MainSimU.R
		  MainSimV.R
		  MainSimEig.R
		  In parallel, as they are estimating the same networks using different community detection approaches
		- When those 3 files have run, use MainSimPostAnalysis.R to view results


2) Networks with random parameters. See the appendix

	Creates 30 1200x1200 networks with random parameters that kind of mimic what we see in the main simulation.
	Runs ACRONYM on all of them.  

	Running our code:
		- Change workspace to /Random1200
		- Run Redo_RandomAcronym30NestedFull2.R 
		  This creates the networks and runs ACRONYM using Eigenvector community detection
		- Once that is complete, you can run 
		  Redo_SingularValuesTest.R
		  Redo_SingularValuesTestV.R
		  In parallel, as they are estimating the same networks using left and right-singular vectors for community detection 
		- When files have run, use postEstimationAnalysisforRandom1200networks.R to view results







3) Subnetworks. Simulation to remove the effect of community detection
	
	Creates 240 200x200 subnetworks generated from 3 different states, and using 4 different H-functions
	20 replicates of each state/H-function combination
	Runs ACRONYM on all subnetworks
	Also runs ACRONYM with only the ``Normal" H-function available, rather than choosing the H-function based on the data
	Also runs ACRONYM with the estimate of sigma set to 0 (Redo_AcronymParallel200sigma0.R)
	Also creates a new set of 240 200x200 subnetworks where the Psi values are drawn from Beta distributions, 
	and runs ACRONYM on those networks (AcronymParallel200Beta.R)


	Running our code:
		- Change workspace to /Subnetworks
		- Run Redo_AcronymParallel200.R
		- Can run AcronymParallel200Beta.R at the same time
		- Once Redo_AcronymParallel200.R is complete, run Redo_AcronymParallel200sigma0.R
		- When all files have run, use postEstimationAnalysisfor200Subnetwork.R to view results



	

