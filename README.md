# BMDA-Simulation

## A. Permuted Real Data Study

1.  Code: code_real.R and false_discovery_code_real.R includes codes for permuting the group labels of the real data, running the permuted data through each of the competing methods, and calculating the number of false discoveries by each method for each permutation. 
2.  Real Data: contains csv file of real data used in "Recurrent urinary tract infection and estrogen shape the taxonomic ecology and functional potential of the postmenopausal urogenital microbiome" by Neugent et al. (2022).

## B. Simulation Study

Simulation code, data, and results are provided in the files for BMDA versus competing methods. The simulation results are provided in the supplementary materials of "Recurrent urinary tract infection and estrogen shape the taxonomic ecology and functional potential of the postmenopausal urogenital microbiome" by Neugent et al. (2022). 

1.  Code/Code contains functions in R and C++ necessary for the simulation.
2.  Code: false_discovery_code.R and power_fdr_code.R recreate the plots in the paper.
3.  Simulated Data: contains the simulated Weiss data
4.  Simulation Results: contains the results of all methods for both the simulated data and the permuted simulated data. 
    - There is a warning message in these folders "Sorry, we had to truncate this directory to 1,000 files. 1,801 entries were omitted from the list."
    - SOLUTION: If you want to see all the files, you need to clone the repository locally. If that does not work, please email Dr. Qiwei Li at qiwei.li@utdallas.edu and he can send you a link to a Dropbox containing files of all of the results. 
