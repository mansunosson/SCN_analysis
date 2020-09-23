# SCN_analysis
Matlab code for analyzing circadian SCN imaging data

run_code_rep1.m : These scripts contain data pre-processing for each of the experimental replicates and calls the function CARmcmc.m 
run_code_rep2.m : lines 58-101 make sure that the initial values are in dense regions of the posterior. This is not required for the analysis 
run_code_rep3.m : but reduces the required number of iterations by shortening the time the chains spend in a transient phase.


CARmcmc.m is the main function that calls:
    CARprior.m:        Takes as input a set of random effects (RE) and variance (tau) and evaluates the joint pdf of RE as
                         specified by Eq. (4) in the main text. Calls function neighbour.m to 
             
    logLikelihoodv2.m: Takes as input parameters of the single-cell TTFL model in Eq. (1) and measurement model in Eq. (2), 
                         along with observed data and calculates likelihood by calling EKBFv2.m.
                         
    Nprior.m:          Evaluates the log normal pdf for "par" with mean "a" and SD "b".
    Uprior.m:          Evaluates the uniform pdf for "par" with upper bound "b" and lower bound "a".
    multNprior.m:      Evaluates the log multivariate normal pdf for "par" with mean "a" and covariance "b".

PPLC.mat is used to calculate the posterior probability of a limit cycle.
SDEsolve.mat is used to simulate ensembles of the model in Eq. (1).




