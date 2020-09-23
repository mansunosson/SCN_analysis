# SCN_analysis
Matlab code for analyzing circadian SCN imaging data


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




