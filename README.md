# SCN_analysis
Matlab code for analyzing circadian SCN imaging data


CARmcmc.mat is the main function that calls:
    CARprior.mat:        Takes as input a set of random effects (RE) and variance (tau) and evaluates the joint pdf of RE as
                         specified by Eq. (4) in the main text. Calls function neighbour.mat to 
             
    logLikelihoodv2.mat: Takes as input parameters of the single-cell TTFL model in Eq. (1) and measurement model in Eq. (2), 
                         along with observed data and calculates likelihood by calling EKBFv2.mat.
                         
    Nprior.mat:          Evaluates the log normal pdf for "par" with mean "a" and SD "b".
    Uprior.mat:          Evaluates the uniform pdf for "par" with upper bound "b" and lower bound "a".
    multNprior.mat:      Evaluates the log multivariate normal pdf for "par" with mean "a" and covariance "b".

LBMGR.mat is used to calculate lugsail batch means effective sample sizes from a serially correlated MCMC chain.
PPLC.mat is used to calculate the posterior probability of a limit cycle.
SDEsolve.mat is used to simulate ensembles of the model in Eq. (1).




