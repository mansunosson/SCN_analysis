%Implementation of Lugsail Batch Means Gelman-Rubin diagnostic from 
% Vats & Knudson (2018)
% INPUT:
% chain   :   MCMC chain (n by m)
% a       :   Number of batches
% b       :   Batch size
% note: a, b have to be specified s.t. a*b=length(chain). Typical choices for b are floor(n^1/2) 
% and floor(n^1/3) 
%
% OUTPUT  :   R, if R<1+epsilon simulation can be terminated.
function [R, ESS] = LBMGR(chain, a, b)
n  = size(chain,1);
m = size(chain,2);
s2 = mean(var(chain));

batchmeans = mean(reshape(chain,b,a));
tau = (b/(a-1))*sum((batchmeans - mean(chain)).^2);


b3 = b/3;
a3 = a*3;

batchmeans3 = mean(reshape(chain,b3,a3));
tau3= (b3/(a3-1))*sum((batchmeans - mean(chain)).^2);


tauL = 2*tau-tau3;
tauL = mean(tauL);
R = sqrt((((n-1)/n)*s2+tauL/n)/s2);
ESS = n*s2/tauL;
end

