%%% CAR prior %%%
% Returns log of normal prior for random effects and hyper-parameters
% Calls parfor, neighbours-function and logmvnpdf-function
%
% Input:
%       RE         : Matrix of random effects (padded)
%       tau        : hyper-variance
%       rho        : optional parameter (REMOVED) for strength 
%                    of spatial correlation that guarantees
%                    propriety if -1<rho<1

function [prior] = CARprior(RE, tau)
  nrow = size(RE,1);
  ncol = size(RE,2);
  [rows, columns] = ndgrid(2:nrow-1, 2:ncol-1);
  L = size(rows,1)*size(rows,2);
  indices = reshape(sub2ind([nrow ncol], rows, columns),L,1); %indices for non-padding
  priors  = NaN(L,1);
    parfor i=1:L %size(indices,1) changed to L
    NB        = neighbours(indices(i), RE);
    nn        = length(NB(~isnan(NB)));  %CAR-scaling of variance
    a         = nanmean(NB);             %Prior mean
    priors(i) = normpdf(RE(indices(i)),a,sqrt(tau/nn)); 
    end
  prior = sum(log(priors));
end
