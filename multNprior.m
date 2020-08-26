function [prior] = multNprior(par,a,b)
prior=logmvnpdf(par,a,b); %par has to be column vector (1 x n)
end
