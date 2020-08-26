function [prior] = Nprior(par,a,b)
prior= logmvnpdf(par,a,b^2);
end
