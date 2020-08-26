function [prior] = Uprior(par,a,b)
if (par<a) || (par>b) 
    prior=-Inf;
else
    prior=1;
end



