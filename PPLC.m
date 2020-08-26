%%% Function that calculates Posterior probability of limit cycle from
%%% array of MCMC output

%INPUT
%Chain  : Array of markov chains from CARmcmc
%burnin : Number of iterations to treat as burn in
%n      : Calculate the PPLC of every nth component of chains (thinning)

function [prob, xstar, meanGarray, sdGarray, R0array, Karray, narray, muarray, kappaarray, sigmaearray] = PPLC(chains, burnin, ss)
syms X

height   = size(chains.thetaest,1);
width    = size(chains.thetaest,2);

meanGarray = NaN(height,width,chains.B);
sdGarray   = NaN(height,width,chains.B);
R0array    = NaN(height,width,chains.B);
Karray     = NaN(height,width,chains.B);
narray     = NaN(height,width,chains.B);
muarray    = NaN(height,width,chains.B);
kappaarray = NaN(height,width,chains.B);
sigmaearray= NaN(height,width,chains.B);

for i = 1:height
for j = 1:width
       meanGarray(i,j,:) = chains.thetaest(i,j, 1:end, 1);
       sdGarray(i,j,:)   = chains.thetaest(i,j, 1:end, 2);
       R0array(i,j,:)    = chains.thetaest(i,j, 1:end, 3);
       Karray(i,j,:)     = chains.thetaest(i,j, 1:end, 4);
       narray(i,j,:)     = chains.thetaest(i,j, 1:end, 5);
       muarray(i,j,:)    = exp(chains.psiest(i,j,1:end,1));
       kappaarray(i,j,:) = exp(chains.psiest(i,j,1:end,2));
       sigmaearray(i,j,:)= exp(chains.psiest(i,j,1:end,3));
end
end 

%initialise estimated probability matrix
N=length(burnin+1:ss:chains.B);
prob   = NaN(height,width, N);
xstar  = NaN(height,width,N);
maxval = NaN(height,width, N);
EV     = cell(height,width, N);
HC     = zeros(height,width, N);

%Loop over locations
for i = 1:height
for j = 1:width
    IND = 1;
    %if location is sampled  
    for u = burnin+1:ss:chains.B
     
       p          = meanGarray(i,j,u)./(sdGarray(i,j,u).^2);
       a          = (meanGarray(i,j,u).^2)./(sdGarray(i,j,u).^2);       
       R0         = R0array(i,j,u);
       K          = Karray(i,j,u);
       n          = narray(i,j,u);
       mu         = muarray(i,j,u);
   
       %Construct Jacobian
       xstar(i,j,IND)      = vpasolve(R0./(1+ (X./K).^n)-mu.*X == 0, X,[0 5000]);
       Fprime     = (-n*R0*(xstar(i,j,IND)/K)^n)./(xstar(i,j,IND)*(1+(xstar(i,j,IND)/K)^n)^2);
       ar= round(a);
       pr= (ar./(a)).*p;
       J = diag(repelem(-pr,ar+1));
       J(1,1) = -mu;
       J(1,ar+1) = Fprime;
       od = [1:ar; 2:ar+1];
       J(2:(ar+2):(ar+1)^2) = pr; 
       if max(real(eig(J)))<0
       prob(i,j,IND) = 0;
       else
       prob(i,j,IND) = 1;  
       end 
       %Solve using characteristic eqn here
       %roots= vpasolve((X+mu)*(X+p)^(a)-Fprime*(-p)^a == 0, X,[4, -2+ 4i]);
       %
       %Store max(re(lambda))
       if -1*det(-J)<0 
           HC(i,j,IND) = 1;
       end
       
       EV{i,j,IND} =  eig(J);
     
       IND = IND+1;
    end  
end 
 i  %print iteration when completing a row    
end

end