%%% SDEsolve %%%
%Numerical solutions to CLE using Euler-Maruyama scheme
%Inputs:
%T:  Solution endpoint, solution lies on [0, T]
%dt: Discretization interval
%aP,...,deg: Parameters of CLE
%hist: initial data, column vector containing X(t), t=dt,...,24
%dWmean: allow for non-zero mean noise process
function X = SDESolve(T, dt,aP, bP, R0, K, n, deg, hist)
timegrid = dt:dt:24;
weightsM =bP^aP/gamma(aP)*(timegrid).^(aP-1).*exp(-bP*(timegrid));
weightsMN=(weightsM./sum(weightsM));

N=T/dt;
ndel=24/dt; %length of delay 
X=NaN(N,1); %Pre-allocate X
X(1:length(hist)) = hist;
dW=normrnd(0,1,[N,1]); %Draw innovations
for t=ndel+1:N
%nu    = R0/(1 + (weightsMN*X(t-ndel:t-1)/K).^n);
nu    = R0/(1 + (weightsMN*X(t-1:-1:t-ndel)/K).^n);
mu    = nu - deg.*X(t-1);
sigma = sqrt(nu+deg.*X(t-1));

%add reflective boundary at X=0
X(t) = X(t-1) + mu*dt + sigma*sqrt(dt)*dW(t);
end

end