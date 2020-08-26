%%% Simulate spatial distribution of parameters and data %%%
function [X, theta, RE, psi, W, Wtil] = simCAR(T, dt, s, rho, tau, thetamean, psimean, psivar, connection)
          %T  : number of frames
          %s  : grid size (s by s)
          %rho: propriety parameter, |rho|<1 that ensures existance
          %of joint distribution of theta.

%Generate index matrix for neighbourhood structure
L = 1:s*s;
P = reshape(L,s,s);

if connection == 4
%4-connected neigbourhood adjacency matrix                         
diagVec1 = repmat([ones(s-1, 1); 0], s, 1);                                            
diagVec1 = diagVec1(1:end-1);               
diagVec2 = ones(s*(s-1), 1);                                                         
adj = diag(diagVec1, 1)+diag(diagVec2, s);   
W = adj+adj.';                            
end

if connection == 8
%8-connected neigbourhood adjacency matrix
diagVec1 = repmat([ones(s-1, 1); 0], s, 1);                                   
diagVec1 = diagVec1(1:end-1);               
diagVec2 = [0; diagVec1(1:(s*(s-1)))];                                                
diagVec3 = ones(s*(s-1), 1);                                                          
diagVec4 = diagVec2(2:end-1);                                                         
adj = diag(diagVec1, 1)+...                 
      diag(diagVec2, s-1)+...
      diag(diagVec3, s)+...
      diag(diagVec4, s+1);
W = adj+adj.';                           
end

%Construct Sigma_theta
invDw = inv(diag(sum(W)));
Wtil  = invDw*W;
S  = NaN(s*s,s*s,length(thetamean));
for i = 1:length(thetamean)
S(:,:,i) = inv(eye(s*s)-rho*Wtil)*tau(i)*invDw; %Covariance exists for |rho|<1
end
%Draw theta parameters
theta = NaN(s,s,length(thetamean)); %allocate empty array for theta

for i=1:length(thetamean)
temp = mvnrnd(repelem(0,size(S,1)),S(:,:,i));
temp = temp -mean(temp); %sum to zero constraint
RE(:,:,i) = reshape(temp, s,s);
end

for i=1:length(thetamean)
theta(:,:,i) = thetamean(i) + RE(:,:,i);
end

%Draw psi parameters
psi = NaN(s,s, length(psimean));
for i=1:length(psimean)
psi(:,:,i) = normrnd(psimean(i),sqrt(psivar(i)),s,s);
end

%Simulate data
X = NaN(s,s,T/dt); %allocate array for simulated data



for i=1:s
for j=1:s
aP = exp(theta(i,j,1)).^2./exp(theta(i,j,2)).^2;
bP = exp(theta(i,j,1))./exp(theta(i,j,2)).^2;
X(i,j,:) = exp(psi(i,j,2)).*SDESolve(T, dt,aP, bP, exp(theta(i,j,3)),exp(theta(i,j,4)),exp(theta(i,j,5)), exp(psi(i,j,1)), repelem(100,24/dt)) + normrnd(0,exp(psi(i,j,3)),T/dt,1);
end
end

end