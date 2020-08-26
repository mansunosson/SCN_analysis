% MCMC for CAR model 
% Input:
%
% X   : s by s by T array, time series data on grid with s^2 locations
% dt  : discretization grid for states
% samp: sampling grid
% B   : Length of chain
% a   : prior mean/lower bound for [rho(lb) tau(mu1 mu2 mu3 mu4 mu5) lnTheta(mu1 mu2 mu3 mu4 mu5) Psi(mu1 mu2 mu3)]
% b   : prior sd/upper bound for   [rho(ub) tau(sd1 sd2 sd3 sd4 sd5)                              Psi(sd1 sd2 sd3)]
function [tau, thetamean, RE, theta, psi, tauprop, REprop, thetameanprop,psiprop, LIK0, tauaccept, REaccept, thetameanaccept] = carMCMC(X,dt,samp,B,a,b, thetainit, REinit, psiinit, tauinit, k0, OUT)
    
    
    %Initialize
    height = size(X,1);
    width  = size(X,2);
    L = height*width;
    T = size(X,3);
    X = reshape(X,L,T);
    P         = [5 4 3 2 1]; %order of updates (Fixed scan Gibbs)
    c1 = 0.02;
    %reshape vector of prior means and variances for psi
    psimean   = reshape(repelem(a(12:14),L)',L,3); %change subscripts if rho is omitted
    psivar    = reshape(repelem(b(12:14),L)',L,3);
    %initialise IF iteration counter is at 1
    if k0==1   
    X = reshape(X,L,T);
    LIK0 = NaN(L,B);
    %Allocate empty arrays for hyperparameters
    %rho = NaN(B,1); 
    %rho1= NaN(1);
    %rhoaccept = zeros(B,1);
    tau = NaN(B,5);
    tau1= NaN(5,1);
    tauaccept = zeros(5,B);
    %empty height by width by B by #param arrays
    thetamean = NaN(5,B);      
    thetamean1= NaN(5,1);  
    thetameanaccept=zeros(B,1);
    RE        = NaN(height+2,width+2,B,5);  %random effect array with NaN-padding 
    %RE1       = NaN(5,1);
    
    REaccept  = zeros(5,B);
    psi       = NaN(L,B,3);
    psi1      = NaN(L,3);    
 
    
     %Quantities for subdivision of random effects
     %temp = reshape(1:D,sqrt(D),sqrt(D));
     %IND = kron(ones(height/sqrt(D),width/sqrt(D)),temp);
    
    
     %Initial proposal log(st.dev) that is tuned acc. to accprob
     %rhoprop       = log(repelem(0.02, B));
     tauprop       = log(repelem(0.01,5,1, B));  
     
     REprop(1,:)   = log(repelem(exp(-5)/(L),1, B));
     REprop(2,:)   = log(repelem(exp(-5)/(L),1, B));
     REprop(3,:)   = log(repelem(exp(-5)/(L),1, B));
     REprop(4,:)   = log(repelem(exp(-7)/(L),1, B));
     REprop(5,:)   = log(repelem(exp(-5)/(L),1, B));
     
     psiprop       = log(repelem(0.04, L,  B, 3));
     thetameanprop = log(repelem(exp(-11), B));
     c             = 0.02; %initial proposal variance adaption
     
    %Initialize chains with draw from priors
        XLIK0 = NaN;
        thetalik0 = NaN(L,1);
        while isnan(XLIK0) | ~isfinite(XLIK0) %while loop ensuring finite likelihood at initialization
        %Draw Hyper-parameters
        %rho(1,1) = unifrnd(b(1),a(1)); %from positive side of prior
        %rho(1,1)  = rhoinit; %initialize close to 0
            for i=1:size(tau,2)
            tau(1,i)  = tauinit(i); 
            end
     
            %Draw parameters|hyper-parameters
            for i=1:5
            thetamean(i,1)  = thetainit(1,1,i);
            RE(2:height+1,2:width+1,1,i)     = REinit(:,:,i);    %This indexing gives NaN-padding for neighbour-function
            end
       
            for i=1:size(psi,3)
             psi(:,1,i) = reshape(psiinit(:,:,i),L,1); 
            end
            
            thetapar = NaN(height,width,5);
            for i=1:5           
            thetapar(:,:,i) = exp(thetainit(:,:,i) + reshape(RE(2:height+1,2:width+1,1,i),height,width));
            end
      
            psipar   = reshape(exp(psi(:,1,:)),L,3);
            
            %Calculate full data likelihood of initialization          
            thetapar = reshape(thetapar, L,5);
            parfor l=1:L 
            thetalik0(l)= logLikelihoodv2(dt,reshape(X(l,:),1,T),thetapar(l,3),thetapar(l,4),thetapar(l,5),psipar(l,1),thetapar(l,1),thetapar(l,2),24,psipar(l,2),psipar(l,3),samp); 
            end
            XLIK0 = sum(thetalik0);
            thetapar = reshape(thetapar, height,width,5);
        end    
    

    LIK0(:,1)    = thetalik0;
    
    psiprior0 = NaN(L,3);
    for r=1:3
    for l=1:L
    psiprior0(l,r) = Nprior(psi(l,1,r),psimean(l,r),sqrt(psivar(l,r))); %replace this with log(normpdf) and change to stdev for speed
    end
    end
    
    %initialisation if iteration counter is > 1
    else
     k0                  = OUT.k0;
     c                   = OUT.c;
     tau                 = OUT.tau;
     thetamean           = OUT.thetamean;
     thetapar            = OUT.thetapar;
     RE                  = OUT.RE;
     psi                 = OUT.psi;
     tauprop             = OUT.tauprop;
     REprop              = OUT.REprop;
     thetameanprop       = OUT.thetameanprop;
     psiprop             = OUT.psiprop;
     LIK0                = OUT.LIK0;
     tauaccept           = OUT.tauaccept;
     REaccept            = OUT.REaccept;
     thetameanaccept     = OUT.thetameanaccept; 
     psiprior0           = OUT.psiprior0;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%MAIN LOOP BEGINS HERE%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = k0+1:B  %change 2 to current k
        
    psipar        = reshape(exp(psi(:,k-1,:)),L,3);    
    RE(:,:,k,:)   = RE(:,:,k-1,:);  %(A)
    thetamean(:,k)= thetamean(:,k-1);
    
        %%%Loop over theta parameters%%%
        for t = P 
          %Propose tau
          tau1 = tau(k-1,t) +normrnd(0,exp(tauprop(t,k-1))); %tune this proposal variance
          
          %Evaluate hyper-variance likelihood
          taulik0 = CARprior(RE(:,:,k,t), exp(tau(k-1,t)));  
          taulik1 = CARprior(RE(:,:,k,t), exp(tau1));  
          
          accprob = exp(taulik1+Nprior(tau1,a(1+t),b(1+t))-taulik0-Nprior(tau(k-1,t),a(1+t),b(1+t))); %replace Nprior with log(normpdf) and change to stdev for speed
         
          
          %Accept/reject
          if unifrnd(0,1)< accprob && isfinite(accprob)
            tau(k,t) = tau1;
            tauaccept(t,k) = 1;
          else
            tau(k,t) = tau(k-1,t);
          end
            
          %Adapt proposal variance
          %if accprob<0.5 
          if mean(tauaccept(t,max(1,k-49):k))<0.5 
              tauprop(t,k) = tauprop(t,k-1) - c;
          else
              tauprop(t,k) = tauprop(t,k-1) + c;
          end 
          
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%% Block-sampling of Random effects %%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
       
              %%MH step to sample random effects              
              RE1      = RE(:,:,k,t);                  %make sure this picks up recently updated REs by (A)
              REprior0 = CARprior(RE1, exp(tau(k,t))); %This does not have to be recalculated!       
                                                                         
              RE1(2:height+1,2:width+1) = RE1(2:height+1,2:width+1)+normrnd(0,sqrt(exp(REprop(t,k-1))),[height, width]);
              %Evaluate prior for proposed value | neighbours
              REprior1 = CARprior(RE1, exp(tau(k,t))); %this does not allow for partial updates of random effects!!!
                  
              %Transform into actual parameter values
              thetapar(:,:,t) = exp(thetamean(t,k) + RE1(2:height+1,2:width+1));
              
              %Evaluate likelihood for random effects
              LIK1=NaN(L,1); 
              thetapar = reshape(thetapar, L,5);
              parfor l=1:L    %this loop contains broadcast variables, might be very slow with large L
              LIK1(l)   = logLikelihoodv2(dt,reshape(X(l,:),1,T),thetapar(l,3),thetapar(l,4),thetapar(l,5),psipar(l,1),thetapar(l,1),thetapar(l,2),24,psipar(l,2),psipar(l,3),samp);          
              end
              thetapar = reshape(thetapar, height,width,5);
              accprob     = exp(sum(LIK1) + REprior1- sum(LIK0(:,k-1)) - REprior0);
              
              %%accept/reject
              if unifrnd(0,1)<accprob && isfinite(accprob)
               RE(:,:,k,t)   = RE1; 
               LIK0(:,k-1) = LIK1;
               REaccept(t,k) = 1;
              else
               thetapar(:,:,t) = exp(thetamean(t,k) + RE(2:height+1,2:width+1,k,t));   %throw away proposal if rejected
              end  
              
              %%adapt proposal variance
              if mean(REaccept(t,max(1,k-49):k))<0.23
               REprop(t,k) = REprop(t,k-1) - c;
              else
               REprop(t,k) = REprop(t,k-1) + c;
              end
             %Sum-to-zero constraint
              thetamean(t,k)               = thetamean(t,k) + mean(reshape((RE(2:height+1,2:width+1,k,t)),L,1));
              RE(2:height+1,2:width+1,k,t) = RE(2:height+1,2:width+1,k,t) - mean(reshape((RE(2:height+1,2:width+1,k,t)),L,1)); 
        end    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Block-sampling of Theta mean %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %start proposing correlated increments at iteration 2500
            propcov = 0.1;
            Seye =diag([1 60 60 60 35]);
            %if k>2500 %change this for large L?
            %propcov = cov(thetamean(:,k-2400:k-1)');
            %Seye =diag([1 1 1 1 1]);
            %end
            
            thetamean1 = thetamean(:,k) + mvnrnd(zeros(5,1),0.9*exp(thetameanprop(k-1))*propcov*Seye + 0.1*exp(thetameanprop(k-1))*Seye)';          
            
            thetaprior0 = 0;
            thetaprior1 = 0;
            for t=1:5
              if t == 1 | t == 2
                thetaprior0   = thetaprior0 + Uprior(exp(thetamean(t,k)), b(6+t), a(6+t)) + thetamean(t,k); %Jacobians
                thetaprior1   = thetaprior1 + Uprior(exp(thetamean1(t)),  b(6+t), a(6+t)) + thetamean1(t);  %Jacobians
                thetapar(:,:,t) = exp(thetamean1(t) + RE(2:height+1,2:width+1,k,t)); %transform to parameter values
              else
                thetaprior0   = thetaprior0 + Nprior(thetamean(t,k), a(6+t), b(6+t));
                thetaprior1   = thetaprior1 + Nprior(thetamean1(t),  a(6+t), b(6+t)); 
                thetapar(:,:,t) = exp(thetamean1(t) + RE(2:height+1,2:width+1,k,t)); %transform to parameter values
              end    
            end                    
       
            LIK1  = NaN(L,1); %reset proposed thetamean likelihood 
            thetapar = reshape(thetapar,L,5);
            parfor l = 1:L    
                LIK1(l)  = logLikelihoodv2(dt,reshape(X(l,:),1,T),thetapar(l,3),thetapar(l,4),thetapar(l,5),psipar(l,1),thetapar(l,1),thetapar(l,2),24,psipar(l,2),psipar(l,3),samp); 
            end
            thetapar = reshape(thetapar,height,width,5);
            
                accprob    = exp(sum(LIK1) + thetaprior1- sum(LIK0(:,k-1)) - thetaprior0);               
                %accept/reject
                if unifrnd(0,1)<accprob && isfinite(accprob)
                 thetamean(:,k) = thetamean1;
                 LIK0(:,k-1) = LIK1;
                 thetameanaccept(k) = 1; 
                else
                 for t=1:5
                  thetapar(:,:,t) = exp(thetamean(t,k) + RE(2:height+1,2:width+1,k,t)); %throw away proposal if rejected
                 end    
                end
                
                %adapt proposal variance
                if mean(thetameanaccept(max(1,k-49):k))<0.23
                   thetameanprop(k) = thetameanprop(k-1) - c;
                else
                   thetameanprop(k) = thetameanprop(k-1) + c;
                end
        
        
             
          psi(:,k,:) =  psi(:,k-1,:); 
          thetapar = reshape(thetapar,L,5);
          %Block proposal of Psi-parameters
          for l = 1:L        
              psi1(l,:)   = reshape(psi(l,k,:),3,1);
              psi1(l,:)   = psi1(l,:)+mvnrnd([0 0 0],diag(exp(reshape(psiprop(l,k-1,:),3,1)).^2)); 
              psipar(l,:) = exp(psi1(l,:));
          end
              LIK1 = NaN(L,1);
              parfor l = 1:L
              LIK1(l)   = logLikelihoodv2(dt,reshape(X(l,:),1,T),thetapar(l,3),thetapar(l,4),thetapar(l,5),psipar(l,1),thetapar(l,1),thetapar(l,2),24,psipar(l,2),psipar(l,3),samp); 
              end
 
          for l = 1:L  
              
              psiprior1(1)  = Nprior(psi1(l,1),psimean(l,1),sqrt(psivar(l,1)));
              psiprior1(2)  = Nprior(psi1(l,2),psimean(l,2),sqrt(psivar(l,2))); 
              psiprior1(3)  = Nprior(psi1(l,3),psimean(l,3),sqrt(psivar(l,3)));
              
              accprob     = exp(LIK1(l) + sum(psiprior1) -LIK0(l,k-1)-sum(psiprior0(l,:)));          
              %adapt proposal variance
              if accprob<0.36
                   psiprop(l,k,:) = psiprop(l,k-1,:) - c;
              else
                   psiprop(l,k,:) = psiprop(l,k-1,:) + c;
              end
              %accept/reject
              if unifrnd(0,1)<accprob && isfinite(accprob)
                psi(l,k,:)     = psi1(l,:);
                psiprior0(l,:) = psiprior1;
                LIK0(l,k-1)   = LIK1(l);                      
              end
          end
        thetapar = reshape(thetapar, height, width, 5);           
        %print iteration
        k
        
        %set current likelihood to most recent
        LIK0(:,k) = LIK0(:,k-1);
        %diminishing adaption
        if k<25000
        c=c-c1/(25000); 
        end
        if k==25001
         c=0
        end
        
      
     %Save output every 500 iterations   
     if mod(k,200)==0   
      k0=k;
     save('OUT.mat', 'k0', 'c', 'tau', 'thetamean', 'thetapar', 'RE', 'psi', 'tauprop', 'REprop', 'thetameanprop', 'psiprop', 'LIK0', 'tauaccept', 'REaccept', 'thetameanaccept', 'psiprior0', '-v7.3')
     end  
     
    %end main loop   
    end
    
    %transform parameters
    theta = NaN(height,width,B,5);
    for t=P
    for i=1:B
    theta(:,:,i,t) = exp(thetamean(t,i) + RE(2:height+1,2:width+1,i,t));
    end
    end
    psi = reshape(psi,height,width,B,3);
    LIK0 = reshape(LIK0,height,width,B); 
    

end