function [ths,thp,meanout,varout,Pxxs,Pxxp, feedback] = EKBFv2(dt,data,R0,Kp,np,degr,meanG,sdG,maxd,kappa,sigmae,samp)

%dt is the discretization time step (hours);
%m is the number of observations to predict before relinearization of unobserved;
%data are the observations;
%R0,Kp,np are parameters of the Hill function;
%degr is the degradation rate;
%meanG and sdG are the mean and standard deviation of the delay distribution;
%maxd is the maximum delay time and length of the initial condition (hours);
%kappa is a scaling factor for the noise;
%sigmae is the measurement error st.dev;
%samp is sampling time interval (in hours)



n=length(data);
integr=samp/dt; %number of states gridpoints in observation
delaydt=maxd/dt; %number of states gridpoints in delay

initY = data(1:(maxd/samp));
initX = initY./kappa;
%initX = interp1(samp:samp:maxd,initY,0:dt:maxd-dt)./kappa; %interp1 is
%very expensive, only use this if discretization of process is finer than data availability
%initX(1:integr) = initX(integr+1);

y=data;

meanout=zeros(1,n); %predicted mean of the observed states y
meanout(1:delaydt/integr) = initY;
varout=zeros(1,1,n); %predicted variance of the observed states y

%matrices for the unobserved states x
PxxI=zeros(integr*n,integr*n); %current variance/covariance matrix
Pxxs=zeros(integr*n,integr*n); %updated variance/covariance matrix
Pxxp=zeros(1,integr*n); %predicted variance estimate
ths=zeros(1,integr*n); %updated mean estimate
thp=zeros(1,integr*n); %predicted mean estimate
feedback=zeros(1,integr*n); %filtered feedback

thp(1:delaydt)=initX;
ths(1:delaydt)=initX; %?

%degradation quantities
A=1-degr*dt;

%Create measurement matrix
FF = eye(maxd/samp);
vecf = ones(1,integr)*dt./samp;
FF = kappa.*kron(FF,vecf);

vv=delaydt; 
M=delaydt; %this is a counter for the number of unbserved states currently predicted
ycount=maxd/samp; %this is a counter for the number of observations currently predicted

t=0;
%%obtain weights for the delay
timegrid=0:dt:(maxd-dt);
timegridH=(timegrid+timegrid+dt)./2;

%get delay parameters for the gamma density
aP=(meanG./sdG).^2;
bP=meanG./sdG.^2;

evalp    =maxd-timegridH; %points at which the density is evaluated
weightsM =bP.^aP./gamma(aP)*(evalp).^(aP-1).*exp(-bP.*(evalp));
weightsMN=(weightsM./sum(weightsM))';

for s=((maxd/samp)+1):n
    th=ths;
    thI=th;
    PxxI(vv-delaydt+1:vv,vv-delaydt+1:vv)=Pxxs(vv-delaydt+1:vv,vv-delaydt+1:vv);
    
for j=M:(M+(integr)-1) %loop until the next observation we want to predict
    
       %evaluate weighted transcription function from 0 to tau_max 
        thInp=thI((j-delaydt+1):j);
           
        PFI=thInp*weightsMN; 
        trfun=R0./(1+(PFI/Kp)^np);
        trfunDerC=-(R0.*np.*(PFI./Kp).^(np-1))./(Kp.*(1+(PFI/Kp)^np).^2);

        %predict mean of process
        thI(j+1)=A*thI(j)+dt.*(trfun); 
        
       %predict variance of process 
       HH=(1-2*degr*dt)*PxxI(vv,vv);
       covmV=dt.*PxxI(vv,vv-delaydt+(1:delaydt))*(trfunDerC.*weightsMN);
       covCC=PxxI(vv-delaydt+(1:delaydt),vv-delaydt+(1:delaydt))*(trfunDerC.*weightsMN);
       DD=dt.*(degr.*thI(j)+trfun); %%diffusion matrix
       PxxI(vv+1,vv+1)=HH+2*sum(covmV)+DD; %%variance
    
       %predict covariance of process
       PxxI(vv-delaydt+(1:delaydt),vv+1)=PxxI(vv-delaydt+(1:delaydt),vv)*A'+dt.*covCC; 
       PxxI(vv+1,vv-delaydt+(1:delaydt))=PxxI(vv-delaydt+(1:delaydt),vv+1)';

       feedback(j+1) = (PFI/Kp)^np;
    
    t=t+1;
    vv=vv+1;
end

     M=M+integr;
     
    %covariance of process and observations
    Pxy=PxxI(vv-delaydt+1:vv,vv-delaydt+1:vv)*(FF(end,:))';
     
    %variance of observations
    Pyy=(FF(end,:))*PxxI(vv-delaydt+1:vv,vv-delaydt+1:vv)*(FF(end,:))';

     thvec(1:integr)=thI(1,(M+1-integr):M);
     th((vv-integr+1):vv)=thvec';     
     thp((vv-integr+1):vv)=thvec';
     Pxxp((vv-integr+1):vv)=diag(PxxI((vv-integr+1):vv,(vv-integr+1):vv));

     
     %calculate E[Yt|X_{t-1}]
     meanout(:,s)=FF(1,1:integr)*thvec';     
     %calculate V[Yt|X_{t-1}]
     varout(:,:,s)=Pyy+sigmae.^2;
     
    %%Kalman update
    ycount=ycount+1;
    ths(vv-delaydt+1:vv)=th(vv-delaydt+1:vv)'+Pxy/(Pyy+sigmae.^2)*(y(ycount)'-FF(1,1:integr)*thvec');
    Pxxs(vv-delaydt+1:vv,vv-delaydt+1:vv)=PxxI(vv-delaydt+1:vv,vv-delaydt+1:vv)-Pxy/(Pyy+sigmae.^2)*Pxy'; 

       
end

end
