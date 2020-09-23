clear all; 


%Initialise using parameter estimates from Exp3
load('OUT - exp3 COMBINED.mat')
load('kappaprior.mat')


clear notrend_reCRY1
%% Data import %%
vData = VideoReader('cry1 exp4 composite  crop.avi');
SCN = read(vData);
zdata = im2double(SCN(:,:,:,:));
cry1 = reshape(zdata(:,:,3,1:288), size(zdata,1), size(zdata,2), size(zdata,4)); %Blue light
T  = size(cry1,3);

clear zdata SCN vData 
%%

%Image resizing
resfun = @(block_struct) mean2(block_struct.data); 
for t =1:T
reCRY1(:,:,t) = blockproc(reshape(cry1(:,:,t),  size(cry1,1) , size(cry1,2)), [4 4], resfun);
end

%DETREND DATA and find timechange%
timechange = NaN(size(reCRY1,1),size(reCRY1,2));
notrend_reCRY1 = NaN(size(reCRY1,1),size(reCRY1,2),T);
for i = 1:size(reCRY1,1)
for j = 1:size(reCRY1,2)
    if ~isnan(reCRY1(i,j,1));
    [M I] = max(reCRY1(i,j,1:48));
    timechange(i,j) = I*0.5;
    deg = M/max(reCRY1(i,j,T-48:end));
    
    notrend_reCRY1(i,j,:) = detrend(reshape(reCRY1(i,j,:),T,1)).*linspace(1, deg, T)'+mean(reshape(reCRY1(i,j,:),T,1));
    
    if min(notrend_reCRY1(i,j,:))<0
        notrend_reCRY1(i,j,:) = notrend_reCRY1(i,j,:)-min(notrend_reCRY1(i,j,:));
    end
    
    end
end
end
clear M I deg a cry1 reCRY1 i j resfun timechange yLocsd yLocse xLocsd xLocse


%imagesc(notrend_reCRY1(5:100,1:50,24)) %view data
dshift = 0;
rshift = 0;
height = 90;
width  = 50;

data = notrend_reCRY1(1+dshift:height+dshift,1+rshift:width+rshift,:);
L=height*width;
kappainit = kpriormean(1+dshift:height+dshift,1+rshift:width+rshift);

init = NaN(height,width,8);
for i=1:height
for j=1:width
    if size(PAR1{i+dshift,j+rshift})~= 0
      par1 = PAR1{i+dshift,j+rshift};
       for p=1:2
       init(i,j,p) = mean(log(par1(p,5001:end))); 
      end
      for p=3:8
       init(i,j,p) = mean(par1(p,5001:end)); 
      end
    end
end
end

%Take initial values from independent runs
meanGinit=nanmean(reshape(init(:,:,1),L,1));
sdGinit = nanmean(reshape(init(:,:,2),L,1));
R0init = nanmean(reshape(init(:,:,3),L,1));
Kinit = nanmean(reshape(init(:,:,4),L,1));
ninit = nanmean(reshape(init(:,:,5),L,1));
deginit=nanmean(reshape(init(:,:,6),L,1));
kappainit=nanmean(reshape(init(:,:,7),L,1));
sigmaeinit=nanmean(reshape(init(:,:,8),L,1));

REinit=NaN(height,width,5);
scale = 0.1;
REinit(:,:,1) = fillmissing(fillmissing(init(:,:,1)-meanGinit,'linear',2),'linear',1);
REinit(:,:,1) = (REinit(:,:,1) - mean(mean(REinit(:,:,1))))*scale;
REinit(:,:,2) = fillmissing(fillmissing(init(:,:,2)-sdGinit,'linear',2),'linear',1);
REinit(:,:,2) = (REinit(:,:,2) - mean(mean(REinit(:,:,2))))*scale;
REinit(:,:,3) = fillmissing(fillmissing(init(:,:,3)-R0init,'linear',2),'linear',1);
REinit(:,:,3) = (REinit(:,:,3) - mean(mean(REinit(:,:,3))))*scale;
REinit(:,:,4) = fillmissing(fillmissing(init(:,:,4)-Kinit,'linear',2),'linear',1);
REinit(:,:,4) = (REinit(:,:,4) - mean(mean(REinit(:,:,4))))*scale;
REinit(:,:,5) = fillmissing(fillmissing(init(:,:,5)-ninit,'linear',2),'linear',1);
REinit(:,:,5) = (REinit(:,:,5) - mean(mean(REinit(:,:,5))))*scale;

thetainit = NaN(height,width,5);
thetainit(:,:,1) = meanGinit*ones(height,width);
thetainit(:,:,2) = sdGinit*ones(height,width);
thetainit(:,:,3) = R0init*ones(height,width);
thetainit(:,:,4) = Kinit*ones(height,width);
thetainit(:,:,5) = ninit*ones(height,width);

%Prior means and variances

    %rho    tau(meanG)   tau(sdG)       tau(R0)      tau(Kd)     tau(n)       ub meanG     ub sdG       R0         Kd          n            degr         kappa      sigmae      
 a = [1        0             0           0           0           0            24          20          0          0            0             -0.55           0          -5.3]; %prior ub/means for log(param)
    %rho    std(tau(1))    std(tau(2))  std(tau(3))  std(tau(4))  std(tau(5)) lb(meanG)    lb(sdG)    std(R0)   std(Kd)      std(n)         var(degr)  var(kappa)  var(sigmae) 
 b = [0        10            10          10         10           10             0           0          10         10          10            0.25^2         100       0.17^2];  %prior lb/stdev 

samp = 0.5;

%rhoinit = 0.01;
tauinit = log([5 5 5 5 5]);
psiinit(:,:,1) = deginit*ones(height,width);
psiinit(:,:,2) = kappainit*ones(height,width);
psiinit(:,:,3) = sigmaeinit*ones(height,width);

%Subset again, take odd rows and columns
data = data(1:2:end,1:2:end,:);

%remove artefacts from location 42
data = reshape(data,size(data,1)*size(data,2),288);
for i=1:288
if data(42,i)>0.75
    data(42,i) = (data(42,i-1) + data(42,i+1))./2;
end   
end    
data = reshape(data,45,25,288);

thetainit = thetainit(1:2:end,1:2:end,:);
REinit    = REinit(1:2:end,1:2:end,:);
psiinit   = psiinit(1:2:end,1:2:end,:);

%rng(10)
rng(13)
%grid for likelihood approximation
dt=0.5;
%number of mcmc iterations
B =40000;

clear deginit dshift height i init j kappainit Kinit kpriormean kpriorsd L meanGinit ninit notrend_reCRY1 p par1 PAR1 R0init rshift scale sdGinit sigmaeinit t T width

%OUT.k0 = 1;
OUT = load('OUT')  

tic
[tauest, thetameanest, REest, thetaest, psiest,tauprop, REprop, thetameanprop,psiprop, LIK0, tauaccept, REaccept, thetameanaccept] = carMCMC(data,dt,samp,B,a,b, thetainit, REinit, psiinit, tauinit, OUT.k0, OUT);
toc  
save('OUT CRY1 Exp4.mat','tauest', 'thetaest', 'REest', 'thetameanest', 'psiest', 'a', 'b', 'B', 'dt', 'samp','LIK0','tauprop','REprop','psiprop','thetameanprop', 'tauaccept', 'REaccept', 'thetameanaccept', '-v7.3');

