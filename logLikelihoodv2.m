%LogLikelihood given by EKBF for Per Transcription
function logLik = logLikelihoodv2(dt,data,R0,Kp,np,degr,meanG,sdG,maxd,kappa,sigmae,samp)
                              
if (R0<0)|| (Kp<0) || (np<0) || (degr<0) || (meanG<0) || (meanG>23) || (sdG<0) || (sdG>20) ||  (kappa<0)|| (sigmae<0) 
        logLik=-Inf;     
 else
        [~,~,meanout,varout,~,~] = EKBFv2(dt,data,R0,Kp,np,degr,meanG,sdG,maxd,kappa,sigmae,samp);
        vv=reshape(varout,length(data),1);    
            if isreal(vv)==1 && sum(sum(isnan(vv)))==0 && sum(vv<0)==0
            logLik=sum(log(normpdf(data(49:end)',meanout(49:end)',sqrt(vv(49:end))))); %maxd/samp + 1  = 49 
            else
               logLik=NaN;
            end
    end
end
