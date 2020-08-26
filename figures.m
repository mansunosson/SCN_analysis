%Simulation study Figures

%% %Thetamean trace plots
fig1 = figure(1) 
subplot(2,3,1)
plot(exp(thetameanest(1,:)))
hold on
plot([1 B], [mean(mean(param.thetaest(:,:,1,1))) mean(mean(param.thetaest(:,:,1,1)))])
title('Delay mean')
subplot(2,3,2)
plot(exp(thetameanest(2,:)))
hold on
plot([1 B], [mean(mean(param.thetaest(:,:,1,2))) mean(mean(param.thetaest(:,:,1,2)))])
title('Delay SD')
subplot(2,3,3)
plot(exp(thetameanest(3,:)))
hold on
plot([1 B], [mean(mean(param.thetaest(:,:,1,3))) mean(mean(param.thetaest(:,:,1,3)))])
title('$R_0$', 'Interpreter', 'Latex')
subplot(2,3,4)
plot(exp(thetameanest(4,:)))
hold on
plot([1 B], [mean(mean(param.thetaest(:,:,1,4))) mean(mean(param.thetaest(:,:,1,4)))])
subplot(2,3,5)
plot(exp(thetameanest(5,:)))
hold on
plot([1 B], [mean(mean(param.thetaest(:,:,1,5))) mean(mean(param.thetaest(:,:,1,5)))])




%% 


subplot(2,3,1)
plot(exp(thetamean(1,:)))
hold on
plot([1 B], [9 9])
title('Delay mean')
subplot(2,3,2)
plot(exp(thetamean(2,:)))
hold on
plot([1 B], [3 3])
title('Delay SD')
subplot(2,3,3)
plot(exp(thetamean(3,:)).^2./exp(thetamean(4,:)))
hold on
plot([1 B], [60 60])
title('$R_0$', 'Interpreter', 'Latex')
subplot(2,3,4)
plot(exp(thetamean(3,:))./exp(thetamean(4,:)))
hold on
plot([1 B], [120 120])
subplot(2,3,5)
plot(exp(thetamean(5,:)))
hold on
plot([1 B], [5 5])



