%Error convergence
figure
loglog(taum,Errors1,taum,Errors2,taum,Errors3)
legend('Hybrid Euler','DPG2','DPG3','Location','SouthEast')
set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',16)

%Compute convergence rates
x1=log(Errors1(2:end)./Errors1(1:end-1));
x2=log(Errors2(2:end)./Errors2(1:end-1));
x3=log(Errors3(2:end)./Errors3(1:end-1));
y=log(taum(2:end)./taum(1:end-1));
ConvRates_RK1=x1./y
ConvRates_RK2=x2./y
ConvRates_RK3=x3./y