%%Plot solution%%
figure
for i=1:steps+1
    plot(xsol,Usol(:,i))
    title('t=',t(i))
    pause
end

% figure
% [Xgrid,Tgrid] = meshgrid(xsol,t);
% surf(Xgrid,Tgrid,Usol');
% xlabel({'$x$'},'interpreter','latex')
% ylabel({'$t$'},'interpreter','latex')
% colormap('jet')
% set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',16)