%Burgers' equation in 1D+time with multistage DPG method%
mkdir(pwd,'Results_Burgers_DPG_UpW')
mypath = fullfile(pwd,'Results_Burgers_DPG_UpW');

r = 10; %Number of elements in each space direction: 2^r
m = 8; %Number of time steps: 2^m

%Data equation
x1 = 0;
x2 = 1;
T = 3;
u0 = @(x) (1/(4*pi))*sin(2*pi*x);

%Mesh and initial condition
nel = 2^r;
xsol = linspace(x1,x2,nel+1);
deltax = xsol(2)-xsol(1);
U0 = u0(xsol(1:end-1)');
dimx = size(U0,1);

%Time grid and step size
steps = 2^m;
tau = T/steps;
t = 0:tau:T;
t = t';

%Initialize energy vectors
En_HEuler = zeros(steps+1,1);
En_DPG2 = zeros(steps+1,1);
En_DPG3 = zeros(steps+1,1);

%Hybrid exponential Euler method
fprintf('Running hybrid exponential Euler for r=%d and %d time steps\n\n',r,steps)
U = zeros(dimx,steps+1);
Usol = zeros(dimx+1,steps+1);
U(:,1) = U0;
Usol(:,1) = [U0;U0(1)];
En_HEuler(1) = sum((0.5*Usol(:,1).^2)*deltax);
U0int = zeros(dimx,steps);
for i = 1:steps
    i
    [Jn,Fun] = compute_lin(U(:,i),dimx,deltax);
    Gn = Fun-Jn*U(:,i);
    U0int(:,i) = U(:,i)+tau*phiB(tau*Jn,[Fun zeros(dimx,1)]);
    U(:,i+1) = U(:,i)+tau*Jn*U0int(:,i)+tau*Gn; 
    Usol(:,i+1) = [U(:,i+1);U(1,i+1)];
    En_HEuler(i+1)= sum((0.5*Usol(:,i+1).^2)*deltax);
end

%Two-stage DPG method
fprintf('Running two-stage DPG method  for r=%d and %d time steps\n\n',r,steps)
U = zeros(dimx,steps+1);
Usol = zeros(dimx+1,steps+1);
U(:,1) = U0;
Usol(:,1) = [U0;U0(1)];
En_DPG2(1) = sum((0.5*Usol(:,1).^2)*deltax);
Un2 = zeros(dimx,steps);
for i = 1:steps
    i
    [Jn,Fun] = compute_lin(U(:,i),dimx,deltax);
    Gn = Fun-Jn*U(:,i);
    Un2(:,i) = U(:,i)+tau*phiB(tau*Jn,[Fun zeros(dimx,1)]);
    [Jn2,Fun2] = compute_lin(U(:,i),dimx,deltax);
    Gn2 = Fun2-Jn*Un2(:,i);
    U(:,i+1) = U(:,i)+tau*phiB(tau*Jn,[8*Gn2-8*Gn zeros(dimx,1) Fun]);
    Usol(:,i+1) = [U(:,i+1);U(1,i+1)];
    En_DPG2(i+1)= sum((0.5*Usol(:,i+1).^2)*deltax);
end

%Three-stage DPG method
fprintf('Running three-stage DPG method for r=%d and %d time steps\n\n',r,steps)
U = zeros(dimx,steps+1);
Usol = zeros(dimx+1,steps+1);
U(:,1) = U0;
Usol(:,1) = [U0;U0(1)];
En_DPG3(1) = sum((0.5*Usol(:,1).^2)*deltax);
Un2 = zeros(dimx,steps);
Un3 = zeros(dimx,steps);
for i = 1:steps
    i
    [Jn,Fun] = compute_lin(U(:,i),dimx,deltax);
    Gn = Fun-Jn*U(:,i);
    Un2(:,i) = U(:,i)+tau*phiB(tau*Jn,[Fun zeros(dimx,1)]);
    Un3(:,i) = U(:,i)+tau*Jn*Un2(:,i)+tau*Gn;
    %Evaluations
    [Jn2,Fun2] = compute_lin(Un2(:,i),dimx,deltax);
    Gn2 = Fun2-Jn*Un2(:,i);
    [Jn3,Fun3] = compute_lin(Un3(:,i),dimx,deltax);
    Gn3 = Fun3-Jn*Un3(:,i);
    %Correction factor
    dGn2 = Jn2-Jn;
    Cn2 = -(1/4)*dGn2*(Un3(:,i)-2*Un2(:,i)+U(:,i));
    %Update
    U(:,i+1) = U(:,i)+tau*phiB(tau*Jn,[-48*Gn2+12*Gn3+36*Gn-48*Cn2 16*Gn2-2*Gn3-14*Gn+16*Cn2 zeros(dimx,1) Fun]);
    Usol(:,i+1) = [U(:,i+1);U(1,i+1)];
    En_DPG3(i+1) = sum((0.5*Usol(:,i+1).^2)*deltax);
end

%Save energy
Tab_energy = table(t,En_HEuler,En_DPG2,En_DPG3,'VariableNames',{'t','Energy_HEuler','Energy_DPG2','Energy_DPG3'});
writetable(Tab_energy,fullfile(mypath,['EnergyBurgers_DPG_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');
%plot_solution

%Function to compute the Jacobian and the nonlinear remainder 
function [Jn,Fun] = compute_lin(U,dimx,deltax)
    Jn = zeros(dimx,dimx);
    Fun = zeros(dimx,1);
    if U(1)>0
        Fun(1) = -U(1).*(U(1)-U(end))/deltax;
        Jn(1,1) = 2*U(1)-U(end);
        Jn(1,end) = -U(1);
    else
        Fun(1) = -U(1).*(U(2)-U(1))/deltax;
        Jn(1,1) = -2*U(1)+U(2);
        Jn(1,end) = U(1);
    end
    for j = 2:dimx-1
        if U(j)>0
            Fun(j) = -U(j).*(U(j)-U(j-1))/deltax;
            Jn(j,j) = 2*U(j)-U(j-1);
            Jn(j,j-1) = -U(j);
        else
            Fun(j) = -U(j).*(U(j+1)-U(j))/deltax;
            Jn(j,j) = -2*U(j)+U(j+1);
            Jn(j,j+1) = U(j);
        end
    end
     if U(end)>0
        Fun(end) = -U(end).*(U(end)-U(end-1))/deltax;
        Jn(end,end) = 2*U(end)-U(end-1);
        Jn(end,end-1) = -U(end);
    else
        Fun(end) = -U(end).*(U(1)-U(end))/deltax;
        Jn(end,end) = -2*U(end)+U(1);
        Jn(end,1) = U(end);
     end
    Jn = sparse(-Jn*(1/deltax));

end

%Function to compute the linear combination of actions with expmv
function out = phiB(A,B)
    q=size(B,2);
    if q == 1
        At = [A, B; zeros(q, length(A)) 0];
        if issparse(A)
            At = sparse(At);
        end  
    else
        if issparse(A)
            At = sparse([A, B; sparse(zeros(q,length(A))), sparse(diag(ones(q-1,1),1))]);
        else
            At = [A, B; (zeros(q,length(A))), (diag(ones(q-1,1),1))];
        end
    end
    e = zeros(length(At),1); e(length(A)+q) = 1;
    out = expmv(1,At,e);
    out = out(1:length(A));
end