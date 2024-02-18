%Allen-Cahn equation in 1D+time with multistage DPG method%
mkdir(pwd,'Results_AllenCahn_DPG_FD')
mypath = fullfile(pwd,'Results_AllenCahn_DPG_FD');

r = 6; %Number of elements in each space direction: 2^r
m = 100; %Numer of time steps

%Data Burgers' equation
eps = 0.01;
x1 = -1;
x2 = 1;
T = 50;
u0 = @(x) 0.53*x+0.47*sin(-1.5*pi*x);

%Mesh, FD matrix and initial condition
nelx = 2^r;
xsol = linspace(x1,x2,nelx+1);
deltax = xsol(2)-xsol(1);
ex = ones(nelx-1, 1);
A = spdiags([ex -2*ex ex], -1:1, nelx-1, nelx-1);
A = -eps*A/deltax^2;
vx = xsol(2:end-1)';
U0 = u0(vx)-vx;
dimx = size(U0,1);

%Time grid and step size
steps = m;
tau = T/steps;
t = 0:tau:T;
t = t';

%Initialize energy vectors
En_HEuler = zeros(steps+1,1);
En_DPG2 = zeros(steps+1,1);
En_DPG3 = zeros(steps+1,1);

%Hybrid exponential Euler method
fprintf('Running hybrid Exponential Euler for r=%d and %d time steps\n\n',r,steps)
U = zeros(dimx,steps+1);
Usol = zeros(dimx+2,steps+1);
U(:,1) = U0;
Usol(:,1)=[-1;U0+vx;1];
En_HEuler(1) = compute_energy(Usol(:,1),eps,deltax);
U0int = zeros(dimx,steps);
for i = 1:steps
    Jn = spdiags(1-3*(U(:,i)+vx).^2,0,dimx,dimx)-A;
    Fun = ((U(:,i)+vx)-(U(:,i)+vx).^3)-A*U(:,i);
    Gn = Fun-Jn*U(:,i);
    U0int(:,i) = U(:,i)+tau*phiB(tau*Jn,[Fun zeros(dimx,1)]);
    U(:,i+1) = U(:,i)+tau*Jn*U0int(:,i)+tau*Gn; 
    Usol(:,i+1)=[-1;U(:,i+1)+vx;1];
    En_HEuler(i+1) = compute_energy(Usol(:,i+1),eps,deltax);
end

%Two-stage DPG method
fprintf('Running two-stage DPG method for r=%d and %d time steps\n\n',r,steps)
U = zeros(dimx,steps+1);
Usol = zeros(dimx+2,steps+1);
U(:,1) = U0;
Usol(:,1)=[-1;U0+vx;1];
En_DPG2(1) = compute_energy(Usol(:,1),eps,deltax);
Un2 = zeros(dimx,steps);
for i = 1:steps
    Jn = spdiags(1-3*(U(:,i)+vx).^2,0,dimx,dimx)-A;
    Fun = ((U(:,i)+vx)-(U(:,i)+vx).^3)-A*U(:,i);
    Gn = Fun-Jn*U(:,i);
    Un2(:,i) = U(:,i)+tau*phiB(tau*Jn,[Fun zeros(dimx,1)]);
    Fun2 = ((Un2(:,i)+vx)-(Un2(:,i)+vx).^3)-A*Un2(:,i);
    Gn2 = Fun2-Jn*Un2(:,i);
    U(:,i+1) = U(:,i)+tau*phiB(tau*Jn,[8*Gn2-8*Gn zeros(dimx,1) Fun]);
    Usol(:,i+1)=[-1;U(:,i+1)+vx;1];
    En_DPG2(i+1) = compute_energy(Usol(:,i+1),eps,deltax);
end

% Three-stage DPG method method
fprintf('Running three-stage DPG method for r=%d and %d time steps\n\n',r,steps)
U = zeros(dimx,steps+1);
Usol = zeros(dimx+2,steps+1);
U(:,1) = U0;
Usol(:,1)=[-1;U0+vx;1];
En_DPG3(1) = compute_energy(Usol(:,1),eps,deltax);
Un2 = zeros(dimx,steps);
Un3 = zeros(dimx,steps);
for i = 1:steps
    Jn = spdiags(1-3*(U(:,i)+vx).^2,0,dimx,dimx)-A;
    Fun = ((U(:,i)+vx)-(U(:,i)+vx).^3)-A*U(:,i);
    Gn = Fun-Jn*U(:,i);
    Un2(:,i) = U(:,i)+tau*phiB(tau*Jn,[Fun zeros(dimx,1)]);
    Un3(:,i) = U(:,i)+tau*Jn*Un2(:,i)+tau*Gn;
    %Evaluations
    Fun2 = ((Un2(:,i)+vx)-(Un2(:,i)+vx).^3)-A*Un2(:,i);
    Gn2 = Fun2-Jn*Un2(:,i);
    Fun3 = ((Un3(:,i)+vx)-(Un3(:,i)+vx).^3)-A*Un3(:,i);
    Gn3 = Fun3-Jn*Un3(:,i);
    %Correction factor
    Jn2 = spdiags(1-3*(Un2(:,i)+vx).^2,0,dimx,dimx)-A;
    dGn2 = Jn2-Jn;
    Cn2 = -(1/4)*dGn2*(Un3(:,i)-2*Un2(:,i)+U(:,i));
    %Update
    U(:,i+1) = U(:,i)+tau*phiB(tau*Jn,[-48*Gn2+12*Gn3+36*Gn-48*Cn2 16*Gn2-2*Gn3-14*Gn+16*Cn2 zeros(dimx,1) Fun]);
    Usol(:,i+1) = [-1;U(:,i+1)+vx;1];
    En_DPG3(i+1) = compute_energy(Usol(:,i+1),eps,deltax);
end
%plot_solution

%Save energy
Tab_energy = table(t,En_HEuler,En_DPG2,En_DPG3,'VariableNames',{'t','Energy_HEuler','Energy_DPG2','Energy_DPG3'});
writetable(Tab_energy,fullfile(mypath,['EnergyAllenCahn_DPG_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');

%Function to compute the energy at each time step
function En = compute_energy(U,eps,deltax)
    dim = size(U,1);
    En = 0;
    for i=2:dim
        En = En + ((eps/2)*((U(i)-U(i-1))/deltax)^2+(1/4)*(U(i)^2-1)^2)*deltax;
    end
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
