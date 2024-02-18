%Hochburck-Osterman Equation in 2D+time with exponential Rosebrock methods%
mkdir(pwd,'Results_HochOster_Rosenbrock_FD')
mypath = fullfile(pwd,'Results_HochOster_Rosenbrock_FD');

rmax = 4; %Maximum elements in each space direction: 2^r
mmax = 8; %Maximum time step size: 2^m

%Data of the problem 
eps = 1;
u0 = @(x,y) x.*(1-x).*y.*(1-y);
x1 = 0;
x2 = 1;
y1 = 0;
y2 = 1;
T = 1;
uexact = @(x,y,t) exp(t).*x.*(1-x).*y.*(1-y);
source = @(x,y,t) uexact(x,y,t)+2*exp(t).*(x.*(1-x)+y.*(1-y))-1./(1+exp(2*t).*(x.*(1-x).*y.*(1-y)).^2);
dtsource = @(x,y,t) uexact(x,y,t)+2*exp(t).*(x.*(1-x)+y.*(1-y))+(2*exp(2*t).*(x.*(1-x).*y.*(1-y)).^2)./(1+exp(2*t).*(x.*(1-x).*y.*(1-y)).^2).^2;

for r = rmax:rmax 
    
    %Initialize vectors for the results
    Errors1 = zeros(mmax,1);
    Errors2 = zeros(mmax,1);
    Errors3 = zeros(mmax,1);
    times_Euler = zeros(mmax,1);
    times_RK2 = zeros(mmax,1);
    times_RK3 = zeros(mmax,1);
    taum = zeros(mmax,1);
    nsteps = zeros(mmax,1);
    
    %Meshes and Parameters
    nelx = 2^r;
    nely = 2^r;
    xsol = linspace(x1,x2,nelx+1);
    ysol = linspace(y1,y2,nely+1);
    dx = xsol(2)-xsol(1);
    dy = ysol(2)-ysol(1);

    %Finite Differences matrix for the Laplacian
    ex = ones(nelx-1, 1);
    ey = ones(nely-1, 1);
    Ax = spdiags([ex -2*ex ex], -1:1, nelx-1, nelx-1);
    Ay = spdiags([ey -2*ey ey], -1:1, nely-1, nely-1);
    A = -eps*(kron(Ay, speye(nely-1))/dx^2 + kron(speye(nelx-1), Ax)/dy^2);

    %Initial condition vector
    [xmat,ymat] = meshgrid(xsol(2:end-1),ysol(2:end-1));
    U0 = u0(xmat,ymat);
    U0 = reshape(U0',[],1);
    dimx = size(U0,1);
    
    %Exact solution at time final time
    UexactT = uexact(xmat,ymat,T);
    UexactT = reshape(UexactT',[],1);
    
    for m = 1:mmax %tau=T/2^m time step size  
    
        %Time grid and step size
        steps = 2^m;
        nsteps(m) = steps;
        tau = T/steps;
        taum(m) = tau;
        t = 0:tau:T;
    
        %Assembly the linear source term
        fprintf('Assembly of the RHS \n\n')
        flin = zeros(dimx,steps+1);
        vn = zeros(dimx,steps);
        flin12 = zeros(dimx,steps);
        flin13 = zeros(dimx,steps);
        vn2 = zeros(dimx,steps);
        for i = 1:steps
            flin(:,i) = eval_int(source,xmat,ymat,t(i));
            vn(:,i) = eval_int(dtsource,xmat,ymat,t(i));
            flin12(:,i) = eval_int(source,xmat,ymat,t(i)+tau/2);
            flin13(:,i) = eval_int(source,xmat,ymat,t(i)+tau);
            vn2(:,i) = eval_int(dtsource,xmat,ymat,t(i)+tau/2);
        end
        flin(:,end) = eval_int(source,xmat,ymat,t(end));
            
        %Exponential Euler method
        fprintf('Running Exponential Euler for r=%d and %d time steps\n\n',r,steps)
        U = zeros(dimx+1,steps+1);
        U(1,1) = t(1);
        U(2:end,1) = U0;
        tic
        for i = 1:steps
            Jn = spdiags(-2*U(2:end,i)./(1+U(2:end,i).^2).^2,0,dimx,dimx)-A;
            JnB = [0 zeros(1,dimx); vn(:,i) Jn];
            Fun = [1; 1./(1+U(2:end,i).^2)+flin(:,i)-A*U(2:end,i)];
            U(:,i+1) = U(:,i)+tau*phiB(tau*JnB,Fun);
        end
        times_Euler(m) = toc;
        Errors1(m) = norm(UexactT-U(2:end,end),'inf');

        %RK2 method
        fprintf('Running RK2 for r=%d and %d time steps\n\n',r,steps)
        U = zeros(dimx+1,steps+1);
        U(1,1) = t(1);
        U(2:end,1) = U0;
        Un2 = zeros(dimx+1,steps);
        tic
        for i = 1:steps
            Jn = spdiags(-2*U(2:end,i)./(1+U(2:end,i).^2).^2,0,dimx,dimx)-A;
            JnB = [0 zeros(1,dimx); vn(:,i) Jn];
            Fun = [1; 1./(1+U(2:end,i).^2)+flin(:,i)-A*U(2:end,i)];
            Gn = Fun-JnB*U(:,i);
            Un2(:,i) = U(:,i)+tau*phiB(tau*JnB,Fun);
            Fun2 = [1; 1./(1+Un2(2:end,i).^2)+flin13(:,i)-A*Un2(2:end,i)];
            Gn2 = Fun2-JnB*Un2(:,i);
            U(:,i+1) = U(:,i)+tau*phiB(tau*JnB,[2*Gn2-2*Gn zeros(dimx+1,1) Fun]);
        end
        times_RK2(m) = toc;
        Errors2(m) = norm(UexactT-U(2:end,end),'inf');

     
        %RK3 method
        fprintf('Running RK3 for r=%d and %d time steps\n\n',r,steps)
        U = zeros(dimx+1,steps+1);
        U(1,1) = t(1);
        U(2:end,1) = U0;
        Un2 = zeros(dimx+1,steps);
        Un3 = zeros(dimx+1,steps);
        tic
        for i = 1:steps
            Jn = spdiags(-2*U(2:end,i)./(1+U(2:end,i).^2).^2,0,dimx,dimx)-A;
            JnB = [0 zeros(1,dimx); vn(:,i) Jn];
            Fun = [1; 1./(1+U(2:end,i).^2)+flin(:,i)-A*U(2:end,i)];
            Gn = Fun-JnB*U(:,i);
            %second stage
            Un2(:,i) = U(:,i)+(tau/2)*phiB((tau/2)*JnB,Fun);
            Fun2 = [1; 1./(1+Un2(2:end,i).^2)+flin12(:,i)-A*Un2(2:end,i)];
            Gn2 = Fun2-JnB*Un2(:,i);
            %third stage
            Un3(:,i) = U(:,i)+tau*phiB(tau*JnB,Gn2+JnB*U(:,i));
            Fun3 = [1; 1./(1+Un3(2:end,i).^2)+flin13(:,i)-A*Un3(2:end,i)];
            Gn3 = Fun3-JnB*Un3(:,i);
            %Update
            U(:,i+1) = U(:,i)+tau*phiB(tau*JnB,[-48*Gn2+12*Gn3+36*Gn 16*Gn2-2*Gn3-14*Gn zeros(dimx+1,1) Fun]);
        end
        times_RK3(m) = toc;
        Errors3(m) = norm(UexactT-U(2:end,end),'inf');
    end
    
        %Save errors and computational times
        Tab_times = table(nsteps,times_Euler,times_RK2,times_RK3,'VariableNames',{'Nsteps','Times_Euler','Times_RK2','Times_RK3'});
        writetable(Tab_times,fullfile(mypath,['TimesHochOster_Rosenbrock_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');
        
        Tab_errors = table(taum,Errors1,Errors2,Errors3,'VariableNames',{'tau','Err_Euler','Err_RK2','Err_RK3'});
        writetable(Tab_errors,fullfile(mypath,['ErrHochOster_Rosenbrock_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');

end

%Function to evaluate the source vector in the internal nodes 
function eval = eval_int(source,xmat,ymat,t)
    eval = source(xmat,ymat,t);
    eval = reshape(eval',[],1);
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
