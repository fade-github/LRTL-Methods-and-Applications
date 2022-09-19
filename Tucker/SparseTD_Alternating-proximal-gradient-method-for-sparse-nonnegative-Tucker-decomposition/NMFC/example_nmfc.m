% An example of nonnegative matrix factorization from partial observations
%% Generate problem data
rand('seed', 0); randn('seed', 0);
m = 500; n = 500; % matrix dimension
r = 20; % matrix rank
L = max(0,randn(m,r)); U = max(0,randn(r,n)); % nonnegative factors
Mtrue = L*U; % true matrix
sn = 60; % signal to noise ratio
% -- add noise --
N = max(0,randn(m,n));
M = Mtrue + 10^(-sn/20) * norm(Mtrue,'fro') * N/norm(N,'fro');
sr = 0.4; % percentage of samples
known = randsample(m*n,round(sr*m*n)); % randomly choose samples
data = M(known);

%% Solve problem
opts.maxit = 1000; % max # of iterations
opts.tol = 1e-4;   % stopping tolerance
% exact rank works but a rough rank overestimate is also fine
esr = round(r*1.25); % overestimated rank
t0 = tic;
[X, Y, Out] = nmfc(data,known,m,n,esr,opts);
time = toc(t0);

%% Reporting

fprintf('time = %4.2e, ', time);
fprintf('solution relative error = %4.2e\n\n', norm(X*Y - Mtrue,'fro')/norm(Mtrue,'fro'));

figure;
semilogy(1:Out.iter, Out.hist_obj, 'k-', 'LineWidth', 2);
xlabel('iteration','fontsize',12); ylabel('objective','fontsize',12);

figure;
semilogy(1:Out.iter, Out.hist_rel(2,:), 'k-', 'LineWidth', 2);
xlabel('iteration','fontsize',12); 
ylabel('Relative residual','fontsize',12);

