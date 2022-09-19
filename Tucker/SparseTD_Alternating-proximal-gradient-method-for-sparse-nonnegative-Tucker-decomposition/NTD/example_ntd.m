addpath('E:\Matlab Files\Toolbox\tensor_toolbox-v3.3')
savepath
% An example of nonnegative Tucker decomposition
%% Generate synthetic 3-order tensor
Nway = [50,50,50]; % dimension of tensor
coreNway = [5,5,5]; % dimension of core tensor

% randomly generate core tensor
G = tensor(max(0,randn(coreNway)));
A = cell(1,ndims(G));
% randomly generate factor matrices
for i = 1:ndims(G)
    A{i} = max(0,randn(Nway(i),coreNway(i)));
end
% generate tensor
Mtrue = full(ttensor(G,A)); N = ndims(Mtrue);

sn = 60; % signal to noise ratio
% -- add noise --
Noise = tensor(max(0,randn(Nway)));
M = Mtrue + 10^(-sn/20)*norm(Mtrue)/norm(Noise)*Noise;

%% Solve problem
opts.maxit = 1000; opts.tol = 1e-4;
t0 = tic;
[A,C,Out] = ntd(M,coreNway,opts);
time = toc(t0);

%% Reporting
relerr = norm(full(ttensor(C,A))-Mtrue)/norm(Mtrue);
fprintf('time = %4.2e, ',time);
fprintf('solution relative error = %4.2e\n\n',relerr);

figure;
semilogy(1:Out.iter, Out.hist_obj,'k-','linewidth',2);
xlabel('iteration','fontsize',12);
ylabel('objective value','fontsize',12)

figure;
semilogy(1:Out.iter, Out.hist_rel(2,:),'k-','linewidth',2);
xlabel('iteration','fontsize',12);
ylabel('relative residual','fontsize',12)

