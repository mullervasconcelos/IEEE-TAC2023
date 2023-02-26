%% GDA algorithm convergence for multidimensional state

clear all
clc
addpath(genpath('./utils/'));
warning off

global c d var X

c = 1; % communication cost
d = 1; % jamming cost
var=1; % variance
seed=0; % random seed
rng(seed);

m=10;

batchsize=10000;

X=randn(m,batchsize)*var; % Samples from multivariate Gaussian distribution

n_MonteCarlo = 100;

n_iteration = 2000;

conditions_GDA_1    = zeros(n_MonteCarlo,n_iteration);
conditions_GDA_2    = zeros(n_MonteCarlo,n_iteration);

a1 = 0.1;

a2 = 0.01;



for i_MonteCarlo = 1:n_MonteCarlo
    
    i_MonteCarlo
    xhat_1 = rand(m,1);
    xhat_0 = rand(m,1);
    
    beta  = rand(1);
    alpha = rand(1);
    theta = [alpha;beta];
    
    k=1;
    
    while k<=n_iteration
        
        %% Gradient descent
        xhat = [xhat_0; xhat_1];
        
        g = batch_grad_GD(theta(1),theta(2),xhat_0,xhat_1);
        
        xhat_new = xhat - a2*g;
        
        
        %% Projected Gradient Ascent
        v = theta + a1*batch_grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
        theta_new=max(0,min(v,1));
        
        xhat = xhat_new;
        xhat_0 = xhat(1:m);
        xhat_1 = xhat(m+1:end);
        
        theta = theta_new;
        alpha = theta(1);
        beta = theta(2);
        
        %% Stopping criteria
        [Delta1,Delta2] = batch_FirstNashEquilibriumChecker(theta(1),theta(2),xhat_0,xhat_1);
        
        conditions_GDA_1(i_MonteCarlo,k) = Delta1;
        conditions_GDA_2(i_MonteCarlo,k) = Delta2;
        
        k = k+1;
    end
    
end

%% Plot
figure

[theta' xhat']

figure
semilogy(mean(conditions_GDA_1),'r--')
hold on
semilogy(mean(conditions_GDA_2),'r-')

grid on
legend({'$\|\nabla_{\hat{x}}\tilde{\mathcal{J}}(\hat{x}^{(k)},\varphi^{(k)})\|_2$-GDA','$\max_{\varphi \in[0,1]^2}\langle \nabla_{\varphi}\tilde{\mathcal{J}}(\hat{x}^{(k)},\varphi^{(k)}),\varphi-\varphi^{(k)}\rangle$-GDA'},'Interpreter','latex')
axis([0 n_iteration 1e-6 1e0])


