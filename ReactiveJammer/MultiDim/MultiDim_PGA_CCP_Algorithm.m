%% stochastic PGA_CCP convergence for multidimensional state

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

xhat_1 = rand(m,1);
xhat_0 = rand(m,1);

beta = 0.5;
alpha = 0.5;

conditions = [];

c = 1; % communication cost
d = 1; % jamming cost
var=1; % variance

batchsize=10000;

a = 0.1; % step size

theta = [alpha;beta];


n_MonteCarlo = 20;

n_iteration = 400;

conditions_PGACCP_1 = zeros(n_MonteCarlo,n_iteration);
conditions_PGACCP_2 = zeros(n_MonteCarlo,n_iteration);

for i_MonteCarlo = 1:n_MonteCarlo
    
    i_MonteCarlo
    xhat_1 = rand(m,1);
    xhat_0 = rand(m,1);
    
    beta  = rand(1);
    alpha = rand(1);
    theta = [alpha;beta];
    
    k=1;
    
    while k<=n_iteration
        
        
        %% projected gradient ascent PGA
        %       v = theta + (a/sqrt(k))*batch_grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
        v = theta + a*batch_grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
        
        theta_new=max(0,min(v,1)); % Projection
        
        theta = theta_new;
        
        alpha = theta(1);
        
        beta = theta(2);
        %% convex-concave procedure CCP
        
        A = [2*(1-alpha)*eye(m) zeros(m); zeros(m) 2*(beta+alpha)*eye(m)];
        
        g = batch_grad_CCP(alpha,beta,xhat_0,xhat_1);
        
        xhat_new = pinv(A)*g;
        
        xhat = xhat_new;
        
        xhat_0 = xhat(1:m);
        
        xhat_1 = xhat(m+1:end);
        
        %% Stopping criteria
        [Delta1,Delta2] = batch_FirstNashEquilibriumChecker(alpha,beta,xhat_0,xhat_1);
        
        conditions_PGACCP_1(i_MonteCarlo,k) =  Delta1;
        conditions_PGACCP_2(i_MonteCarlo,k) =  Delta2;
        
        k=k+1;
    end
end
%% Plot
x=1:1:n_iteration;
figure

semilogy(mean(conditions_PGACCP_1),'b--')
hold on
semilogy(mean(conditions_PGACCP_2),'b-')



