%% Monte Carlo simulations for PGA-CCP and GDA in 1 dimension

clear all
clc
addpath(genpath('./utils/'));
warning off

global c d var

c = 1; % communication cost
d = 1; % jamming cost
var=1; % variance

n_MonteCarlo = 50;

n_iteration = 2000;


conditions_PGACCP_1 = zeros(n_MonteCarlo,n_iteration);
conditions_PGACCP_2 = zeros(n_MonteCarlo,n_iteration);

conditions_GDA_1    = zeros(n_MonteCarlo,n_iteration);
conditions_GDA_2    = zeros(n_MonteCarlo,n_iteration);


for i_MonteCarlo = 1:n_MonteCarlo
    
    i_MonteCarlo
    
    x_init = rand(2,1);
    phi_init = rand(2,1);
    
    
    %% PGA-CCP Algorithm
    
    xhat_1 = x_init(1);
    xhat_0 = x_init(2);
    
    beta  =  phi_init(1);
    alpha =  phi_init(2);
    
    a = 0.1; % step size
    
    theta = [alpha;beta];
    
    flag = 0;
    
    k=1;
    
    while k<=n_iteration
        
        
        %% projected gradient ascent PGA
        %   v = theta + (t/sqrt(k))*proj_sub_grad(theta(1),theta(2),xhat_0,xhat_1);
        v = theta + a*grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
        
        theta_new=max(0,min(v,1)); % Projection
        
        theta = theta_new;
        
        alpha = theta(1);
        
        beta = theta(2);
        %% convex-concave procedure CCP
        
        A = [2*(1-alpha) 0; 0 2*(beta+alpha)];
        
        g = grad_CCP(alpha,beta,xhat_0,xhat_1);
        
        xhat_new = pinv(A)*g;
        
        xhat = xhat_new;
        
        xhat_0 = xhat(1);
        
        xhat_1 = xhat(2);
        
        %% Stopping criteria
        [Delta1,Delta2] = FirstNashEquilibriumChecker(alpha,beta,xhat_0,xhat_1);
     
        
        conditions_PGACCP_1(i_MonteCarlo,k) =  Delta1;
        conditions_PGACCP_2(i_MonteCarlo,k) =  Delta2;
        
        k=k+1;
    end
    
    
    %% GDA Algorithm
    
    xhat_1 = x_init(1);
    xhat_0 = x_init(2);
    
    beta  =  phi_init(1);
    alpha =  phi_init(2);
    
    a1 = 0.1;
    
    a2 = 0.01;
    
    theta = [alpha;beta];
    
    k=1;
    
    while k<=n_iteration
        
        %% Gradient descent
        xhat = [xhat_0; xhat_1];
        
        
        g = grad_GD(theta(1),theta(2),xhat_0,xhat_1);
        
        xhat_new = xhat - a2*g;
        
        
        
        %% Projected Gradient Ascent
        v = theta + a1*grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
        theta_new=max(0,min(v,1));
        
        
        xhat = xhat_new;
        xhat_0 = xhat(1);
        xhat_1 = xhat(2);
        
        theta = theta_new;
        alpha = theta(1);
        beta = theta(2);
        
        
        %% Stopping criteria
        [Delta1,Delta2] = FirstNashEquilibriumChecker(theta(1),theta(2),xhat_0,xhat_1);
        
        conditions_GDA_1(i_MonteCarlo,k) = Delta1;
        conditions_GDA_2(i_MonteCarlo,k) = Delta2;
        
        k = k+1;
    end
    
    
    
end
%% Plot
n_iteration = 2000;
x=1:1:n_iteration;
figure

conditions_GDA =max(conditions_GDA_1,conditions_GDA_2);

conditions_PGACCP =max(conditions_PGACCP_1,conditions_PGACCP_2);

patch([x,fliplr(x)],[mean(conditions_GDA)-std(conditions_GDA), ...
    fliplr(mean(conditions_GDA)+std(conditions_GDA))],[1.00,0.4,0.2],'FaceAlpha',0.5,'EdgeColor','none');
hold on
plot(mean(conditions_GDA),'--','color',[1.00,0.4,0.2])

patch([x,fliplr(x)],[mean(conditions_PGACCP)-std(conditions_PGACCP), ...
    fliplr(mean(conditions_PGACCP)+std(conditions_PGACCP))],[0.3 0.7 0.9],'FaceAlpha',0.5,'edgecolor','none');
plot(mean(conditions_PGACCP),'-','color',[0.3 0.7 0.9])

axis([0 500 0 2])
xlabel('Number of Iterations ($k$)','interpreter','latex')
ylabel('$\varepsilon$-FNE','interpreter','latex')

grid on
legend({'GDA std dev','GDA mean','PGA-CCP std dev','PGA-CCP mean'},'Interpreter','latex')