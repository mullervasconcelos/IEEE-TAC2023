%% Optimal jamming probabilities alpha* and beta* as a function of sigma^2

clear all
clc
addpath(genpath('./utils/'));
warning off

global c d var

%% initialization

c = 1; % communication cost
d = 1.2; % jamming cost

n=50;
step=0.1;
Beta=zeros(n,1);
Alpha=zeros(n,1);

epsilon=1e-4;

%% Main
for i=1:n
    i
    var=i*step;
    
    xhat_1 = rand;
    xhat_0 = rand;
    
    beta = 0.5;
    alpha = 0.5;
    
    t = 0.05;
    
    theta = [alpha;beta];
    
    flag = 0;
    

    
    while flag==0
        
       %% projected gradient ascent PGA
        %   v = theta + (t/sqrt(k))*grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
        v = theta + t*grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
        
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
        
        flag=((Delta1<=epsilon)&&(Delta2<=epsilon)); % Check whether this is a epsilon first nash equilibrium
       
    end
    
    Alpha(i)=alpha;
    Beta(i)=beta;  
end

%% Plot
figure

plot([step:step:step*n],Beta,'--','color',[0.85,0.33,0.10])

hold on

plot([step:step:step*n],Alpha,'-','color',[0.00,0.45,0.74])

grid on

set(gca,'xtick',[0.5:0.5:5],'ylim',[0,1])

xlabel('$\sigma^2$','interpreter','latex')

ylabel('Jamming Probability')

legend('$\beta^\star$','$\alpha^\star$','interpreter','latex')


