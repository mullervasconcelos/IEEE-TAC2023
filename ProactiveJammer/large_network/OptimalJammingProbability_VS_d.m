%% Optimal jamming probabilities phi* as a function of d

clc
clear
c=0.01;
var=1;
kappa=0.2;

n=200;
step=0.01;
phi=zeros(n,1);
for i=1:n
    
    d=i*step;
    
    g1=@(x) exp(-x.^2/(2*var))./sqrt(2*pi*var);
    g2=@(x) x.^2.*exp(-x.^2/(2*var))./sqrt(2*pi*var);
    
    flag1=(2*integral(g1,sqrt(c),Inf)>=kappa);
    flag2=(2*integral(g2,sqrt(c),Inf)>=d);
    
    tol=1e-6; % Tolerance
    if flag2==0
        phi(i)=0;
    elseif flag1==0
        f=@(b) 2*integral(g2,sqrt(c/(1-b)),Inf)-d;
        [phi(i),f_opt,stepNum]=goldenOpt(f,0,1,tol);
    elseif flag1==1
        f1=@(l_l) 2*integral(g1,sqrt(l_l),Inf)-kappa;
        f2=@(l_p) 2*integral(g2,sqrt(l_p),Inf)-d;
        f1_opt=100;
        f2_opt=100;
        max_range=1e2;
        while max(abs(f1_opt),abs(f2_opt))>tol
            [l_lambda,f1_opt,~]=goldenOpt(f1,c,max_range,tol);
            [l_phi,f2_opt,~]=goldenOpt(f2,c,max_range,tol);
            max_range=max_range*10;
        end
        if l_lambda>l_phi
            phi(i)=0;
        elseif l_lambda<l_phi
            phi(i)=1-c/l_phi;
        else
            phi(i)=0;
        end
    end
end
%% Plot
figure
plot([step:step:n*step],phi)
xlabel('$d$','interpreter','latex')
ylabel('$\phi^\star$','interpreter','latex')
grid on
set(gca,'xtick',[0.2:0.2:2],'ytick',[0.1:0.1:1])

function [x_opt,f_opt,stepNum] = goldenOpt(f,a,b,Theta_error)

r=(sqrt(5)-1)/2;
a1=b-r*(b-a);
a2=a+r*(b-a);
stepNum=0;
while abs(b-a)>Theta_error
    stepNum=stepNum+1;
    f1=feval(f,a1);
    f2=feval(f,a2);
    if f1>0
        a=a1;
    end
    if f2<0
        b=a2;
    end
    a1=b-r*(b-a);
    a2=a+r*(b-a);
end
x_opt=(a+b)/2;
f_opt=feval(f,x_opt);
end
