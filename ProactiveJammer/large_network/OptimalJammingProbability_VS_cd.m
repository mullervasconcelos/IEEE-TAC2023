%% Optimal jamming probability phi* as a function of c and d

clc
clear

var=3;
kappa=0.5;

n=100;
step=0.03;
phi=zeros(n,n);
for i=1:n
    i
    d=i*step;
    for j=1:n
        c=j*step;
        
        g1=@(x) exp(-x.^2/(2*var))./sqrt(2*pi*var);
        g2=@(x) x.^2.*exp(-x.^2/(2*var))./sqrt(2*pi*var);
        
        flag1=(2*integral(g1,sqrt(c),Inf)>=kappa);
        flag2=(2*integral(g2,sqrt(c),Inf)>=d);
        
        tol=1e-6; % Tolerance
        if flag2==0
            phi(i,j)=0;
        elseif flag1==0
            f=@(b) 2*integral(g2,sqrt(c/(1-b)),Inf)-d;
            [phi(i,j),f_opt,stepNum]=goldenOpt(f,0,1,tol);
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
                phi(i,j)=0;
            elseif l_lambda<l_phi
                phi(i,j)=1-c/l_phi;
            else
                phi(i,j)=0;
            end
        end
    end
end
%% Plot
figure
% imagesc([step:step:n*step],[step:step:n*step],beta')
surf([step:step:n*step],[step:step:n*step],phi')
% colormap(gray)
ax=gca;
ax.YDir = 'normal'
xlabel('$d$','interpreter','latex')
ylabel('$c$','interpreter','latex')
zlabel('$\varphi^\star$','interpreter','latex')
% grid on
set(gca,'xtick',[0:20*step:n*step],'ytick',[0:20*step:n*step],'ztick',[0:0.2:1])
view([135 45]);

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
