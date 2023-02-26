
bound=6;
x0=[-bound:0.1:bound];

N=length(x0);

J=zeros(1,N);

phi=0.5;

c=1;

d=1;
x1=0;

for n=1:N
    f = @(x) min((1-phi)*(x-x0(n)).^2,c).*exp(-x.^2/2)/sqrt(2*pi);

    J(n) = integral(f,-Inf,Inf)+phi*(1-d);
    
end

plot(x0,J,'linewidth',1.5)
axis([-6,6,0.3,1.1])
xlabel('$\hat{x}_0$','interpreter','latex')
ylabel('$\tilde{J}$','interpreter','latex')