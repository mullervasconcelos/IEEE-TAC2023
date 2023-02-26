function [g] = batch_grad_GD(alpha,beta,xhat_0,xhat_1)

global c d X




%% Accelerated version

index=(beta*sum((X-xhat_1).^2)+c-d*beta <= alpha*sum((X-xhat_1).^2) + (1-alpha)*sum((X-xhat_0).^2)-d*alpha);

g=[mean((-2*(1-alpha)*(X-xhat_0)).*(1-index),2);mean((-2*beta*(X-xhat_1)).*index,2)+mean((-2*alpha*(X-xhat_1)).*(1-index),2)];


% Standard version
% s1=0;
% s2=0;
% for i_batch=1:batchsize
% 
% 	x = randn(length(xhat_0),1)*var;
% 	
% 	s1 = s1 + (-2*(1-alpha)*(x-xhat_0))*(beta*norm(x-xhat_1)^2+c-d*beta > alpha*norm(x-xhat_1)^2 + (1-alpha)*norm(x-xhat_0)^2-d*alpha);
% 
% 	s2 = s2 + (-2*beta*(x-xhat_1)).*(beta*norm(x-xhat_1)^2+c-d*beta <= alpha*norm(x-xhat_1)^2 + (1-alpha)*norm(x-xhat_0)^2-d*alpha) +...
%         (-2*alpha*(x-xhat_1)).*(beta*norm(x-xhat_1)^2+c-d*beta > alpha*norm(x-xhat_1)^2 + (1-alpha)*norm(x-xhat_0)^2-d*alpha);
% 
% end
% 
% g = [s1/batchsize; s2/batchsize];