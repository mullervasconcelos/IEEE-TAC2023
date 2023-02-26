function [g] = batch_grad_PGA(alpha,beta,xhat_0,xhat_1)

global c d X


%% Accelerated version

index=(beta*sum((X-xhat_1).^2)+c-d*beta <= alpha*sum((X-xhat_1).^2) + (1-alpha)*sum((X-xhat_0).^2)-d*alpha);

g=[mean((sum((X-xhat_1).^2)-sum((X-xhat_0).^2)-d).*(1-index),2);mean((sum((X-xhat_1).^2)-d).*index,2)];


% Standard version
% s1=0;
% s2=0;
% for i_batch=1:batchsize
% 
% 	x=randn(length(xhat_0),1)*var;
% 
% 	s1 = s1 + (norm(x-xhat_1)^2-norm(x-xhat_0)^2-d)*(beta*norm(x-xhat_1)^2+c-d*beta > alpha*norm(x-xhat_1)^2 + (1-alpha)*norm(x-xhat_0)^2-d*alpha);
% 
% 	s2 = s2 + (norm(x-xhat_1)^2-d)*(beta*norm(x-xhat_1)^2+c-d*beta <= alpha*norm(x-xhat_1)^2 + (1-alpha)*norm(x-xhat_0)^2-d*alpha);
% 
% end
% 
% g = [s1/batchsize; s2/batchsize];