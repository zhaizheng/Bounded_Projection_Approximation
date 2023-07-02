% function [y,S,w,alpha,clusternum] = Horder_CLR(A,c,alpha_k) %(X,c)
% function [y,S,w,alpha,vals] = Horder_CLR(X,W,c,order,knn,alpha_k,ln) %(X,c)
function [y,S,w,alpha,vals] = Horder_CLR(A,c,order,knn,alpha_k,ln) %(X,c)
% X \in R^{d*n} c:classnum
% [dim,num]=size(X);
% ln = 0: l_p norm;  ln = 1: ln() function

if nargin < 6
    ln = 0;
end
%% 构图
% order = length(A); %高阶阶数
% knn = ln;

%% CLR 初始化
alpha = 1/order*ones(1,order);
% alpha = zeros(1,order);
% alpha(1) = 1;
NITER = 5; % max loop
% lambda = 40; % the penalty parameter on graph 
lambda = 0.1;
vals = [];
for iter = 1:NITER
    
    % fix alpha, update S
    if iter ==1
       [y, S,vals] = CLR_wei(vals,alpha,A,c,lambda,ln);
    else
       [y, S,vals] = CLR_wei(vals,alpha,A,c,lambda,ln,S0);
    end
    
    % fix S, update alpha
    if ln == 0
        for o = 1:order
            alpha(o) = 0.5/norm(S-A{o},'fro');
        end
    elseif ln == 1
        for o = 1:order
            alpha(o) = 1/(norm(S-A{o},'fro')^2+1);
        end
    end
    w = alpha;
    alpha = selec_max(alpha,alpha_k);
    
    S0 = S;
    % calculate obj
%     obj = objectfun(S,A,F,lambda);
%     vals = [vals,obj];
    if (iter>1 && abs(vals(end)-vals(end-1)) < 10^-6)
        break;
    end      
end
%     lambda = 100;
%     [y, S] = CLR(alpha,A,classnum,lambda,S0);


end

function vec = selec_max(alpha,k)
    vec = zeros(1,size(alpha,2));
    for i =1:k
        col = find(alpha==max(alpha));
        vec(col) = alpha(col);
        alpha(col) = 0;
    end
end

