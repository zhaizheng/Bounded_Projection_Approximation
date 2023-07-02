function [y, S,vals] = CLR_wei(vals,alpha, A00, c, lambda,ln,S0)
% This function is a modification of the code provided by the following work. 
% Ref:
% Feiping Nie, Xiaoqian Wang, Michael I. Jordan, Heng Huang.
% The Constrained Laplacian Rank Algorithm for Graph-Based Clustering.
% The 30th Conference on Artificial Intelligence (\textbf{AAAI}), Phoenix, USA, 2016.
viewnum = size(alpha,2);

if nargin<7
   S0 = zeros(size(A00{1},1));
   for v = 1:viewnum
       S0 = S0+alpha(1,v)*A00{v};
   end   
end

S0 = S0-diag(diag(S0));
num = size(S0,1);
S10 = (S0+S0')/2;
D10 = diag(sum(S10));
L0 = D10 - S10;

NITER = 50;
zr = 10e-11;

[F0, ~, evs]=eig1(L0, num, 0);
F = F0(:,2:(c+1));

    % caculate objvals
    obj = objectfun(S10,A00,F,L0,lambda,ln);
    vals = [vals,obj];
for iter = 1:NITER
    dist = L2_distance_1(F',F');
    S = zeros(num);
    for i=1:num
        a0 = zeros(1,num);
        for v = 1:viewnum
            temp = A00{v};
            a0 = a0+alpha(1,v)*temp(i,:);
        end
        
        idxa0 = find(a0>0);
        ai = a0(idxa0);
        di = dist(i,idxa0);

        ad = (ai-0.5*lambda*di)/sum(alpha);
        S(i,idxa0) = EProjSimplex_new(ad);
    end;

    A = S;
    A = (A+A')/2;
    D = diag(sum(A));
    L = D-A;
    F_old = F; % store F temporaly
    [F, ~, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;
    % caculate objvals
    obj = objectfun(A,A00,F,L,lambda,ln);
    vals = [vals,obj];
   
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > zr
        lambda = 2*lambda;
    elseif fn2 < zr
        lambda = lambda/2;
        F = F_old;
    else
        break;
    end; 
    
%     if (iter>1 && abs(vals(end)-vals(end-1)) < 10^-10)
%         break;
%     end    
end;
 
%[labv, tem, y] = unique(round(0.1*round(1000*F)),'rows');
[clusternum, y]=graphconncomp(sparse(A)); y = y';
if clusternum ~= c
    sprintf('Can not find the correct cluster number: %d', c)
end;


