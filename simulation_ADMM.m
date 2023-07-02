%
clear

[S, real_A] = data_generation(0.45, 0.7);
k = 2;
B = [0.01,0.001,0.0001,0.00001,0.000001];
Rho = [100,1000,10000,100000,1000000];
Error = zeros(length(B),length(Rho));
Error2 = Error;
result = cell(length(B),length(Rho),2);

for i = 1:length(B)
    for j = 1:length(Rho)
        result{i,j,1} = ADMMn(S, Rho(j), k, B(i));
        result{i,j,2} = ADMMnm(S, Rho(j), k, B(i),0.05);
        [~, P_real] = principal_k(real_A, k);
        Error(i,j,1) = norm(result{i,j,1}.X-P_real,'fro');
        Error2(i,j) = norm(result{i,j,2}.X-P_real,'fro');
    end
end

[~, P] = principal_k(S, k);
fprintf('principal error :%.3f',norm(P_real-P,'fro'));

[r,c] = find(Error2 == min(Error2(:)));
figure(1)

subplot(1,4,1)
imagesc(S);
subplot(1,4,2)
imagesc(result{r,c,2}.X);
subplot(1,4,3)
[~,P] = principal_k(S, k);
imagesc(P)
subplot(1,4,4)
hold on
semilogy(result{r,c,2}.h,'-','linewidth',2,'markersize',4);
%legend(leg);
box on

% figure(2)
% for i = 1
%     semilogy(re1{i}.g,'-','linewidth',2);
%     legend(leg);
%     hold on
% % end
% 
% [~, P_real] = principal_k(real_A, k);
% fprintf('Dis norm real and est: %.3f,real and princ: %.3f\n',norm(P_real-re1{1}.X,'fro'), norm(P_real-P,'fro'))
%%
figure
subplot(1,2,1)
plot(result{r,c,2}.initU(:,2));
subplot(1,2,2)
plot(result{r,c,2}.U(:,2));









     


function [S, A] = data_generation(signal, noise)
    k = 2;
    each_k = 20;
    n = k*each_k;
    S =  zeros(n,n);
    
    label = [];
    for j = 1:k
 %       S((j-1)*each_k+1:j*each_k,(j-1)*each_k+1:j*each_k) = 1;
        SIGNAL = rand(n/2,n/2) > signal;
        SIGNAL_Temp1 = triu(SIGNAL,1);
        SIGNAL_Temp2 = triu(SIGNAL);
        S((j-1)*each_k+1:j*each_k,(j-1)*each_k+1:j*each_k) = SIGNAL_Temp2 +SIGNAL_Temp1';
        label = [label, j*ones(1,each_k)];
    end
    E =  rand(n/2,n/2)> noise;
    E2 = [zeros(n/2),E;E',zeros(n/2)];
    A = [ones(n/2),zeros(n/2);
        zeros(n/2),ones(n/2)];
    S = S+E2;
end





%%   



function [error,U, P, g, h] = rspectral(A, lambda, k, m)

    [U, P] = principal_k(A, k);
    G = sG(P, m);
    [U, nP] = principal_k(lambda*G, k);
    iter = 1;
    error(iter) = norm(nP-P);
    [g(iter),h(iter)] = f(P, A, lambda,m);
    
    while error(iter) >1.e-9 && iter<100
        P = nP;
        G = sG(P, m);
        [U, nP] = principal_k((A+lambda*G), k);
        iter = iter+1;
        error(iter) = norm(nP-P);
        [g(iter),h(iter)] = f(nP, A, lambda, G);
    end

    function [g1,g2] = f(X, A,lambda,G)
        g2 = norm(X-G,'fro')^2; 
        g1 = norm(X-A,'fro')^2+lambda*norm(X-G,'fro')^2; 
    end

%     function NY = sG(X,m)
%         NY = max(X,0);
%     end

    function NY = sG(X,m)
        E = true(size(X));
        T = abs(X(:));
        %T = X(triu(E));
        [~,ind] = sort(T(:),'descend');
        T(ind(m+1:end)) = 0;
        Y = zeros(size(X));
        Y(E) = T(:);
        NY = Y;
    end

%     function NY = sG(X,m)
%         E = true(size(X));
%         T = abs(X(triu(E)));
%         %T = X(triu(E));
%         [~,ind] = sort(T(:),'descend');
%         T(ind(m+1:end)) = 0;
%         Y = zeros(size(X));
%         Y(triu(E)) = T(:);
%         NY = Y+Y'-diag(diag(Y));
%     end

end

    
function [U,P] = principal_k(A, k)
    [U, D]= eig(A);
    [~, ind] = sort(diag(D),'descend');
    U = U(:,ind(1:k));
    P = U*U';
end
