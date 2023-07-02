%
exp_graph(0.6,1,2)

exp_graph(0.8,3,4)



function exp_graph(sparse,fi_k,fi_s)
    k = 2;
    each_k = 20;
    n = k*each_k;
    S =  zeros(n,n);
    label = [];
    for j = 1:k
        S((j-1)*each_k+1:j*each_k,(j-1)*each_k+1:j*each_k) = 1;
        label = [label, j*ones(1,each_k)];
    end
    E =  rand(n/2,n/2)>sparse;
    E2 = [zeros(n/2),E;E',zeros(n/2)];

    mu = 0.5;

    S = mu*S+(1-mu)*E2;%(E+E')/2;

    m = 800;
    Lambda = 1:5:50;
    for i = 1:10
        [error{i},U, P{i}, g{i}, h{i}] = rspectral(S, Lambda(i), k, m);
    end
    figure(fi_k)
    subplot(1,4,1)
    imagesc(S)
    for i = 1:3
        subplot(1,4,i+1)
        imagesc(P{(i-1)*4+1})
    end

    figure(fi_s)
    subplot(1,2,1)
    leg= cell(1,10);
    for i = 1:10
        hold on
        semilogy(error{i},'-','linewidth',3,'markersize',4)
        xlabel('iteration')
        ylabel('$\|X_{n+1}-X_n\|_F^2$','interpreter','latex')
        leg{i} = ['\lambda=',num2str(Lambda(i))];
    end
    legend(leg)
    box on
    subplot(1,2,2)
    for i = 1:10
        hold on
        semilogy(h{i}(2:end),'-','linewidth',3,'markersize',4)
        xlabel('iteration')
        ylabel('$\|X_n-Y_n\|_F^2$','interpreter','latex')
    end
    legend(leg)
    box on
%     subplot(2,3,5)
%     scatter(U(:,1),U(:,2),[],label,'filled');
% 
%     subplot(2,3,6)
%     [U1,P1] = principal_k(S, k);
%     scatter(U1(:,1),U1(:,2),[],label,'filled');
end





%%   



function [error,U, P, g, h] = rspectral(A, lambda, k, m)

    [U, P] = principal_k(A, k);
    G = sG(P, m);
    [U, nP] = principal_k(A+lambda*G, k);
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
