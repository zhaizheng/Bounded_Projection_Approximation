function [error, U, G, P, g, h] = rspectral(A, lambda, k, m)

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
        %[U, nP] = principal_k(G, k);
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
        T = X(:);
        %T = abs(X(:));
        %T = X(triu(E));
        [~,ind] = sort(T(:),'descend');
        T(ind(floor(m+1):end)) = 0;
        Y = zeros(size(X));
        Y(:) = T(:);
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


