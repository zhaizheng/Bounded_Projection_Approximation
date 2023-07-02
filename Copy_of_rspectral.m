function [X, Y, Z, g, h] = ADMM(A, rho, k, m)

    Z = zeros(size(A));
    [~, X] = principal_k(A, k);
    Y = sG(X+Z/rho);
    
    iter = 1;
    error(iter) = norm(Y-X);
    [g(iter),h(iter)] = f(X, A, rho,m);
    
    while error(iter) >1.e-9 && iter<100
        X = principal_k((2.*A)+(rho.*Y)-Z, k);
        Y = sG(X+Z/rho);
        Z = Z + rho*(X-Y);
        iter = iter+1;
        error(iter) = norm(X-Y);
        
        [g(iter),h(iter)] = f(X,Y,Z,A,G);
    end

    function [g1,g2] = f(X,Y,Z,A,G)
        g1 = norm(X-Y);
        g2 = norm(A-X,'fro')^2;
    end


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


end

    
function [U,P] = principal_k(A, k)
    [U, D]= eig(A);
    [~, ind] = sort(diag(D),'descend');
    U = U(:,ind(1:k));
    P = U*U';
end


