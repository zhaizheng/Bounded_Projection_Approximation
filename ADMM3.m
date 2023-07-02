function result = ADMM3(A, rho, k, a, b)

    Z = zeros(size(A));
    [~, initX] = principal_k(A, k);
    X = initX;
    %Y = sG(X+Z/rho, m);
    Y = SG((A+rho*X+Z)/(1+rho), a,b);
    
    iter = 1;
    error(iter) = norm(Y-X);
    [g(iter),h(iter)] = f(X,Y,Z,A);
    r(iter) = norm(X-initX);
    
    while error(iter) >1.e-9 && iter<500
        [U, X] = principal_k(A+(rho.*Y)-Z, k);
        Y = SG((A+rho*X+Z)/(1+rho), a,b);
        Z = Z + rho*(X-Y);
        iter = iter+1;
        error(iter) = norm(X-Y);
        
        [g(iter),h(iter)] = f(X,Y,Z,A);
        r(iter) = norm(X-initX);
    end
    result.U = U;
    result.X = X;
    result.Y = Y;
    result.Z = Z;
    result.g = g;
    result.h = h;
    result.r = r;
    
    function [g1,g2] = f(X,Y,Z,A)
        g1 = norm(X-Y);
        g2 = 0.5*norm(A-X,'fro')^2+0.5*norm(A-Y,'fro')^2;
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

    
    function X = SG(X,a,b)
        X(X>a) = a;
        X(X<b) = 0;
    end


end

    
function [U,P] = principal_k(A, k)
    [U, D]= eig(A);
    [~, ind] = sort(diag(D),'descend');
    U = U(:,ind(1:k));
    P = U*U';
end


