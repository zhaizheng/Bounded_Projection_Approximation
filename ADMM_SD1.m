function re = ADMM_SD1(A, n, K, rho)
    Z = zeros(n);
    Y = zeros(n);
    U = zeros(n);
    V = zeros(n);
    %X = project_onto_affine2(0.5*(Z-U+Y-V+A./rho),n, K);
    X = project_onto_affine(0.5*(rho*Z-U+rho*Y-V+A)/rho,n, K);
    Z = max(X+U/rho,0);
    Y = projection_SD(X+V/rho);
    U = U+rho*(X-Z);
    V = V+rho*(X-Y);
    iter = 1;
    h(iter) = -X(:)'*A(:)+(X(:)-Z(:))'*U(:)+(X(:)-Y(:))'*V(:)+rho*(norm(X-Z,'fro')^2+norm(X-Y,'fro')^2)/2;
    
    while norm(X-Z)+norm(X-Y)>10^(-4)
        %X = project_onto_affine(0.5*(Z-U+Y-V+A./rho),n, K);
        X = project_onto_affine(0.5*(rho*Z-U+rho*Y-V+A)/rho,n, K);
        %Z = max(X+U,0);
        Z = max(X+U/rho,0);
        Y = projection_SD(X+V/rho);
        U = U+rho*(X-Z);
        V = V+rho*(X-Y);
        h(iter) = -X(:)'*A(:)+(X(:)-Z(:))'*U(:)+(X(:)-Y(:))'*V(:)+rho*(norm(X-Z,'fro')^2+norm(X-Y,'fro')^2)/2;
        iter = iter+1;
    end
    re.X = X;
    re.Z = Z;
    re.Y = Y;
    re.U = U;
    re.V = V;
    re.h = h;
end


function SX = projection_SD(X)
    [U,D] = eig(X);
    SX = U*max(D,0)*U';
end

function RE = project_onto_affine(Y, n, K)
    b = [ones(n,1)*2*(n/K-1);ones(n,1)];
    LY = [2*(sum(Y,2)-diag(Y)); diag(Y)];
    LL = zeros(2*n);
    LL(1:n, 1:n) = 1/(2*n-4)*(eye(n)-ones(n)./(2*n-2));
    LL(n+1:2*n,n+1:2*n) = eye(n);
    ve = LL*(LY-b);
    RE = Y-(ve(1:n)+ve(1:n)'-2*diag(ve(1:n))+diag(ve(n+1:2*n)));
end


