function re = ADMMnm(A, lambda, k, a,b, initA)

    n = size(A,1);

    rho = (3*sqrt(2))*lambda;
    Z = zeros(size(A));

    [initU, initX] = principal_k(initA, k);
    X = initX;
    
    r = rho/(2*lambda+rho);
    Temp = X + (Z/rho);
    Y = min(max(Temp, r.*Temp+(1-r).*a),r.*Temp+(1-r).*b);
    
    Z =  Z + rho*(X-Y);
    
    
    iter = 1;
    error(iter) = norm(Y-X);
    [g(iter),h(iter)] = f(X,Y,Z,A,a,b,lambda,rho);
    
    while error(iter) >1.e-9 && iter<100
        [U, X] = principal_k((2.*A)+(rho.*Y)-Z, k);
        r = rho/(2*lambda+rho);
        Temp = X + (Z/rho);
        Y = min(max(Temp, r.*Temp+(1-r).*a),r.*Temp+(1-r).*b);
        
        Z = Z + rho*(X-Y);
        iter = iter+1;
        error(iter) = norm(X-Y);
        
        [g(iter),h(iter)] = f(X,Y,Z,A,a,b,lambda,rho);
        r(iter) = norm(X-initX);
    end
    re.X = X;
    re.Y = Y;
    re.Z = Z;
    re.g = g;
    re.h = h;
    re.r = r;
    re.U = U;
    re.initU = initU;

    function [g1,g2] = f(X,Y,Z,A,a,b,lambda,rho)
        g1 = norm(X-Y);
        %Tem = min(X-a,0);
        Tem = min(Y-a,0);
        Tem2 = min(b-Y,0);
        XY = X-Y;
        %g2 = norm(A-X,'fro')^2+(lambda*norm(Tem,'fro')^2);
        g2 = norm(A-X,'fro')^2+(lambda*norm(Tem,'fro')^2+lambda*norm(Tem2,'fro')^2)+norm(X-Y,'fro')^2*rho/2+Z(:)'*XY(:);
        
    end


end

    
% function [U,P] = principal_k(A, k)
%     [U, D]= eig(A);
%     [~, ind] = sort(diag(D),'descend');
%     U = U(:,ind(1:k));
%     P = U*U';
% end

function [U,P] = principal_k(A, k)
    [U, ~]= eigs(A,k);
    P = U*U';
end



