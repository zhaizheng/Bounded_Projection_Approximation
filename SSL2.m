function [Z,U] = SSL2(fname,A,eta,c,theta,tau)
%fprintf('SSL process by zyzhang');
n = size(A,1);
A = A*(sqrt(c)/norm(A,'fro'));
A_abs = abs(A);

Z = A; err = 1; iter = 0;
while err>tau 
    iter = iter+1;
    Z_old  = Z;
    
    % update U
    [U,Lambda] = eig(Z); lambda = real(diag(Lambda));
    [slambda,J] = sort(lambda); 
    U(:,J(1:n-c)) = [];
    if slambda(n-c+1)-slambda(n-c)<1.e-8
        warning('equal eigenvalues of Z')
    end
    
    % update Z
    W = U*U'-A;
    switch fname
        case 'ell_F^2'
            H = (theta/(1+theta))*W;
            Z = A+H;
            D = abs(Z);
            
        case 'ell_1'
            H = sign(W).*max(abs(W)-1/(2*theta),0);
            Z = A+H;
            D = A_abs-abs(H)+theta*((A+W).^2-(abs(H)-abs(W)).^2);
    end
    Z = trunc(Z,D,eta);    
    err = norm(Z-Z_old,'fro');

%     switch fname
%         case 'ell_F^2'
%             fa = norm(A-Z,'fro')^2;
%         case 'ell_1'
%             fa = norm(A-Z,1);
%     end
%     figure(2),spy(Z),pause

end
% 
    function Z = trunc(Z,D,eta)
        nz = size(Z,2);
        eta = (eta-nz)/2;
        D = triu(D,1);
        [sd,I] = sort(D(:),'descend');
%         if sd(eta)-sd(eta+1)<1.e-9
%             warning('equal entries in D')
%         end
        
        a = zeros(nz*nz,1);
        a(I(1:eta)) = 1;
      %  z = vec2mat(a,nz)';
        z = reshape(a,[nz,nz]);
        z = z+z'+eye(nz);
        Z = Z.*z;
    end

end
% function X = trunc(A,B,eta)
%     X = zeros(size(A));   
%     [s,~] = sort(B(:),'descend');
%     X(B>=s(eta)) = A(B>=s(eta));
% end