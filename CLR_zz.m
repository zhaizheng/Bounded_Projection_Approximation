


% project_onto_positive_simplex([0.8,0.7,0.2,-0.3 0 ])
% s = [25 25 25];
% [S, A, label] = data_generation_unbal(0.51, 0.8, s);
% lambda = 1; K = 3; n_iter = 100;
% [S_re, err, f_value]= CLR_zzz(A, lambda, K, n_iter);
% subplot(1,4,1)
% imagesc(S);
% subplot(1,4,2)
% imagesc(S_re+S_re');
% subplot(1,4,3)
% semilogy(err)
% subplot(1,4,4)
% plot(f_value)


function [S,error,f_value] = CLR_zz(A, lambda, K, n_iter)
    error = [];
    f_value = [];
    S = A;
    m = size(A,1);
    for i = 1:m
        S(i,:) = project_onto_positive_simplex(A(i,:));
    end
    for i = 1:n_iter
        T = (S+S')/2;
        L = diag(sum(T,1))-T;
        [U, Singular, ~] = svd(L);
        F = U(:,m-K+1:m);
        DIST = sum(F.^2,2)*ones(1,m)+ones(m,1)*sum(F.^2,2)' - 2* F*F';
        for i = 1:m
            S_new(i,:) = project_onto_positive_simplex(A(i,:)-(lambda/2)*DIST(i,:));
        end
        if norm(S_new - S) < 1.e-5
            break;
        end
        error = [error, norm(S_new-S)];
        f_value= [f_value, norm(S-A,'fro')^2+2*lambda*sum(diag(Singular(m-K+1:m,m-K+1:m)))];
        S = S_new;
    end
end


function z = project_onto_positive_simplex(a)
    b = sort(a,'descend');
    for i = 2:length(b)
        test = sum(b(1:i-1)-b(i));
        if test > 1
            break;
        end
    end
    c = (sum(b(1:i-1))-1)/(i-1);
    z = max(a - c, 0);
end

