function val = objectfun(S,A,F,L,lambda,ln)
    % caculate order graph
if ln == 0
    order = length(A);
    s0 = 0;
    for o = 1:order
        s0 = s0 + norm(S-A{o},'fro');
    end
    val = s0 + 2*lambda*trace(F'*L*F);
elseif ln == 1
    order = length(A);
    s0 = 0;
    for o = 1:order
        s0 = s0 + log(norm(S-A{o},'fro')^2+1);
    end
    val = s0 + 2*lambda*trace(F'*L*F);
end

%     order = size(alpha,2);
%     s0 = 0;
%     for o = 1:order
%         s0 = s0 + log(norm(S-A{o},'fro')+1);
%     end
%     
%     D = diag(sum(S));
%     L = D-S;
%     val = s0 + 2*lambda*trace(F'*L*F);
% end