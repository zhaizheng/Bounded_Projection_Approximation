clear
name = {'COIL10','COIL20','DIGIT5','DIGIT10'}
for i = 1:4
addpath('./Data/')
Data = load(name{i});
X = Data.fea;
%X = X(1:1000,:);
lab = Data.gnd;
%lab = lab(1:1000,:);


NX = diag(1./sqrt(sum(X.^2,2)))*X;
n = length(lab);


s = sum(X.^2,2);
D = s*ones(1,n)+ones(n,1)*s'-2*X*X';
ave = sum(D(:))/(n^2);
M = exp(-D/ave);

K = length(unique(lab));
result1 = ADMM_SD1(M, n, K, 5);
result2 = ADMM_SD2(M, n, K, 5);
%%
class = K;  repk = 50;

[U1,~] = principal_k(result1.X, K);
[Acc1, Nmi1] = rep_kmeans(U1, class, lab, repk);
%%
[U2,~] = principal_k(result2.X, K);
[Acc2, Nmi2] = rep_kmeans(U2, class, lab, repk);

[U3,P3] = principal_k(M, K);
[Acc3, Nmi3] = rep_kmeans(U3, class, lab, repk);


%%
fname = 'ell_F^2'; eta = floor(n^2/K);  theta = 1; tau = 0.001;
[Z,U] = SSL2(fname,M,eta, K, theta, tau); 
[U4,~] = principal_k(Z, K);
[Acc4, Nmi4] = rep_kmeans(U4, class, lab, repk);



lambda = 0.01; n_iter = 100;
[A, ~, ~] = CLR_zz(M, lambda, K, n_iter);
[U5,~] = principal_k(A+A', K);
[Acc5, Nmi5] = rep_kmeans(U5, class, lab, repk);


%%
result6 = ADMMnm(A+A', 1000000, K, 0, 1.01*(K/n), U5*U5');
[U6,~] = principal_k(result6.X, K);
[Acc6, Nmi6] = rep_kmeans(U6, class, lab, repk);
fprintf([name{i},'\n']);

format_convert([Acc1,Acc2,Acc3,Acc4,Acc5,Acc6, Nmi1,Nmi2,Nmi3,Nmi4,Nmi5,Nmi6])

t = tiledlayout(1,7,'TileSpacing','Compact');

% nexttile
% imagesc(M)
% 
% nexttile
% imagesc(result1.X)
% 
% nexttile
% imagesc(result2.X)
% 
% nexttile
% imagesc(P3)
% 
% nexttile
% imagesc(Z)
% 
% nexttile
% imagesc(A+A')
% 
% nexttile
% imagesc(result6.X)

end



function [U,P] = principal_k(A, k)
    [U, D]= eig(A);
    [~, ind] = sort(diag(D),'descend');
    U = U(:,ind(1:k));
    P = U*U';
end



function [Acc, Nmi] = rep_kmeans(F, class, lab, repk)
    F = diag(1./(sqrt(sum(F.^2,2))+eps))*F;
    lab = round(lab);
    for i = 1:repk
        idx = round(kmeans(F, class));
        [acc(i),~] = calAC(lab',idx');
        %[acc(i),~,~] = AccMeasure(lab',idx');
        nmii(i) = calMI(idx',lab');
        %nmii(i) = nmi(lab', idx');
    end
    Acc = mean(acc);
    Nmi = mean(nmii);
end
