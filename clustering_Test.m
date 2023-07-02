clear
A = load('mfeat-pix'); 
%A = extract();
n = 1000; dim = 5; class = n/200; repk=100; 
%C = build_graph2(A);
C = build_graph(A(1001:2000,:));
lab = label_build(n);
[F,e] = spectral_projection(C, dim);
theta = 1:3:10;
beta = 0.70:0.1:1.1;
spar = 2*10^5; 
%m = 1.8*10^5; 
[Acc1, Nmi1] = rep_kmeans(F, class, lab', repk);
fprintf('spectral acc=%f, nmi=%f\n', Acc1, Nmi1);
for j = 1:length(beta)
    m = beta(j)*spar;
    rho = 10; k = 5;
    result = ADMM(C, rho, k, m);
    [Acc2(j), Nmi2(j)] = rep_kmeans(result.U, class, lab', repk);
    fprintf('beta=%f.3,Acc=%f.3\n', beta(j), Acc2(j));
end
%fprintf('theta=%f.3,beta=%f.3,Acc=%f.3\n',theta(i), beta(j), Acc2(i,j));


%%
%[Acc3, Nmi3] = rep_kmeans(F2, class, lab', repk);
%fprintf('truncated ACC=%f, NMI=%f\n', Acc3, Nmi3);

%%
f1 = non_diag(F*F', lab);
f2 = non_diag(result.U*result.U', lab);
fprintf('Eigen Refine non-diag=%f, original=%f\n', f2, f1);
%%
figure(2)
subplot(1,2,1); imagesc(F*F')
subplot(1,2,2); imagesc(result.U*result.U')
%subplot(1,3,3); imagesc(abs(result.U*result.U'-F*F'))

%%
f_o = non_diag(F*F', lab);
for i = 1:6
    for j = 1:10
        f(i,j) = non_diag(P{i,j}, lab)/f_o;
    end
end


function A = extract()
    data = load('mnist.mat');
    A = [];
    for i = 1:5
        Temp = data.trainX(data.trainY == i-1,:);
        A = [A; Temp(1:200,:)];
    end
end


function f = non_diag(A, lab)
    a = unique(lab);
    L = zeros(length(lab));
    for i = 1:length(a)
        L(lab==a(i),lab==a(i)) = 1;
    end
    B = A.*(ones(length(lab))-L);
    f = norm(B,'fro')^2;
end


function [Acc, Nmi] = rep_kmeans(F, class, lab, repk)
    for i = 1:repk
        idx = kmeans(F, class);
        [acc(i),~,~] = AccMeasure(lab',idx');
        nmii(i) = nmi(lab', idx');
    end
    Acc = mean(acc)/100;
    Nmi = mean(nmii);
end


function lab = label_build(n)
    label = zeros(1000,1);
    for i = 1:5
        label((i-1)*200+1:i*200) = i;
    end
    lab = label(1:n);
end


function C = build_graph(A)
    data = double(A);
    data = data-mean(data,1);
    sd = sum(data.^2,2);
    e = ones(size(A, 1),1);
    dis = sd*e'+e*sd'-2*data*data';
    r = mean(dis(:));
    C = exp(-dis/r);
end


function C = build_graph2(A)
    data2 = zeros(size(A));
    data2(A>0) = 1;
    %data = double(A);
    %datan = diag(1./sqrt(sum(data2.^2,2)))*data2;
    C = data2*data2';
end

function [F,e] = spectral_projection(C, dim)
    [V,D] = eig(C);
    [~,ind] = sort(diag(D),'descend');
    F = V(:,ind(1:dim));
    e1 = diag(D);
    e = e1(ind);
end

