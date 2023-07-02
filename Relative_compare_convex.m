%[S, real_A] = data_generation(0.4, 0.7,[30,50]);


addpath('./funs/')
%s= [30,20,30];
s = [25 25 25]
noise = [0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8];
signal = [0.42,0.45,0.48,0.51,0.54,0.57,0.6,0.63];
repk = 100;


% t = tiledlayout(1,5,'TileSpacing','Compact');
% nexttile
% xlabel('fadf')
% imagesc(S)
% 
% nexttile
% xlabel('fadf')
% imagesc(result1.X)
% 
% nexttile
% imagesc(result2.X)
% 
% nexttile
% imagesc(P)
% 
% nexttile
% imagesc(result3.X)

%class = 2; lab = [ones(1,20),2*ones(1,20)]';
RESULT = zeros(10,10);
rep = 10;
IMPROVE = zeros(10,10);
for i = 1:8
    Result = zeros(10,10);
    Intermedia_result = zeros(rep,10);
    Intermedia = cell(1,5);
    for j = 1:rep
        [S, A, label] = data_generation_unbal(signal(i), noise(i), s);

        n = sum(s);
        %S = S;
        K = 3; class = K; lab = label;

        result1 = ADMM_SD1(S, n, K, 5);
        Intermedia{1} = result1.X;
        [U1,~] = principal_k(result1.X, K);
        [Acc1, Nmi1] = rep_kmeans(U1, class, lab, repk);

        result2 = ADMM_SD2(S, n, K, 5);
        Intermedia{2} = result2.X;
        [U2,~] = principal_k(result2.X, K);
        [Acc2, Nmi2] = rep_kmeans(real(U2), class, lab, repk);
        
        [U,P] = principal_k(S, K);
        Intermedia{3} = S;
        [Acc3, Nmi3] = rep_kmeans(U, class, lab, repk);

        fname = 'ell_F^2'; eta = sum(s.^2); c = 3; theta = 1; tau = 0.001;
        [Z,U] = SSL2(fname,S,eta,c,theta,tau); 
        [U4,~] = principal_k(Z, K);
        [Acc4, Nmi4] = rep_kmeans(U4, class, lab, repk);
        Intermedia{4} = Z;


        lambda = 1; n_iter = 1000;
        [A, ~, ~] = CLR_zz(S, lambda, K, n_iter);
        [U5,~] = principal_k(A+A', K);
        [Acc5, Nmi5] = rep_kmeans(U5, class, lab, repk);
        Intermedia{5} = A+A';
        
        for ss = 1:5
            MM = ADMMnm(Intermedia{ss}, 100000, K, 0, 0.04, Intermedia{ss});  %U*U'+result1.X+Z
            WER = MM.X;
            [UMM,~] = principal_k(WER+WER', K);
            [Intermedia_result(j,ss), Intermedia_result(j,ss+5)] = rep_kmeans(UMM, class, lab, repk);
        end

        resultACC = [Acc1,Acc2,Acc3,Acc4, Acc5];
        resultNmi = [Nmi1,Nmi2,Nmi3,Nmi4, Nmi5];
        Result(j,:) = [resultACC, resultNmi];
        %[resultACC,resultNmi]
    end
    mean(Result,1)
    mean(Intermedia_result,1)
    RESULT(i,:) = mean(Result,1);
    IMPROVE(i,:) = mean(Intermedia_result,1);
    RESULT
    IMPROVE
end
FINAL = [];
for i = 1:10
    FINAL  = [FINAL,[RESULT(:,i),IMPROVE(:,i)]];
end

format_convert(FINAL)





function [Acc, Nmi] = rep_kmeans(F, class, lab, repk)
    F = diag(1./(sqrt(sum(F.^2,2))+eps))*F;
    for i = 1:repk
        idx = kmeans(real(F), class);
        [acc(i),~,~] = AccMeasure(lab',idx');
        nmii(i) = nmi(lab', idx');
    end
    Acc = mean(acc)/100;
    Nmi = mean(nmii);
end




    
function [U,P] = principal_k(A, k)
    [U, D]= eig(A);
    [~, ind] = sort(diag(D),'descend');
    U = U(:,ind(1:k));
    P = U*U';
end