%[S, real_A] = data_generation(0.4, 0.7,[30,50]);


addpath('./funs/')
%s = [30,20,30];
%s= [30,25,30];
s = [25 25 25];
noise = [0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8];
%signal = [0.42,0.51,0.54,0.57,0.6,0.63,0.66,0.69,0.72];
signal = [0.45,0.51,0.57,0.63,0.69];
repk = 100;


%class = 2; lab = [ones(1,20),2*ones(1,20)]';
RESULT = zeros(10,10);
rep = 30;
IMPROVE = zeros(10,10);
IMPROVE2 = zeros(10,10);
for i = 1%:length(signal)
    tic
    Result = zeros(10,10);
    Intermedia_result = zeros(rep,10);
    Intermedia_result2 = zeros(rep,10);
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

        fname = 'ell_F^2'; eta = sum(s.^2); c = 3; theta = 3; tau = 0.001;
        [Z,U] = SSL2(fname,S,eta,c,theta,tau); 
        [U4,~] = principal_k(Z, K);
        [Acc4, Nmi4] = rep_kmeans(U4, class, lab, repk);
        Intermedia{4} = Z;


        lambda = 0.12; n_iter = 100000;
        [A, ~, ~] = CLR_zz(S, lambda, K, n_iter);
        [U5,~] = principal_k(A+A', K);
        [Acc5, Nmi5] = rep_kmeans(U5, class, lab, repk);
        Intermedia{5} = A+A';
        init = zeros(size(S));
        for o = 1:5
            init = init+Intermedia{o};
        end
        rho = 10:30:130;
        rho = 10^4;
        for ss = 1%:5
            %MM = ADMMnm(S, 100000, K, 0.04*0.02, 0.98*0.04, init-Intermedia{2});  %U*U'+result1.X+Z
            %MM = ADMMnm(S, rho(ss), K, 0.05*0.02, 0.98*0.05, init-Intermedia{2}); 
            %%MM = ADMMnm(S, rho(ss), K, 0, 0.04, result1.X); 
            MM = ADMMnm(S, rho(ss), K, 0, 0.04, A+A'); 
            %MM = ADMMnm(S, 1000, K, 0, 0.04, init-Intermedia{2});  %U*U'+result1.X+Z
            WER = MM.X;
            [UMM,~] = principal_k(WER+WER', K);
            [Intermedia_result(j,ss), Intermedia_result(j,ss+5)] = rep_kmeans(UMM, class, lab, repk);

            MM2 = ADMMnm(S, rho(ss), K, 0, 0.04, Intermedia{1}); 
            %MM = ADMMnm(S, 1000, K, 0, 0.04, init-Intermedia{2});  %U*U'+result1.X+Z
            WER2 = MM2.X;
            [UMM2,~] = principal_k(WER2+WER2', K);
            [Intermedia_result2(j,ss), Intermedia_result2(j,ss+5)] = rep_kmeans(UMM2, class, lab, repk);
        end

        resultACC = [Acc1,Acc2,Acc3,Acc4, Acc5];
        resultNmi = [Nmi1,Nmi2,Nmi3,Nmi4, Nmi5];
        Result(j,:) = [resultACC, resultNmi];     

        %[resultACC,resultNmi]
    end
    %%
    if i == 1

        t = tiledlayout(1,7,'TileSpacing','Compact');
        nexttile
        imagesc(S)
        
        nexttile
        imagesc(result1.X)
        
        nexttile
        imagesc(result2.X)
        
        nexttile
        imagesc(P)

        nexttile
        imagesc(Z)

        nexttile
        imagesc(A+A')
        
        nexttile
        imagesc(MM.X)
    end

    RESULT(i,:) = mean(Result,1);
    IMPROVE(i,:) = mean(Intermedia_result,1);
    IMPROVE2(i,:) = mean(Intermedia_result2,1);
    elapsed = toc;
    fprintf('time cost is: %f second per iteration\n',elapsed)
end
%%
fprintf('\n')
FINAL = [];
FINAL = [FINAL,[RESULT(:,1:5), max(IMPROVE2(:,1:5),[],2), max(IMPROVE(:,1:5),[],2)]];
FINAL = [FINAL,[RESULT(:,6:10), max(IMPROVE2(:,6:10),[],2), max(IMPROVE(:,6:10),[],2)]];
% for i = 1:10
%     FINAL  = [FINAL,[RESULT(:,i),IMPROVE(:,i)]];
% end
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