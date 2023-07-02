

function [S, A, label] = data_generation_unbal(signal, noise, s)
    K = length(s);
    n = sum(s);
    S =  zeros(n,n);
    label = [];
    for j = 1:K
 %       S((j-1)*each_k+1:j*each_k,(j-1)*each_k+1:j*each_k) = 1;
        each_k = s(j);
        SIGNAL = rand(s(j),s(j)) > signal;
        SIGNAL_Temp1 = triu(SIGNAL,1);
        SIGNAL_Temp2 = triu(SIGNAL);
       % S((j-1)*each_k+1:j*each_k,(j-1)*each_k+1:j*each_k) = SIGNAL_Temp2 +SIGNAL_Temp1';
       if j>1
           S(sum(s(1:j-1))+1:sum(s(1:j)), sum(s(1:j-1))+1:sum(s(1:j))) = SIGNAL_Temp2 +SIGNAL_Temp1';
       else 
           S(1:sum(s(1:j)), 1:sum(s(1:j))) = SIGNAL_Temp2 +SIGNAL_Temp1';
       end
        label = [label, j*ones(1,each_k)];
    end

    for j = 1:K
        for k = j+1:K
            NOISE = rand(s(j),s(k)) > noise;
            if j>1
                S(sum(s(1:j-1))+1:sum(s(1:j)), sum(s(1:k-1))+1:sum(s(1:k))) = NOISE;
                S(sum(s(1:k-1))+1:sum(s(1:k)),sum(s(1:j-1))+1:sum(s(1:j))) = NOISE';
            else
                S(1:sum(s(1)),sum(s(1:k-1))+1:sum(s(1:k))) = NOISE;
                S(sum(s(1:k-1))+1:sum(s(1:k)),1:sum(s(1))) = NOISE';
            end
        end
    end
    
    A = S;
    
%     E =  rand(n/2,n/2)> noise;
%     E2 = [zeros(n/2),E;E',zeros(n/2)];
%     A = [ones(n/2),zeros(n/2);
%         zeros(n/2),ones(n/2)];
%     S = S+E2;
end