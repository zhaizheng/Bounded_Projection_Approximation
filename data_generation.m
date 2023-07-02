function [S, A] = data_generation(signal, noise)
    k = 2;
    each_k = 20;
    n = k*each_k;
    S =  zeros(n,n);
    
    label = [];
    for j = 1:k
 %       S((j-1)*each_k+1:j*each_k,(j-1)*each_k+1:j*each_k) = 1;
        SIGNAL = rand(n/2,n/2) > signal;
        SIGNAL_Temp1 = triu(SIGNAL,1);
        SIGNAL_Temp2 = triu(SIGNAL);
        S((j-1)*each_k+1:j*each_k,(j-1)*each_k+1:j*each_k) = SIGNAL_Temp2 +SIGNAL_Temp1';
        label = [label, j*ones(1,each_k)];
    end
    E =  rand(n/2,n/2)> noise;
    E2 = [zeros(n/2),E;E',zeros(n/2)];
    A = [ones(n/2),zeros(n/2);
        zeros(n/2),ones(n/2)];
    S = S+E2;
end
