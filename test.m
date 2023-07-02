clear
A = [0.5 0.5 0 0;
    0.5 0.5 0 0;
    0 0 0.5 0.5;
    0 0 0.5 0.5];
B = round(rand(4,4)-0.5,2);
format_convert(A)
format_convert(B+B')
C = ADMMnm(A+B+B',10^9,2,0,0.5);
format_convert(C.X)
[U,~,~] = svd(A+B+B');
Temp = U(:,1:2);
format_convert(Temp*Temp')


