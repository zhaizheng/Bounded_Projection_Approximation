clear

A = load('mfeat-pix'); n = 1000; dim = 5; class = n/200; repk=100; 
lab = label_build(n);



figure(1)
distribution(0, A, lab)

function distribution(k,A,lab)
    %subplot(1,4,k+1)
    n = 1000;
    C = build_graph(A, 1000, 'tail');
    dim = 2;
    [F,~] = spectral_projection(C, dim);
    %scatter(F(:,1),F(:,2),[],lab(1:n));
    %plot(F(:,1),F(:,2),'*');
    %subplot(1,4,k+2)
    s = [0.0002,0.002]*4;
    for i = 1:5:1000
         hold on
        small = reshape(A(i,:),[15,16]);
        image(small*40,'XData',[F(i,1),F(i,1)+s(1)],'YData',[F(i,2),F(i,2)+s(2)])
    end
    axis off
end




function [F,e] = spectral_projection(C, dim)
    [V,D] = eig(C);
    [~,ind] = sort(diag(D),'descend');
    F = V(:,ind(1:dim));
    e1 = diag(D);
    e = e1(ind);
end


function lab = label_build(n)
    label = zeros(2000,1);
    for i = 1:10
        label((i-1)*200+1:i*200) = i;
    end
    lab = label(1:n);
end


function C = build_graph(A, n, part)
    if strcmp(part, 'front')
        data = A(1:n,:);
    elseif strcmp(part, 'tail')
        data = A(n+1:end,:);
    end
    data = data-mean(data,1);
    sd = sum(data.^2,2);
    e = ones(n,1);
    dis = sd*e'+e*sd'-2*data*data';
    r = mean(dis(:));
    C = exp(-dis/r);
end