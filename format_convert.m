function format_convert(A)
 ind = 1:size(A,2);
    for i = 1:size(A,1)
        fprintf(' & %.3f ', A(i,ind));
        fprintf('\n');
    end
end