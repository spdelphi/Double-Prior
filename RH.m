function [S] = RH(x,lambda_s)
[L,N] = size(x);
T = diag(x*x');
S = zeros(size(x));
for i = 1:L
    for j = 1:N
        if T(i) > lambda_s
            S(i,j) = x(i,j);
        else
            S(i,j) = 0;
        end
    end
end

end