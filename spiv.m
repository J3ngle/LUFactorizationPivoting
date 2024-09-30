function [spL, spU] = spiv(A)
n = size(A, 1);
    % Initialize L and U
    spL = eye(n);
    spU = A;
    
    for k = 1:n
        for i = k+1:n
            if spU(k, k) == 0
                error('Zero pivot encountered. LU decomposition is not possible without pivoting.');
            end
            
            factor = spU(i, k) / spU(k, k);
            spL(i, k) = factor;
            spU(i, k:n) = spU(i, k:n) - factor * spU(k, k:n);
        end
    end
end
