function [Lcomp, Ucomp, Pcomp, Qcomp] = completePivotingLU(A)
    n = size(A, 1);
    Lcomp = eye(n);
    Ucomp = A;
    Pcomp = eye(n);
    Qcomp = eye(n);

    for k = 1:n
        % Find the indices of the pivot element (max absolute value)
        [maxVal, idx] = max(abs(Ucomp(k:n, k:n)));
        [pivot_row, pivot_col] = ind2sub(size(idx), idx);
        pivot_row = pivot_row + k - 1;
        pivot_col = pivot_col + k - 1;

        % Swap rows and columns in U
        Ucomp([k, pivot_row], k:n) = Ucomp([pivot_row, k], k:n);
        Ucomp(k:n, [k, pivot_col]) = Ucomp(k:n, [pivot_col, k]);

        % Swap rows in L (below the main diagonal)
        if k > 1
            Lcomp([k, pivot_row], 1:k-1) = Lcomp([pivot_row, k], 1:k-1);
        end

        % Swap columns in Q
        Qcomp(:, [k, pivot_col]) = Qcomp(:, [pivot_col, k]);

        % Swap rows in P
        Pcomp([k, pivot_row], :) = Pcomp([pivot_row, k], :);

        % Perform elimination
        for i = k+1:n
            factor = Ucomp(i, k) / Ucomp(k, k);
            Lcomp(i, k) = factor;
            Ucomp(i, k:n) = Ucomp(i, k:n) - factor * Ucomp(k, k:n);
        end
    end
end
