%Complete Pivot
function [L, U, P, Q] = completePivotingLU(A)
    n = size(A, 1);
    L = eye(n);
    U = A;
    P = eye(n);
    Q = eye(n);

    for k = 1:n
        % Find the indices of the pivot element (max absolute value)
        [maxVal, idx] = max(abs(U(k:n, k:n)));
        [pivot_row, pivot_col] = ind2sub(size(idx), idx);
        pivot_row = pivot_row + k - 1;
        pivot_col = pivot_col + k - 1;

        % Swap rows and columns in U
        U([k, pivot_row], k:n) = U([pivot_row, k], k:n);
        U(k:n, [k, pivot_col]) = U(k:n, [pivot_col, k]);

        % Swap rows in L (below the main diagonal)
        if k > 1
            L([k, pivot_row], 1:k-1) = L([pivot_row, k], 1:k-1);
        end

        % Swap columns in Q
        Q(:, [k, pivot_col]) = Q(:, [pivot_col, k]);

        % Swap rows in P
        P([k, pivot_row], :) = P([pivot_row, k], :);

        % Perform elimination
        for i = k+1:n
            factor = U(i, k) / U(k, k);
            L(i, k) = factor;
            U(i, k:n) = U(i, k:n) - factor * U(k, k:n);
        end
    end
end
