%Partial Pivot Function
function [L, U, P] = partialPivoting(A)
    n = size(A, 1);
    L = eye(n);
    U = A;
    P = eye(n);

    for k = 1:n
        % Partial Pivoting: Find the pivot row
        [~, pivot_row] = max(abs(U(k:n, k)));
        pivot_row = pivot_row + k - 1;

        % Swap rows in U
        U([k, pivot_row], k:n) = U([pivot_row, k], k:n);

        % Swap rows in L (below the main diagonal)
        if k > 1
            L([k, pivot_row], 1:k-1) = L([pivot_row, k], 1:k-1);
        end

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