function pivdecision = decision(A)

    % Check for symmetric matrix (for complete pivoting)
    sym = isequal(A, A');

    % Check for diagonally dominant matrix (for no pivoting)
    diagdom = all(2 * abs(diag(A)) >= sum(abs(A), 2) - 2 * abs(diag(A)));

    % Decide the pivoting strategy based on the above checks
    if sym
        pivdecision = 'Complete Pivoting';
    elseif diagdom
        pivdecision = 'No Pivoting';
    else
        pivdecision = 'Partial Pivoting';
    end

    % Display the recommendation
    fprintf('Use %s for LU decomposition.\n', pivdecision);
end
