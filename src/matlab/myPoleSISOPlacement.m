function K = myPoleSISOPlacement(A, B, poles)
    % Verify that the system is controllable, or poles is wrong or not.

    controllabilityMatrix = ctrb(A,B);
    r = rank(controllabilityMatrix);
    n = size(A,1);
    
    % if we failed to control:
    if r ~= n
        %[Abar, Bbar, Cbar, T, k] = ctrbf(A,B,C);
        [T, A_c, B_c, A_uc, B_uc] = controllabilityDecomposition(A,B);
        %C_ = B';
        %[A_bar, B_bar, C_bar, T, k] = ctrbf(A,B,C_);
        
        %T
        %A_bar
        %B_bar
        
        %n_controllable = sum(k)
        %A_c = A_bar(n-n_controllable+1:n, n-n_controllable+1:n)
        %B_c = B_bar(n-n_controllable+1:n)
        % todo: deal with this. Decomposite.
        % A11 A12
        % 0   A22
        % A_dot = inv(controllabilityMatrix)*A*controllabilityMatrix
        % A = A_dot(1:r,1:r)
        % B = int(controllabilityMatrix)*B
        %error('The system is not controllable and poles cannot be placed as desired.');
    end
    
    % Get the size of the state matrix
    % n = r;
    
    if r ~= n
        K = [myBasePolePlacement(A_c,B_c,poles), zeros(1,n-r)];
        % K=place(A_c,B_c,poles);
        % K = K(end:-1:1);
        %K = [zeros(n-r), K]
        %K =K*T;
        K=K*inv(T);
    else
        K = myBasePolePlacement(A,B,poles);
    end
end

function K = myBasePolePlacement(A, B, poles)
    controllabilityMatrix = ctrb(A,B);
    n = size(A,1);
     % Construct the desired characteristic polynomial from the poles
    charPoly = poly(poles);
    originalPoly = poly(A);
    K_ccf = (charPoly - originalPoly);
    %K_ccf = zeros(1,n)
    K_ccf = K_ccf(end:-1:2);
    % perfect test!
    P_ccf_inv = zeros(n,n);
    for i=1:n-1
        for j = 1:n-i
            P_ccf_inv(i,j)=originalPoly(n-i-j+2);
        end
        P_ccf_inv(i,n-i+1)  = 1;
    end
    P_ccf_inv(n,1)=1;

    T_ccf = controllabilityMatrix * P_ccf_inv;
    % Calculate the gain matrix K using Ackermann's formula

    %K_ccf
    %P_ccf_inv
    %T_ccf
    %if()
    K = K_ccf*(T_ccf^(-1));
    % K = K(end:-1:1);

end

function [T, A_c, B_c, A_uc, B_uc] = controllabilityDecomposition(A, B)
    % Calculate the controllability matrix
    
    
    C = ctrb(A, B);
    
    % Determine the rank of the controllability matrix
    rankC = rank(C);
    
    % Select rankC independent columns from C
    [U, S, V] = svd(C, 'econ');
    independent_cols = V(:, 1:rankC);
    
    % Ensure we have a full set of basis vectors for the state space
    % If system is not fully controllable, complete the basis.
    if rankC < size(A,1)
        % Add additional vectors to span the entire space
        % Find the null space of the controllability matrix
        null_vectors = null(C');
        
        % Construct T by combining the independent columns and the null space
        T = [independent_cols, null_vectors];
    else
        % System is fully controllable
        T = independent_cols;
    end
    
    % Verify that T is full rank (non-singular)
    if rank(T) < size(A,1)
        error('Transformation matrix T is singular. The system may not be controllable.');
    end
    
    % Transform the system matrices
    T_inv = inv(T);
    A_tilde = T_inv * A * T;
    B_tilde = T_inv * B;
    
    % Extract the submatrices for the controllable part
    A_c = A_tilde(1:rankC, 1:rankC);
    B_c = B_tilde(1:rankC, :);
    
    % Extract the submatrices for the uncontrollable part, if any
    A_uc = A_tilde(rankC+1:end, rankC+1:end);
    B_uc = B_tilde(rankC+1:end, :);
end
