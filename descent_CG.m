%% ECE 726      Project     Vinod K. Singla     11/27/2017
% The code presented here will refer to equations and terminlogies given in
% [1] extensively to ensure completeness while maintaining brevity.

% Inputs:

% A = System matrix
% B2 = Control input matrix
% F = initial feedback gain
% P = eqn. (NC-P)
% L = controllability gramain, eqn. (NC-L)
% R = Weighting matrix for inputs
% del_J = objective function gradient
% S = Identified sparsity structure

% Outputs:

% F_cj = descent direction found through Newton's method of conjugate grad

% [1]   Lin, Fu, et al. ?Design of Optimal Sparse Feedback Gains via the
%       Alternating Direction Method of Multipliers.? IEEE Transactions on
%       Automatic Control, vol. 58, no. 9, 2013, pp. 2426?2431.,
%       doi:10.1109/tac.2013.2257618

function F_cj = descent_CG(A, B2, F, P, L, R, del_J, S)

% Number of nonzero elements
q = nnz(S);
Pi = del_J;
Delta = - del_J;
F_tilda = zeros(size(S));
ndel_J = norm(del_J, 'fro');

% Closed-loop A-matrix
Acl = A - B2 * F;

% Denote Z  = R * F - B2' * P
Z = R * F - B2' * P;

% Conjugate gradient scheme
for k = 0 : q
    
    G1 = B2 * Delta * L;
    G2 = -Z' * Delta;
    L_tilda = lyap( Acl, - G1 - G1' );
    P_tilda = lyap( Acl', - G2 - G2');
    H = 2*((( R * Delta - B2' * P_tilda) * L + Z * L_tilda) .* S);
    
    % Negative curvature test form the inner product between H and Delta
    trHDelta = sum(sum( H .* Delta));
    if (trHDelta <= 0) && (k == 0)
        F_cj = - del_J;
        break;
    elseif (trHDelta <= 0) && (k > 0)
        F_cj = F_tilda;
        break;
    end
    
    alpha  = - sum(sum( Pi .* Delta)) / trHDelta;
    F_tilda = F_tilda + alpha * Delta;
    Pi = Pi + alpha * H;
    beta = sum(sum(H .* Pi)) / trHDelta;
    Delta = -Pi + beta * Delta;
    
    % Nocedal and Wright '06 p168
    % Stopping criterion for the conjugate gradient method
    if norm(Pi, 'fro') < min(0.5, sqrt(ndel_J)) * ndel_J
        F_cj = F_tilda;
        break;
    end
end

if k == q
    F_cj = F_tilda;
end
end