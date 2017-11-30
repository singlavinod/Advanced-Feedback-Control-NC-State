%% ECE 726      Project     Vinod K. Singla     11/27/2017
% The code presented here will refer to equations and terminlogies given in
% [1] extensively to ensure completeness while maintaining brevity. The
% function sparselqr solves eqn. (SP) using Alternating Direction Method
% of Multipliers (ADMM) and then solves (SH2) using identified sparsity
% pattern
%   
%               minimize J(F) + gamma*g(F)              (SP)
%
%               minimize J(F)                           (SH2)
%
%    subject to:
%
%               F belonging to identified sparsity pattern

% Inputs:

% A = System matrix
% B1 = Disturbance input matrix
% B2 = Control input matrix
% Q = Weighting matrix for states
% R = Weighting matrix for inputs
% lambda = lagrange multiplier
% rho = a positive scalar
% G = defined by eqn. (12b)
% U = G - lambda/rho
% n = max iterations
% tol = tolerance for stopping criteria on gradient of phi
% Gamma = sparsity promoting parameter

% Outputs:

% F = optimal feedback gains found through solution of eqn. (SP)
% J = optimal value of LQR objective function using F
% N = number of non-zero elements in F

% [1]   Lin, Fu, et al. ?Design of Optimal Sparse Feedback Gains via the
%       Alternating Direction Method of Multipliers.? IEEE Transactions on
%       Automatic Control, vol. 58, no. 9, 2013, pp. 2426?2431.,
%       doi:10.1109/tac.2013.2257618

function [F, J, N] = sparselqr(A, B1, B2, Q, R, rho, n, Gamma)

% Extract input matrix dimensions
[i, j] = size(B2);

% Set the initial value of G to be LQR centralized gain, also equal to
% initial value of F
F = lqr(A, B2, Q, R);
G = F;

% Initialize lagrange multipliers to zero
lambda = zeros(j, i);

% Define tolerances for convergence
tol_ADMM_abs = 1e-4; % absolute tolerance for ADMM
tol_ADMM_rel = 1e-2; % relative tolerance for ADMM
tol_AM = 1e-2; % tolerance for Anderson-Moore method
tol_NT = 1e-3; % tolerance for Newton's method of conjugate gradient

% Define number of Fmin iterations for each k
n_AM = 100;

for k = 1: n
    
    % Perform F minimization using Fmin (Anderson-Moore method)
    U = G - lambda/rho;
    F = Fmin(A, B1, B2, Q, R, U, rho, n_AM, tol_AM);
    
    % Perform G minimization
    V = F + lambda/rho;
    
    % G_new is calculated assuming penalty is cardinality function
    b = sqrt(2 * Gamma / rho);
    G_new = V .* ( abs(V) > b );
    
    %%%%%%%%%%%%%%%%%%% Convergence check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate residuals to check for convergence
    FG_res = norm(F - G_new, 'fro');
    G_res = norm(G - G_new, 'fro');
    
    % Calculate residual stopping criteria
    FG_lim = sqrt(i*j) * tol_ADMM_abs + ...
        tol_ADMM_rel * max(norm(F,'fro'), norm(G_new,'fro'));
    G_lim = sqrt(i*j) * tol_ADMM_abs + tol_ADMM_rel * norm(lambda,'fro');
    
    % Check if stopping criteria is met
    if (FG_res < FG_lim) && (rho*G_res < G_lim)
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update G
    G = G_new;
    
    % Update lagrange multipliers
    lambda = lambda + rho*(F - G);
    
end

% Check if G is stabilizing
if max(real(eig(A - B2*G))) < 0
    F = G;
end

% Number of non-zero elements in F
N = nnz(F);

% Solve the structured H2 norm problem, (SH2)
S = double( F ~= 0 );
[F, J] = solve_SH2(A, B1, B2, Q, R, S, F, tol_NT);
fprintf("\t %.4f \t\t %.2f \t\t %d \n", Gamma, J, N);
end