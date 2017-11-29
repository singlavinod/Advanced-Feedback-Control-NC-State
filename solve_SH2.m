%% ECE 726      Project     Vinod K. Singla     11/27/2017
% The code presented here will refer to equations and terminlogies given in
% [1] extensively to ensure completeness while maintaining brevity.
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
% S = Identified sparsity structure
% tol = tolerance for stopping criteria on gradient of J

% Outputs:

% F = optimal feedback gains found through solution of eqn. (SH2)
% J = optimal value of LQR objective function using F

% [1]   Lin, Fu, et al. ?Design of Optimal Sparse Feedback Gains via the
%       Alternating Direction Method of Multipliers.? IEEE Transactions on
%       Automatic Control, vol. 58, no. 9, 2013, pp. 2426?2431.,
%       doi:10.1109/tac.2013.2257618

function [F, J] = solve_SH2(A, B1, B2, Q, R, S, F0, tol)
% Max iterations for Newton's method
n = 100;

% Check if initial condition is stabilizing
if max(real(eig(A - B2*F0))) > 0
    error('Closed loop system is unstable with F0');
end

% Initialize F
F = F0;

% Armijo rule variables
alpha = 0.3;
beta = 0.5;

% Solve eqn. (NC-L)
L = lyap(A - B2*F, B1*B1');

% Solve eqn. (NC-P)
P = lyap((A - B2*F)', Q + F'*R*F);

% Initial value of objective function, J
J = trace(L*(Q + F'*R*F));

for i = 1:n
    
    % Initial stepsize, s for line search
    step = 1;
    
    % Gradient of J
    del_J = 2*(R*F - B2'*P)*L.*S;
    
    % Check if stopping criteria is met or not
    if norm(del_J, 'fro') < tol
        break;
    end
    
    % Use conjugate gradient to compute descent direction
    F_cj = descent_CG(A, B2, F, P, L, R, del_J, S);
    
    while True
        
        F_update = F + step*F_cj;
        
        % Extract the closed loop eigen values
        Acl_eig = real(eig(A - B2*F_update));
        
        % Check if F_update is stabilizing or not
        if max(Acl_eig) > 0
            % Set objective function to infinity if Acl is not Hurwitz
            J_update = NaN; % infinity is defined as NaN aka not-a-number
        else
            % Calculate values for the update step
            L_update = lyap(A - B2*F_update, B1*B1');
            P_update = lyap((A - B2*F_update)', Q + F_update'*R*F_update);
            J_update = trace(L_update*(Q + F_update'*R*F_update));
        end
        
        % Use line search method (Armijo rule), for determining step size
        % while maintaining closed loop stability and convergence to a
        % stationary point of objective function, phi
        if ((J - J_update) > step*alpha*trace( - F_cj'*del_J))...
                && ~isnan(J_update)
            break;
        end
        
        % Update stepsize
        step = step*beta;
        
        % Error if stepsize is too small
        if step < 1.e-14
            error('step size is too small');
        end
    end
    
    % update values
    F = F_update;
    L = L_update;
    P = P_update;
    J = J_update;
end

if (i == n)
    disp('Max iterations reached for Newton"s method');
    disp('Gradient norm = '); sprintf('%.4f',norm(del_J, 'fro'));
end

end