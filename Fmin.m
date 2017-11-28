%% ECE 726      Project     Vinod K. Singla     11/27/2017
% The code presented here will refer to equations and terminlogies given in
% [1] extensively to ensure completeness while maintaining brevity. The
% function Fmin performs minimization of eqn. (12a) using Anderson-Moore
% method.

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
% Output, minimizing feedback gain, F

% [1]   Lin, Fu, et al. “Design of Optimal Sparse Feedback Gains via the
%       Alternating Direction Method of Multipliers.” IEEE Transactions on
%       Automatic Control, vol. 58, no. 9, 2013, pp. 2426–2431.,
%       doi:10.1109/tac.2013.2257618

function F = Fmin(A, B1, B2, Q, R, U, rho, n, tol)

% Initial feedback gain, F0
F0 = lqr(A,B2,Q,R);

% Solve eqn. (NC-L)
L = lyap(A - B2*F, B1*B1');

% Solve eqn. (NC-P)
P = lyap((A - B2*F)', Q + F'*R*F);

% Initial value of objective function, phi(F)
phi = trace(L*(Q + F'*R*F)) + (rho/2)*norm(F - U, 'fro')^2;

% Armijo rule variables
alpha = 0.3;
beta = 0.5;

for i = 1:n
    
    % Initial stepsize, s for line search
    step = 1;
    
    % Solve sylvester eqn. (NC-F) for F_bar
    F_bar = lyap(rho*inv(R), 2*L, -inv(R)*(2*B2'*P*L + rho*U));
    
    % Descent direction of phi(F), F_tilda
    F_tilda = F_bar - F;
    
    % Gradient of phi
    del_phi = 2*(R*F - B2'*P)*L + rho*(F - U);
    
    % Check if stopping criteria is met or not
    if norm(del_phi, 'fro') < tol
        break;
    end
    
    % Verify if F_tilda is a descent direction
    if trace(F_tilda'*del_phi) > 1.e-9
        error('F_tilda is not a descent direction');
    end
    
    while True
        
        F_update = F + step*F_tilda;
        
        % Extract the closed loop eigen values
        Acl_eig = real(eig(A - B2*F_update));
        
        % Check if F_update is stabilizing or not
        if max(Acl_eig) > 0
            % Set objective function to infinity if Acl is not Hurwitz
            phi_update = NaN; % infinity is defined as NaN aka not-a-number
        else
            % Calculate values for the update step
            L_update = lyap(A - B2*F_update, B1*B1');
            P_update = lyap((A - B2*F_update)', Q + F_update'*R*F_update);
            phi_update = trace(L_update*(Q + F_update'*R*F_update)) + ...
                (rho/2)*norm(F_update - U, 'fro')^2;
        end
        
        % Use line search method (Armijo rule), for determining step size
        % while maintaining closed loop stability and convergence to a
        % stationary point of objective function, phi
        if ((phi - phi_update) > step*alpha*trace( - F_tilda*del_phi))...
                && ~isnan(phi_update)
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
    phi = phi_update;
end

if (i == n)
    disp('Max iterations reached for Anderson-Moore method');
    disp('Gradient norm = '); sprintf('%.4f',norm(del_phi, 'fro'));
    
end