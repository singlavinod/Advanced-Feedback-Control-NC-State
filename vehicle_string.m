%% ECE 726      Project     Vinod K. Singla     11/28/2017
% Sparse feedback gains for string of moving vehicles
clear; clc; close all; % Housekeeping

% State-space representation of the vehicular formation with N vehicles
N = 10;
A = zeros(2*N - 1);
C = A;
Q = A;
B = zeros(2*N - 1, N);
i = 1;
j = 1;
k = 1;
l = 1;

while (i <= 2*N - 1)
    A(i,j) = -1;
    C(i,j) = 1;
    B(i,k) = 1;
    if i < 2*N - 1
        A(i+1,j) = 1;
        A(i+1,j+2) = -1;
        Q(i+1,j+1) = 10; % Can be increased to speed up response
    end
    i = i + 2;
    
    j = j + 2;
    k = k + 1;
    l = l + 1;
end

B1 = B; % Disturbance input matrix
B2 = B; % Control input matrix
R = eye(N);
rho = 100; 
n = 100; % Number of ADMM iterations
J = zeros(1,50); % Initialize objective function value
cardF = J; % Initialize cardinality
gamma = linspace(0, 2, 50); % Range of gamma values for promoting sparsity
fprintf("\t Gamma \t\t\t H2_norm \t Number of non-zero elements in F\n");

tic
for i = 1:50
    [F, J(i), cardF(i)] = sparselqr(A, B1, B2, Q, R, rho, n, gamma(i));
end
toc

% Plots

% Final sparsity pattern
figure;
spy(F, 15);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
xlabel('Sparsity pattern of F')

% Cardinality vs Gamma
figure;
yyaxis left;
plot(gamma, cardF,'-','LineWidth',2);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
xlab = xlabel('Gamma, \gamma','interpreter', 'tex');
ylab = ylabel('Cardinality(F)');
set(xlab, 'FontName', 'cmmi10', 'FontSize', 18)
set(ylab, 'FontName', 'cmmi10', 'FontSize', 18)

yyaxis right;
plot(gamma, cardF.*(100/(N*(2*N-1))),':','LineWidth',2);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
ylab = ylabel('Percent of non-zero elements');
set(ylab, 'FontName', 'cmmi10', 'FontSize', 18)
legend('Cardinality','% of non-zero elements');

% Objective function J(F) vs Gamma
figure;
[Ff, P] = lqr(A,B2,Q,R);
Jc = trace(P*(B1*B1'));
yyaxis left;
plot(gamma,J,'-','LineWidth',2);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
xlab = xlabel('Gamma, \gamma','interpreter', 'tex');
ylab = ylabel('Objective, J(F)');
set(xlab, 'FontName', 'cmmi10', 'FontSize', 18)
set(ylab, 'FontName', 'cmmi10', 'FontSize', 18)

yyaxis right;
plot(gamma,J.*(100/Jc),':','LineWidth',2);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
ylab = ylabel('Percent of original objective');
set(ylab, 'FontName', 'cmmi10', 'FontSize', 18)
legend('J(F)','% of orig. J(F)');

% Compare steady state performance
X0 = zeros(2*N - 1, 1); X0(5) = 2.1;
sys1 = ss((A - B2*Ff), B, C(1:N,:), 0); % For full-order feedback, Ff
sys2 = ss((A - B2*Ff), B, C(N+1:end,:), 0); % For full-order feedback, Ff
figure; initial(sys1, X0, 10);
figure; initial(sys2, X0, 10);

sys3 = ss((A - B2*F), B, C(1:N,:), 0); % For sparse feedback, F
sys4 = ss((A - B2*F), B, C(N+1:end,:), 0); % For sparse feedback, F
figure; initial(sys3, X0, 10);
figure; initial(sys4, X0, 10);