%% ECE 726      Project     Vinod K. Singla     11/28/2017
% Sparse feedback gains for string of moving vehicles

% State-space representation of the vehicular formation with N = 50 masses
N = 50;
I = eye(N,N);
Z = zeros(N,N);
T = toeplitz([2 -1 zeros(1,N-2)]);
A = [Z I; -T Z]; 
B1 = [Z; I]; 
B2 = [Z; I]; 
Q = eye(2*N); 
R = 10*I;