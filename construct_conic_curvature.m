function [Q1, Q2, W, v0, err] = construct_conic_curvature(Xs, dXs, ddXs)
% Code for generating points along the solution space given specified 
% nodes, and infinitesimal or finite motions. If z is not 0 or 1,
% default to infinitesimal
% 
% Inputs
% Xs:       2 x 3 matrix of specified node initial positions
% dXs:      2 x 3 matrix of specified node motions or final positions
% dXs:      2 x 3 matrix of specified node motions or final positions
%
% Outputs
% Q:    m+1 x m+1       conic matrix for the m-dimensional solution space
% W:    2d+1 x m        nullspace matrix of homogenous solutions
% v0:   2d+1 x 1        vector for the non-homogeneous solution
% err:  1 x 1           distance of particular solution to solution space

% Initialize Values and Matrices
d = 2;                          % Dimensions
n = 3;                          % Nodes
p1 = [zeros([3*d, 1]); 1; 0];   % Inidicator vector
p2 = [zeros([3*d, 1]); 0; 1];   % Inidicator vector

% Construct matrices
M = [[dXs'   Xs'     zeros(3,2) -ones(3,1)  zeros(3,1)];...
     [ddXs'  2*dXs'  Xs'        zeros(3,1)  -ones(3,1)]];
O1 = zeros(8); 
O1(1:2,3:4) = eye(2)/2; 
O1(3:4,1:2) = eye(2)/2;
O2 = zeros(8); 
O2(3:4,3:4) = eye(2); 
O2(1:2,5:6) = eye(2)/2; 
O2(5:6,1:2) = eye(2)/2;
b = sum([(Xs.*dXs)  dXs.^2+Xs.*ddXs])';

% Find particular solution
v0 = pinv(M)*b;
err = norm(M*v0-b);

% Find homogeneous solution
W = null(M);
W = rref(W')';
% W = fliplr(flipud(rref(flipud(W)')'));
% Solution 1
AP1 = W'*O1*W;
BP1 = W'*(2*O1*v0 - p1)/2;
CP1 = (v0'*O1 - p1')*v0;
Q1 = [AP1,  BP1;...
      BP1', CP1];
% Solution 2
AP2 = W'*O2*W;
BP2 = W'*(2*O2*v0 - p2)/2;
CP2 = (v0'*O2 - p2')*v0;
Q2 = [AP2,  BP2;...
      BP2', CP2];
end