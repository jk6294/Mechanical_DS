function [Q, W, v0, err] = construct_conic(Xs, Xsdot)
% Code for generating points along the solution space given specified nodes
% Able to construct conic in any dimension for any number of nodes
% 
% Inputs
% Xs:       d x n matrix of specified node positions
% Xsdot:    d x n matrix of specified node motions
%
% Outputs
% Q:    m+1 x m+1       conic matrix for the m-dimensional solution space
% W:    2d+1 x m        nullspace matrix of homogenous solutions
% v0:   2d+1 x 1        vector for the non-homogeneous solution
% err:  1 x 1           distance of particular solution to solution space

% Initialize Values and Matrices
d = size(Xs,1);                         % Dimensions
n = size(Xs,2);                         % Nodes
M = sym([Xsdot' Xs' -ones(n,1)]);       % Linearized matrix M from main
O = [zeros(d) eye(d) zeros(d,1);...     % O from supplement 
     eye(d)   zeros(d,d+1);...
     zeros(1,2*d+1)]/2;
b = sum(sym(Xs.*Xsdot))';
p = [zeros([2*d, 1]); 1];               % Inidicator vector

% Find homogeneous and non-homogenous solutions to linear problem
W = null(M);           % Nullspace gives homogeneous solutions
v0 = pinv(M)*b;        % Particular solution v* in main text
err = norm(M*v0 - b);  % Error in existence of partincular solution
if(err > 1e-12)
    disp('Error in reconstructing particular solution');
end

% Construct Conic Matrix by Applying Nonlinear Constraint
AP = W'*O*W;
BP = W'*(2*O*v0 - p)/2;
CP = (v0'*O - p')*v0;
Q = [AP,  BP;...
     BP', CP];

% Reconvert to Double
Q = double(Q);
W = double(W);
v0 = double(v0);
end