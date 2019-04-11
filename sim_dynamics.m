function [XC, XV, t, KE, PE] = sim_dynamics(X, Xdot, conn, M, K, T, NV)
% Function for simulationg the motion of frames
% Inputs
%       X:          N x d matrix of node coordinates
%       Xdot:       N x d matrix of node velocities
%       conn:       m x 2 matrix of connections
%       M:          N x 1 vector of node masses
%       K:          m x 1 vector of spring constants
%       T:          1 x 1 scalar of end time
%       NV:         1 x 1 scalar for resample interpolation
%
% Outputs
%       XC:          N x 3 x T matrix of mode positions over T time steps
%       fC:         1 x T matrix of motion errors


%% Parameters
% Dimension
d = size(X,2);                                  % Dimensions
N = size(X,1);                                  % Nodes
m = size(conn,1);                               % Edges
o = ones(d*N,1);                                % Vector of 1

% Reshape Node Matrix to Vector and Extend Connectivity Vector
X = X(:);                                       % dN x 1
Xdot = Xdot(:);                                 % dN x 1
dL = repmat(([1:d]-1)*N, [m,1]);                % dk x 1
conn = repmat(conn,[d,1])+dL(:);                % dk x 1
M = repmat(M,[d,1]);

% Create Indication Matrix
S = ones(d*N)*(d*m+1);
for i = 1:(d*m)
    S(conn(i,1),conn(i,2)) = i;
    S(conn(i,2),conn(i,1)) = i;
end


%% Functions
% Edge Length Function: Input: dN x 1. Output: k x 1.
fL = @(x) sqrt(sum(reshape((x(conn(:,1)) - x(conn(:,2))).^2,[m,d]),2));

% Spring Equilibrium Length
Leq = fL(X);                                    % k x 1

% Potential Function: Input: dN x 1. Output: 1 x 1
fV = @(x) 0.5 * sum(K.*(Leq - fL(x)).^2);    

% Kinetic Energy Function: Input: dN x 1. Output; 1 x 1
fE = @(xdot) 0.5 * sum(xdot.^2 .* M);

% Potential Gradient: Scaling Factor: k x 1
fK = @(x) repmat(K.*(1-Leq./fL(x)),[d,1]);

% Selection Function
select = @(k,m) k(m);

% Potential Gradient: Input: dN x 1. Output: dN x 1
fdU = @(x) sum(select([fK(x); 0],S) .* (x*o' - o*x'),2);

% State Evolution
fA = @(y) [y([1:d*N]+d*N); -fdU(y(1:d*N))./M];


%% Simulation
% Simulation Parameters
y0 = [X; Xdot];
tspan = [0, T];
options = odeset('RelTol',1e-10, 'AbsTol', 1e-10);

% Simulate
[t, y] = ode45(@(t,y) fA(y), tspan, y0, options);
y = y';
% Interpolate
tq = linspace(0,T,NV);
y = interp1(t,y',tq)';
% Reshape
XC = reshape(y(1:d*N,:),[N,d,length(tq)]);
XV = reshape(y([1:d*N]+d*N,:),[N,d,length(tq)]);


disp('Iteration');
PE = zeros(1,length(tq));
KE = zeros(1,length(tq));
for i = 1:length(tq)
    PE(i) = fV(y(1:d*N,i));
    KE(i) = fE(y([1:d*N]+d*N,i));
end

end