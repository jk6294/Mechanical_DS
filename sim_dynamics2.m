function [XC, XV, t, KE, PE] = sim_dynamics2(X, Xdot, conn, M, K, T, NV)

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


%% Functions
XS = sym('x', [N, 1]); assume(XS, 'real');
YS = sym('y', [N, 1]); assume(YS, 'real');

% Edge Length Function: Input: dN x 1. Output: k x 1.
fL = @(x) sqrt(sum(reshape((x(conn(:,1)) - x(conn(:,2))).^2,[m,d]),2));

% Spring Equilibrium Length
Leq = fL(X);                                    % k x 1

% Kinetic Energy Function: Input: dN x 1. Output; 1 x 1
fE = @(xdot) 0.5 * sum(xdot.^2 .* M);

% Potential Function: Input: dN x 1. Output: 1 x 1
fV = @(x) 0.5 * sum(K.*(Leq - fL(x)).^2);   

% Potential Energy Function
sU = 0.5*sum(K.*(Leq-fL([XS;YS])).^2);

% Potential Gradient
sdU = jacobian(sU,[XS;YS])';
fdU = matlabFunction(sdU, 'Optimize', false, 'Vars', {[XS;YS]});

% State Evolution
fA = @(y) [y([1:d*N]+d*N); -fdU(y(1:d*N))./M];


%% Simulation
% Simulation Parameters
y0 = [X; Xdot];
tspan = [0, T];
options = odeset('RelTol',1e-10, 'AbsTol', 1e-10);

% Simulate
fA
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
