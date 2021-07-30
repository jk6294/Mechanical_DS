function [XC, fC] = sim_motion8(Xs, Xu, conn, delS, n, X0, pV)
% Function for simulating the motion of frames in 2-dimensions
% with one degree of freedom
% 
% Inputs
% Xs:     2 x ns        matrix of coordinates for specified nodes
% Xu:     2 x nu        matrix of coordinates for unspecified nodes
% conn:   k x 2         matrix of k edges connecting node i to node j
% delS:   1 x 1         instantaneous length of motion
% n:      1 x 1         number of simulation time steps
% X0:     2 x N         vector of desired initial direction of motion
% pV:     1 x 1         scalar to determine plotting. 0 for no plot
%
% Outputs
% XC:     2 X N X n     matrix of node motions across n time steps
% fC:     1 X n         matrix of motion errors given as energy


%% Lengths
ns = size(Xs,2);        % Number of specified nodes
nu = size(Xu,2);        % Number of unspecified nodes
X = [Xs Xu];            % Combined node positions
N = ns + nu;            % Number of total nodes
% All edge lengths
LVal = sqrt((X(1,conn(:,1)) - X(1,conn(:,2))).^2 +...
            (X(2,conn(:,1)) - X(2,conn(:,2))).^2)';
        
        
%% Symbolic 
% Define symbolic variables
XS = sym('x', [N, 1]); assume(XS, 'real');
YS = sym('y', [N, 1]); assume(YS, 'real');

% Construct rigidity matrix
I = zeros(size(conn,1), N);
for i = 1:size(conn,1)
    I(i,conn(i,1)) =  1;
    I(i,conn(i,2)) = -1;
end
J = [(XS(conn(:,1)) - XS(conn(:,2))) .* I,...
     (YS(conn(:,1)) - YS(conn(:,2))) .* I];
Jf = matlabFunction(J, 'Optimize', false, 'Sparse', true);

% Energy
E = sum((LVal - sqrt((XS(conn(:,1)) - XS(conn(:,2))).^2 + (YS(conn(:,1)) - YS(conn(:,2))).^2)).^2);
Ef = matlabFunction(E, 'Optimize', false);

% Rigid Body Motions
r1 = [ones([N,1]); zeros([N,1])];           % x-translation
r2 = [zeros([N,1]); ones([N,1])];           % y-translation
fR3 = @(xP, yP) [sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2)] .*...
                [-sin(atan2(yP, xP)); cos(atan2(yP, xP))];
            
% Error threshhold 
eTh = 1e-14;


%% RK4 Approximation
% Initial Conditions
xC = zeros([N, n]); xC(:,1) = X(1,:)';
yC = zeros([N, n]); yC(:,1) = X(2,:)';
fC = zeros([1, n]); fC(1) = 0; X0 = X0';
delX = zeros([2*N, n]); delX(:,1) = X0(:);

% Print Total Progress
fprintf([repmat('.',1,100) '\n\n']);
nTh = 0;

tic
nV = [];
for i = 2:n
    % Indicate Progress
    if(i/n >= nTh)
        fprintf(repmat('\b',1,length(nV)));
        fprintf('\b=\n');
        nV = num2str(floor(toc*(n-i)/i));
        fprintf(nV);
        nTh = nTh+.01;
    end
    
    % Step 1
    XP1 = xC(:,i-1); 
    YP1 = yC(:,i-1);
    XAC = num2cell([XP1; YP1]);
    R = [r1 r2 fR3(XP1, YP1)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k1 = double(V*null(R'*V));
    if(size(k1,2)>1)
        k1 = k1*(k1'*delX(:,i-1));
    end
    k1 = sign(k1'*delX(:,i-1)) * k1 / sqrt(k1'*k1);
    
    % Step 2
    XP2 = xC(:,i-1) + 1/4*delS*k1(1:N); 
    YP2 = yC(:,i-1) + 1/4*delS*k1((1:N)+N);
    XAC = num2cell([XP2; YP2]);
    R = [r1 r2 fR3(XP2, YP2)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k2 = double(V*null(R'*V));
    if(size(k2,2)>1)
        k2 = k2*(k2'*delX(:,i-1));
    end
    k2 = sign(k2'*delX(:,i-1)) * k2 / sqrt(k2'*k2);
    
    % Step 3
    XP3 = xC(:,i-1) + 3/32*delS*(k1(1:N)+3*k2(1:N)); 
    YP3 = yC(:,i-1) + 3/32*delS*(k1((1:N)+N)+3*k2((1:N)+N));
    XAC = num2cell([XP3; YP3]);
    R = [r1 r2 fR3(XP3, YP3)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k3 = double(V*null(R'*V));
    if(size(k3,2)>1)
        k3 = k3*(k3'*delX(:,i-1));
    end
    k3 = sign(k3'*delX(:,i-1)) * k3 / sqrt(k3'*k3);
    
    % Step 4
    XP4 = xC(:,i-1) + 12/2197*delS*(161*k1(1:N)-600*k2(1:N)+608*k3(1:N)); 
    YP4 = yC(:,i-1) + 12/2197*delS*(161*k1((1:N)+N)-600*k2((1:N)+N)+608*k3((1:N)+N));
    XAC = num2cell([XP4; YP4]);
    R = [r1 r2 fR3(XP4, YP4)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k4 = double(V*null(R'*V));
    if(size(k4,2)>1)
        k4 = k4*(k4'*delX(:,i-1));
    end
    k4 = sign(k4'*delX(:,i-1)) * k4 / sqrt(k4'*k4);
    
    % Step 5
    XP5 = xC(:,i-1) + 1/4104*delS*(8341*k1(1:N)-32832*k2(1:N)+29440*k3(1:N)-845*k4(1:N)); 
    YP5 = yC(:,i-1) + 1/4104*delS*(8341*k1((1:N)+N)-32832*k2((1:N)+N)+29440*k3((1:N)+N)-845*k4((1:N)+N)); 
    XAC = num2cell([XP5; YP5]);
    R = [r1 r2 fR3(XP5, YP5)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k5 = double(V*null(R'*V));
    if(size(k5,2)>1)
        k5 = k5*(k5'*delX(:,i-1));
    end
    k5 = sign(k5'*delX(:,i-1)) * k5 / sqrt(k5'*k5);
    
    % Step 6
    XP6 = xC(:,i-1) + delS*(-(8/27)*k1(1:N)+2*k2(1:N)-(3544/2565)*k3(1:N)+(1859/4104)*k4(1:N)-(11/40)*k5(1:N)); 
    YP6 = yC(:,i-1) + delS*(-(8/27)*k1((1:N)+N)+2*k2((1:N)+N)-(3544/2565)*k3((1:N)+N)+(1859/4104)*k4((1:N)+N)-(11/40)*k5((1:N)+N)); 
    XAC = num2cell([XP6; YP6]);
    R = [r1 r2 fR3(XP6, YP6)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k6 = double(V*null(R'*V));
    if(size(k6,2)>1)
        k6 = k6*(k6'*delX(:,i-1));
    end
    k6 = sign(k6'*delX(:,i-1)) * k6 / sqrt(k6'*k6);
    
    % Update
    delX(:,i) = 1/5*((16/27)*k1+(6656/2565)*k3+(28561/11286)*k4-(9/10)*k5+(2/11)*k6);
    xC(:,i) = xC(:,i-1) + delS * delX(1:N,i);
    yC(:,i) = yC(:,i-1) + delS * delX((1:N)+N,i);
    XAC = num2cell([xC(:,i); yC(:,i)]);
    fC(i) = Ef(XAC{:});
    
    if(fC(i) > eTh)
        disp('error over threshhold');
        break;
    end
end
toc

XC = zeros(2, N, n);
XC(1,:,:) = reshape(xC, [1, N, n]);
XC(2,:,:) = reshape(yC, [1, N, n]);


%% Plot
% Plot Parameters
ms = 4;         % Marker Size
lw = 1;         % Line Width
ea = .5;        % Edge Transparency
C_SN = [255 100 100]/255;       % Color of Specified Node
C_UN = [100 100 255]/255;       % Color of Unspecified Node

if(pV ~= 0)
    hold on;
    % Edges
    line([XC(1,conn(:,1),1); XC(1,conn(:,2),1)],...
         [XC(2,conn(:,1),1); XC(2,conn(:,2),1)],...
         'linewidth', lw, 'color', [0 0 0 ea]);
    % Node Movements Over Time
    for i = 1:n
        plot(XC(1,1:ns,i), XC(2,1:ns,i), 'r.', 'linewidth', ms, 'markersize', ms, 'color', C_SN)
        plot(XC(1,(1:nu)+ns,i), XC(2,(1:nu)+ns,i), 'b.', 'linewidth', ms, 'markersize', ms, 'color', C_UN);
    end
    % Specified Nodes
    plot(XC(1,1:ns,1), XC(2,1:ns,1), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN)
    % Unspecified Nodes
    plot(XC(1,(1:nu)+ns,1), XC(2,(1:nu)+ns,1), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN);
    hold off;
    set(gca,'visible',0);
    set(gcf,'color','w');
end

end