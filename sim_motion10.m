function [XC, fC] = sim_motion10(Xs, Xu, conn, delS, n, X0, pV)
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

%% Coefficients for RK10
% ak
p = load('RK10_ak.txt');
ak = zeros(size(p,1),1);
ak(p(:,1)+1) = p(:,2);
% bk
p = load('RK10_ck.txt');
ck = zeros(size(p,1),1);
ck(p(:,1)+1) = p(:,2);
% beta_jk
p = load('RK10_Bkj.txt');
Bkj = zeros(max(p(:,1))+1);
for i = 1:size(p,1)
    Bkj(p(i,1)+1,p(i,2)+1) = p(i,3);
end


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


%% RK10 Approximation
% Define variables
xC = zeros([2*N, n]); xC(:,1) = [X(1,:) X(2,:)]';   % Real states
fC = zeros([1, n]); fC(1) = 0; X0 = X0';            % Error eval
delX = zeros([2*N, n]); delX(:,1) = X0(:);          % Derivative direction
XK = zeros([2*N, length(ak)]);                      % Intermediary values

% Indicate Progress
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

    XK(:,1) = xC(:,i-1);
    F = zeros([2*N,length(ak)]);
    for k = 2:length(ak) + 1
        % Evaluate derivative
        R = [r1 r2 fR3(XK(1:N,k-1), XK(N+1:end,k-1))];
        XAC = num2cell(XK(:,k-1));
        K = full(Jf(XAC{:}));
        V = null(K);
        kp = double(V*null(R'*V));
        if(size(kp,2)>1)
            kp = kp*(kp'*delX(:,i-1));
        end
        F(:,k-1) = sign(kp'*delX(:,i-1)) * kp / sqrt(kp'*kp);
        
        if(k <= length(ak))
            FB = zeros([2*N,1]);
            for j = 1:k
                % Evaluate derivative
                FB = FB + Bkj(k,j)*F(:,j);
            end
            XK(:,k) = XK(:,1) + delS*FB;
        end
    end
    
    for k = 1:length(ak)
        delX(:,i) = delX(:,i) + ck(k)*F(:,k);
    end
    xC(:,i) = xC(:,i-1) + delS*delX(:,i);
    
    % Evaluate energy
    XAC = num2cell(xC(:,i));
    fC(i) = Ef(XAC{:});
end
fprintf('\b');
XC = zeros(2, N, n);
XC(1,:,:) = reshape(xC(1:N,:), [1, N, n]);
XC(2,:,:) = reshape(xC(N+1:end,:), [1, N, n]);


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