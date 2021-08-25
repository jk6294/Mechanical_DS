%% Define unit
params_fig;
sq = sqrt(3)/2;
L = 1;                                          % Initial length
LF = 1.7;                                       % Final length
ds1 = sqrt(3);
ds2 = LF*sqrt(3);

Xs1 = [-sq 0  sq;...
      -.5 1 -.5]*L;
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Distance ranges
DLin = linspace(2*sq-.01,ds2+.01,nT);           % Distance range, d1,d2
DLin2 = linspace(1.52,3.97,nT);                 % Distance range, c


%% Tesselate quadrifolium
% Parameters
a = 16.2;                                   % Quadrifolium radius
f = @(tho)[a*sin(2*tho) .* cos(tho);...     % Quadrifolium equation
           a*sin(2*tho) .* sin(tho)];

% Trace quadrifolium
nv = 1000;
thv = linspace(0,2.1*pi,nv) + pi/4;         % Sampled angles
xv = f(thv);                                % Sampled points

% Initialize 
xP = zeros(2,113);                          % Triangle tesselation points
xP(:,1) = f(thv(1))*1.084543;               % Initial triangle point
thop = thv(1);

% Tesselate quadrifolium with isoscelese triangles
options = optimset('TolFun',1e-15,'TolX',1e-15);
for i = 1:size(xP,2)
    [thop, fV] = fzero(@(t) (sqrt(sum((f(t)-xP(:,i)).^2))-ds2/2),...
                             thop+[.01,.2],options);
    xP(:,i+1) = xP(:,i) + 2*(f(thop)-xP(:,i));
end

% Bin triangle distances
Dc = sqrt(sum((xP(:,1:end-2)-xP(:,3:end)).^2));     % Triangle distances: c
[C,~,IC] = unique(round(Dc,3),'stable');            % Bin to 3rd decimal
Xuc = zeros(2,2,length(C));                         % Unique node positions
Xucf = zeros(2,2,length(C));
Xfp = zeros(2,3,length(C));

% Construct units for each unique bin
Xu0a = [-3.5*sq 3.5*sq; 2.2 2.2];                   % Initial conditions
Xu0b = [-sq+.1 -sq+.1; -2 2.5];
Xu0c = [-sq*2 sq*2;  -3.3 -3.3];
for i = 1:length(C)
    % Different initial conditions based on distance c
    if(C(i)<2.2)
        Xu0 = Xu0a;
    elseif(C(i)>4)
        Xu0 = Xu0c;
    else
        Xu0 = Xu0b;
    end
    Xfp(:,:,i) = [[-C(i)/2; -sqrt(ds2^2 - (C(i)/2)^2)] [0;0]...
                  [C(i)/2; -sqrt(ds2^2 - (C(i)/2)^2)]] + [0;Xs1(2,2)];
    Xup = construct_network(Xs1,Xfp(:,:,i),Xu0,conn,0,1);
    Xuc(:,:,i) = Xup(1:2,:);
    Xucf(:,:,i) = Xup(3:4,:);
end

% Construct chain for quadrifolium network at the start position
[Xscca,Xucca,conncc,CSSc] = network_chain_x(Xs1(:,1:2),Xuc,IC',...
                                          interp1(DLin2,CP2,Dc));
                                      
Ns = size(Xscca,2);               % Number of non-added nodes in quad example
XCc0 = zeros(2,size([Xscca Xucca],2));
XCc0(1,1) = -1;
                                      

%% For loop iteration
nV = 0.00022;
nIt = 200;

load quad_disorder2.mat
% XCcc = zeros([size(XCc0) 2000 nIt]);
% fCcc = zeros(2000,nIt);
% parfor i = 1:nIt
%     disp(i);
%     % Add noise
%     Xscc = Xscca + (rand(size(Xscca))-.5)*nV*i;
%     Xucc = Xucca + (rand(size(Xucca))-.5)*nV*i;
%     [XCccp,fCccp,~,~] = sim_motion10(Xscc,Xucc,conncc,1,2000,XCc0,0);
%     XCcc(:,:,:,i) = XCccp;
%     fCcc(:,i) = fCccp;
% end
% save('quad_disorder.mat', 'XCcc','fCcc');


%% Distribution of bond lengths
% Center trace
xPc = xP - mean(xP,2);

% Iterate
nIter = size(XCcc,4);
nStep = size(XCcc,3);
BLena = sqrt(sum((Xscca(:,conncc(:,1)) - Xucca(:,conncc(:,2)-Ns)).^2));
DMat = zeros([nIter, nStep]);
BLenM = zeros([nIter,size(conncc,1)]);

for j = 1:nIter
    xCcc = XCcc(:,:,:,j);
    BLenM(j,:) = sqrt(sum((xCcc(:,conncc(:,1),1) - xCcc(:,conncc(:,2),1)).^2));
    for i = 1:nStep
        xprox = xCcc(:,1:Ns,i) - mean(xCcc(:,1:Ns,i),2);
        [u,s,v] = svd(xPc * xprox');
        Rp = u*v';
        thp = atan2d(Rp(2,1),Rp(2,2));
        Rz = rotz(thp); Rz = Rz(1:2,1:2);
        XDiff = Rz * xprox - xPc;
        DMat(j,i) = mean(sqrt(sum(XDiff.^2)));
    end
end

BLenDel = mean(abs((BLenM - BLena) ./ BLena),2);
DMatMin = min(DMat,[],2);