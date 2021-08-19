%% Define unit
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
[Xscc,Xucc,conncc,CSSc] = network_chain_x(Xs1(:,1:2),Xuc,IC',...
                                          interp1(DLin2,CP2,Dc));


%% Simulate quadrifolium
Ns = 114;               % Number of non-added nodes in quad example

XCc0 = zeros(2,size([Xscc Xucc],2));
XCc0(1,1) = -1;
[XCcc,fCcc] = sim_motion10b(Xscc,Xucc,conncc,1.5,1157,XCc0,0);
disp(['mean simulation error: ' num2str(mean(fCcc))]);

% Correct rotation
XCcc = XCcc - XCcc(:,Ns,:);
XCcc(1,:,:) = -XCcc(1,:,:);
for i = 1:size(XCcc,3)
    Rz = rotz(atan2d(diff(XCcc(2,[-2 0]+Ns,i)),diff(XCcc(1,[-2 0]+Ns,i))));
    XCcc(:,:,i) = Rz(1:2,1:2)'*XCcc(:,:,i);
end

% Uncomment to save
save('quadrifolium.mat','Xscc','Xucc','conncc','XCcc','fCcc','CSSc');