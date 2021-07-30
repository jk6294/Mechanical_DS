%% Main module
% L:                            side length:            1 x 1
% Xs:                           start coordinates       2 x 4
% conn:                         connectivity            4 x 2
% Xsa:                          tesselated unit
% conna:                        tesselated connectivity
% XC:                           unit trajectory
% XCa:                          network trajectory
% D:                            unit distances
% Da:                           network distances
% DLin:                         color interpolation

L = 1;
% Single unit
Xs = [0  0  0  2*L;...
      L -L  0  0];
conn = [1 3; 1 4; 2 3; 2 4];
% Double unit
Xsc = [0  0  0  2*L  L  L;...
       L -L  0  0    0 -2*L];
connc = [1 3; 1 4; 2 3; 2 4; 3 5; 3 6; 4 5; 4 6];
Rz = rotz(45); Rz = Rz(1:2,1:2);
% Tesselate
[Xsa, conna] = tesselate_network(Rz*Xsc, connc, [sqrt(2);0], [6;1]);

% Simulate network
% Initial conditions
Xs0 = [-1 -1  0  0;...
        0  0  0  0];
Xsa0 = -Xsa; Xsa0(1,end-1) = 1;

% Simulate networks
[XC ,fC ] = sim_motion10(Xs, [],conn, .01,77 ,Xs0 ,0);
[XCa,fCa] = sim_motion10(Xsa,[],conna,.050045,150,Xsa0,0);

% Translate networks to align with origin
for i = 1:size(XC,3)
    XC(:,:,i) = XC(:,:,i) - XC(:,4,i);
end
for i = 1:size(XCa,3)
    XCa(:,:,i) = XCa(:,:,i) - XCa(:,1,i);
end

% Simulation errors
disp(['mean simulation error: '...
     [num2str(mean(fC)) '  ' num2str(mean(fCa))]]);

% Calculate distances between nodes
D = squeeze(sqrt(sum(diff(XC,1,2).^2)));                % 1 unit
D = D([1,3],:);
Da = squeeze(sqrt(sum(diff(XCa,1,2).^2)));              % CHain
Da = Da(1:2:end,:);
DLin = linspace(sqrt(2)-.01,2+.01,nT);