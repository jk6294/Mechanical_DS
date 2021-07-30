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

% Simulate network
% Initial conditions
Xs0 = [-1 -1  0  0;...
        0  0  0  0];

% Simulate networks
[XC ,fC ] = sim_motion10(Xs, [],conn, .01,77 ,Xs0 ,0);

% Translate networks to align with origin
for i = 1:size(XC,3)
    XC(:,:,i) = XC(:,:,i) - XC(:,4,i);
end

% Simulation errors
disp(['mean simulation error: ' num2str(mean(fC))]);

% Calculate distances between nodes
D = squeeze(sqrt(sum(diff(XC,1,2).^2)));                % 1 unit
D = D([1,3],:);
DLin = linspace(sqrt(2)-.01,2+.01,nT);