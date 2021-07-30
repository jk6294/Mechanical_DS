%% Panel b: 4-node linkage
% Xf:                           end coordinates:    2 x 4
% conn:                         node pairs:         4 x 2
% L:                            bar length:         1 x 1

L = 1;
Xf = [0  0  0  2;...
      1 -1  0  0];
conn = [1 3; 1 4; 2 3; 2 4];


%% Panel c: rigid body motions
% sq:                           lengths
% Xs:                           start position
% XC:                           trajectory
% fC:                           error
% D:                            inter-node distances
% DLin:                         distance interpolation

sq = sqrt(2)*L/2;
Xs = [-sq -sq 0 sq*2;...
       sq -sq 0 0];
   
% Simulate network
[XC,fC] = sim_motion10(Xs,[],conn,.01,76,[0 0 0 0; 1 -1 0 0],0);
D = sqrt(squeeze(sum(diff(XC(:,:,:),1,2).^2))); D = D([1 3],:);
disp(['mean simulation error: ' num2str(mean(fC))]);

% Gradient: Distance. Interpolate between 3 colors
DLin = linspace(sqrt(2)-.01,2+.01,nT);

