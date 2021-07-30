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

%% Node positions
sq = sqrt(3)/2;
Ls = sqrt(3);
Le = 1.7*Ls;

% Module
th = 20;
Xs = [-Ls/2 0 Ls/2;...
      -Ls*sq 0 -Ls*sq];
Xf = Le*[[-sind(th);-cosd(th)] [0;0] [sind(th);-cosd(th)]];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
Xu = construct_network(Xs,Xf,[-Ls/2 .566; -Ls/2 -2.203]',conn,0,1);
Xu = Xu(1:2,:);

% Tesselate network
Xsc = [-Ls/2 0 Ls/2 Ls;...
       -Ls*sq 0 -Ls*sq 0];
Xuc = [Xu [Xu(1,:)+Ls/2;-Xu(2,:) - Ls*sq]];
connc = [1 5; 2 5; 3 5; 1 6; 2 6; 3 6; 2 7; 3 7; 4 7; 2 8; 3 8; 4 8];
nRep = 8; ns = nRep+2;
[Xsa,Xua,conna,~] = network_chain_x(Xs(:,1:2),Xu,ones(1,nRep));
conna = sortrows(conna,2);

% Simulate networks
[XC,fC] = sim_motion10(Xs,Xu,conn,.019,100,[Xs Xu],0);
[XCa,fCa] = sim_motion10([Xsa Xua],[],conna,.091327,100,[Xsa Xua],0);

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
 
L = sqrt(squeeze(sum((XC(:,conn(:,1),1) - XC(:,conn(:,2),1)).^2)));
DLin = linspace(Ls-.01,Le+.01,nT);
DLin2 = linspace(1.52,3.97,nT);
D = sqrt(squeeze(sum(diff(XC(:,:,:),1,2).^2))); D = D([1 2],:);
Da = sqrt(squeeze(sum(diff(XCa(:,[1:max(conna(:,1))],:),1,2).^2))); Da = Da(1:nRep+1,:);

