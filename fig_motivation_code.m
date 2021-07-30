%% Panel b-c: quadrifolium
% conncc:                       quadrifolium connectivity
% Xscc:                         initial specified node positions
% Xucc:                         initial unspecified node positions
% XCcc:                         full network trajectory
% CSSc:                         unspecified node color
% fCcc:                         energy cost

% Load data
load quadrifolium.mat;
Ns = 114;               % Number of non-added nodes in quad example

% Correct rotation
XCcc = XCcc - XCcc(:,Ns,:);
XCcc(1,:,:) = -XCcc(1,:,:);
for i = 1:size(XCcc,3)
    Rz = rotz(atan2d(diff(XCcc(2,[-2 0]+Ns,i)),diff(XCcc(1,[-2 0]+Ns,i))));
    XCcc(:,:,i) = Rz(1:2,1:2)'*XCcc(:,:,i);
end


%% Panel d: mechanical AND gate
% Define Unit
s = sqrt(3);
Xs = [-s/2 0 s/2;...
      -1/2 1 -1/2];
Xu = [-s/2 -s/2; sqrt(s^2-s^2/4) -sqrt(s^2-s^2/4)];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Define branches
nR2 = 12;
[Xsa2,Xua2,conna2] = network_chain_x(Xs(:,1:2), Xu, ones(1,nR2));
[Xsa3,Xua3,conna3] = network_chain_x(Xs(:,1:2), Xu, ones(1,nR2-1));
% Center
Xua2 = Xua2 - Xsa2(:,end);
Xsa2 = Xsa2 - Xsa2(:,end);
Xua3 = Xua3 - Xsa3(:,end);
Xsa3 = Xsa3 - Xsa3(:,end);
% % Flip
Xua3(2,:) = -Xua3(2,:);
Xsa3(2,:) = -Xsa3(2,:);
% Rotate
Rz = rotz(-60);
Xua3 = Rz(1:2,1:2)*Xua3;
Xsa3 = Rz(1:2,1:2)*Xsa3;

% Define network
[Xsa4,Xua4,conna4] = tesselate_network_old([Xsa2 Xsa3], [Xua2 Xua3],...
                     [conna2(:,1) conna2(:,2)+max(conna3(:,1));...
                      conna3(:,1)+max(conna2(:,1)) conna3(:,2)+max(conna2(:,2))], [0;0], [1;1]);

% Simulate
X041 = zeros(2,size(Xsa4,2)+size(Xua4,2));
X042 = X041;
X041(1,1) = -1;
X042(:,11) = [-1/2;sqrt(3)/2];
[XCa4a,fCa4a] = sim_motion10(Xsa4,Xua4,conna4,.1,116,X041,0);
[XCa4b,fCa4b] = sim_motion10(Xsa4,Xua4,conna4,.1,117,X042,0);
[XCa4c,fCa4c] = sim_motion10(XCa4a(:,1:size(Xsa4,2),end),...
                             XCa4a(:,(1:size(Xua4,2))+size(Xsa4,2),end),...
                             conna4,.125,319,X042,0);
disp(['mean simulation error: ' [num2str(mean(fCa4a)) '  ' num2str(mean(fCa4c))]]);
XCa4 = cat(3,XCa4a, XCa4c(:,:,2:end));

% Correct position and rotation
% Initial
XP0 = XCa4b(:,:,1); XP0 = XP0 - XP0(:,size(Xsa4,2));
% Top extension
XP1 = XCa4b(:,:,end); XP1 = XP1 - XP1(:,size(Xsa4,2));
Rz = rotz(atan2d(diff(XP1(2,[1 3])),diff(XP1(1,[1 3]))));
XP1 = Rz(1:2,1:2)'*XP1;
% Left extension
XP2 = XCa4a(:,:,end); XP2 = XP2 - XP2(:,size(Xsa4,2));
Rz = rotz(atan2d(diff(XP2(2,[1 3])),diff(XP2(1,[1 3]))));
XP2 = Rz(1:2,1:2)'*XP2;
% Full extension
XP3 = XCa4(:,:,end); XP3 = XP3 - XP3(:,size(Xsa4,2));
Rz = rotz(atan2d(diff(XP3(2,[1 3])),diff(XP3(1,[1 3]))));
XP3 = Rz(1:2,1:2)'*XP3;

% Distances
DLin = linspace(sqrt(3)-.01,3+.01,nT);


