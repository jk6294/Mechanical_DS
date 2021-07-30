%% Node positions
sq = sqrt(3)/2;
Ls = sqrt(3);
Le = 1.7*Ls;

% Module
th = 30;
Xs = [-Ls/2  0  Ls/2;...
      -Ls*sq 0 -Ls*sq];
Xf = Le*[[-sind(th);-cosd(th)] [0;0] [sind(th);-cosd(th)]];
Xu0 = [-Ls/2  sqrt(Le^2/3-Ls^2/4)-Ls/sqrt(3);...
       -Ls/2 -sqrt(Le^2/3-Ls^2/4)-Ls/sqrt(3)]';
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Design positions
Xu = construct_network(Xs,Xf,Xu0,conn,0,1);
Xu = Xu(1:2,:);

% Simulate network
[XC,fC] = sim_motion10(Xs,Xu,conn,.01844,100,[Xs Xu],0);
D = sqrt(squeeze(sum(diff(XC(:,:,:),1,2).^2))); D = D([1 2],:);
Dc = sqrt(squeeze(sum(diff(XC(:,[1 3],:),1,2).^2)));

% Correct trajectory and rotation
XC = XC - XC(:,2,:);
for i = 1:size(XC,3)
    thv = atan2d(diff(XC(2,[1 3],i)),diff(XC(1,[1 3],i)));
    Rz = rotz(thv); Rz = Rz(1:2,1:2);
    XC(:,:,i) = Rz'*XC(:,:,i);
end

L = sqrt(squeeze(sum((XC(:,conn(:,1),1) - XC(:,conn(:,2),1)).^2)));
DLin = linspace(Ls-.01,Le+.01,nT);
DLin2 = linspace(1.52,3.97,nT);