%% Node positions
sq = sqrt(3)/2;
L}s = sqrt(3);
Le = 1.7*Ls;

% Unit 1
th = 30;
Xs = [-Ls/2  0  Ls/2;...
      -Ls*sq 0 -Ls*sq];
Xf = Le*[[-sind(th);-cosd(th)] [0;0] [sind(th);-cosd(th)]];
Xu0 = [-Ls/2  sqrt(Le^2/3-Ls^2/4)-Ls/sqrt(3);...
       -Ls/2 -sqrt(Le^2/3-Ls^2/4)-Ls/sqrt(3)]';
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Unit 2
th2 = 15;
Xs2 = [-Ls/2  0  Ls/2;...
      -Ls*sq 0 -Ls*sq];
Xf2 = Le*[[-sind(th2);-cosd(th2)] [0;0] [sind(th2);-cosd(th2)]];
Xu02 = [ 0.866  0.760;...
         0.000 -2.100]';

% Design positions
Xu = construct_network(Xs,Xf,Xu0,conn,0,1);
Xu = Xu(1:2,:);
Xu2 = construct_network(Xs2,Xf2,Xu02,conn,0,1);
Xu2 = Xu2(1:2,:);

% Simulate network
[XC,fC] = sim_motion10(Xs,Xu,conn,.018265,100,[Xs Xu],0);
D = sqrt(squeeze(sum(diff(XC(:,:,:),1,2).^2))); D = D([1 2],:);
[XC2,fC2] = sim_motion10(Xs2,Xu2,conn,.01775,100,[Xs2 Xu2],0);
D2 = sqrt(squeeze(sum(diff(XC2(:,:,:),1,2).^2))); D2 = D2([1 2],:);

% figure(6); clf; plot(D(1,:)-D(2,:));
% D(1,end)-D(2,end);

% Correct trajectory and rotation
XC = XC - XC(:,2,:);
for i = 1:size(XC,3)
    thv = atan2d(diff(XC(2,[1 3],i)),diff(XC(1,[1 3],i)));
    Rz = rotz(thv); Rz = Rz(1:2,1:2);
    XC(:,:,i) = Rz'*XC(:,:,i);
end
XC2 = XC2 - XC2(:,2,:);
for i = 1:size(XC2,3)
    thv = atan2d(diff(XC2(2,[1 3],i)),diff(XC2(1,[1 3],i)));
    Rz = rotz(thv); Rz = Rz(1:2,1:2);
    XC2(:,:,i) = Rz'*XC2(:,:,i);
end


L = sqrt(squeeze(sum((XC(:,conn(:,1),1) - XC(:,conn(:,2),1)).^2)));
L2 = sqrt(squeeze(sum((XC2(:,conn(:,1),1) - XC2(:,conn(:,2),1)).^2)));
DLin = linspace(Ls-.01,Le+.01,nT);
DLin2 = linspace(1.52,3.97,nT);