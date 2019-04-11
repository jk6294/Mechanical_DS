%% Prepare Space
clear; clc;


%% Module
s = sqrt(3);
Xs = [ 0.0  0.5  1.0  1.5  2.0  2.5;...
       0.0 -s/2  0.0  s/2  0.0 -s/2];
conn = [1 3; 1 4; 2 3; 2 4; 3 5; 3 6; 4 5; 4 6];

figure(1); clf;
subplot(2,8,1);
visualize_network(Xs,[],conn);

subplot(2,8,[2:8]);
[Xsa, conna] = tesselate_network(Xs,conn,[2;0], [20,1]);
visualize_network(Xsa,[],conna);
axis([0 max(max(Xsa)) -1.2 1.2]);

subplot(2,8,[2:8]+8);
[XC, fC] = sim_motion(Xsa,[],conna,.03,3000,-Xsa,0);
visualize_network(XC(:,:,end),[],conna);
axis([0 max(max(Xsa)) -1.2 1.2]);


%% Rigidity
% Second Conformation
Xsa2 = XC(:,:,end);

figure(2); clf;
nVa = size(Xsa,2);

% Initial
LVal = sqrt(squeeze(sum((Xsa(:,conna(:,1))-Xsa(:,conna(:,2))).^2)))';
R = (1./LVal).*rigidity(Xsa,conna);
Q = R'*R;
[V, U] = eig(Q);
u = diag(U);
subplot(2,3,1);
imagesc(V);

% Final
R2 = (1./LVal).*rigidity(Xsa2,conna);
Q2 = R2'*R2;
DM = [zeros(size(Q2)) eye(size(Q2));...
      -Q2             zeros(size(Q2))];

[V2, U2] = eig(Q2);
u2 = diag(U2);
subplot(2,3,4);
imagesc(V2);

% Rigidity Propagation
% r = zeros(size(XC,3),1);
% for i = 1:size(XC,3)
%     RP = rigidity(XC(:,:,i),conna);
%     r(i) = det(RP*RP');
% end

% Eigenvectors of minimum y displacement
subplot(2,3,2);
v = sum((V([1:2:nVa]+nVa,:)).^2);
vInd = find(v==min(v));
imagesc(V([1:2:nVa]+nVa,:));

subplot(2,3,3);
plot(V([1:2:nVa],vInd));

subplot(2,3,5);
% plot(r);

subplot(2,3,6);
pInd = 45;
Xsdot = zeros(size(Xsa)); Xsdot(1,2)=.4;
% Xsdot(:,[1:20]+41) = [V([1:20],pInd) V([1:20]+82,pInd)]';
Xsdot = [V([1:nVa],pInd)'; V([1:nVa]+nVa,pInd)']*.1;
M = ones(size(Xsa,2),1);
K = ones(size(conna,1),1);
[XCd, fCd, t, KE, PE] = sim_dynamics(Xsa',Xsdot',conna,M,K*10,100,5000);

figure(7); clf;
XCv = diff(XCd,1,3);
XCvA = reshape(XCv,[2*nVa,size(XCd,3)-1]);
XCVel = squeeze(XCv(:,1,:).^2 + XCv(:,2,:).^2);
imagesc(XCVel);


%% Animate
XCP = permute(XCd, [2 1 3]);
fig = figure(9); clf;
fName = 'animation.gif';
nSV = 1;
dT = 0.03;
nV = [1:3];
nS = 2;
for i = 1:50:size(XCP,3)
    cla;
    visualize_network(XCP(:,:,i), [], conna);
    axis([min(min(min(XCP(1,:)))) max(max(max(XCP(1,:)))) min(min(min(XCP(2,:))))  max(max(max(XCP(2,:))))]);
    drawnow;

%     % Capture the plot as an image
%     frame = getframe(fig);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if nSV == 1
%       imwrite(imind,cm,fName,'gif', 'Loopcount',inf,'DelayTime',dT);
%       nSV = 0;
%     else
%       imwrite(imind,cm,fName,'gif','WriteMode','append','DelayTime',dT);
%     end
end