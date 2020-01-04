% Phase: Phase diagram
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% FP ==> FP
s = sqrt(3);
% Generate map
figure(1); clf;
Xs = [-s/2 s/2 0;...
      -1/2 -1/2 1];
Xf = 1.2*Xs;

conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_conic_finite(Xs,Xf,[-2 2; -2 2],[200;200],0,1,1);

% Get conic coordinates
[Q,W,v0,err] = construct_conic(Xs,Xf,1);
d = 2;
P = [W([1:d],:) v0(1:d);...
     zeros(1,d) 1]^-1;
Q = P'*Q*P;
CB1 = sample_conic(Q,[-2 2 -2 2],720,96);

% Plot
figure(1);
hold on;
scatter(CB1(1,:),CB1(2,:),[],parula(size(CB1,2)));
hold off;
drawnow;


% Phase diagram
% Evalulate slope at initial position
Cm1 = size(CB1,2);
% Distances
d1 = sqrt((Xs(:,1)-Xs(:,3))' * (Xs(:,1)-Xs(:,3)));
d2 = sqrt((Xs(:,2)-Xs(:,3))' * (Xs(:,2)-Xs(:,3)));
% Change in distances
d1dot1 = zeros(Cm1);
d2dot1 = zeros(Cm1);

UssM1 = zeros([size(Xs),Cm1,Cm1]);
L1 = size(CB1,2);
tic
fprintf([repmat('.',1,L1) '\n\n']);
parfor i = 1:L1
    fprintf('\b=\n');
    for j = 1:L1
        if i ~= j
            [Uss,Uus,err] = construct_motion(Xs,Xs,CB1(:,[i,j]),conn,0,0);
            UssM1(:,:,i,j) = Uss;
            d1dot1(i,j) = (Xs(:,1)-Xs(:,2))' * (Uss(:,1)-Uss(:,2));
            d2dot1(i,j) = (Xs(:,2)-Xs(:,3))' * (Uss(:,2)-Uss(:,3));
        end
    end
end
toc
d1dot1 = d1dot1/d1;
d2dot1 = d2dot1/d2;


%% Plot
figure(2); clf;
subplot(1,2,1);
dRat1 = d2dot1./d1dot1;
colormap(parula(2^11));
imagesc(tanh(abs(dRat1)));
% caxis([-1,1]*5);
hold on;
contour(abs(dRat1),[1,1],'color','r','linewidth',1);
hold off;

subplot(1,2,2);
dRat1T = (abs(dRat1) > 1)*2;
imagesc(dRat1T);
hold on;
contour(abs(dRat1),[1,1],'color','r','linewidth',1);
hold off;


%% Limit cycle maps
figure(1); clf;

s = sqrt(3);
nxSh = s/4;
nySh = 1/2;
Xs = [-s/2 -nxSh  s/2;...
      -1/2  nySh   -1/2];
Xf = [-s/2  nxSh  s/2;...
      -1/2  nySh   -1/2];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_conic_finite(Xs,Xf,[-2 2; -2 2],[200;200],0,1,1);

% Get conic coordinates
[Q,W,v0,err] = construct_conic(Xs,Xf,1);
d = 2;
P = [W([1:d],:) v0(1:d);...
     zeros(1,d) 1]^-1;
Q1 = P'*Q*P;
% Final
P2 = [W([1:d]+d,:) v0([1:d]+d);...
     zeros(1,d) 1]^-1;
Q2 = P2'*Q*P2;
[CB1 sInd] = sample_conic(Q1,[-1 1 -1 1] - [0 0 .5 .5],1000,-90);
% Remap
S = [W v0] * P * [CB1; ones(1,size(CB1,2))];
CB2 = S([1:d]+d,:);
% CB2 = sample_conic(Q2,[-2 2 -2 2]*1,1000,-90);


% Plot
figure(1);
hold on;
scatter(CB1(1,:),CB1(2,:),[],parula(size(CB1,2)));
scatter(CB2(1,:),CB2(2,:),[],parula(size(CB2,2)));
hold off;
drawnow;


% Phase diagram
% Evalulate slope at initial position
Cm1 = sInd;
Cm2 = size(CB1,2)-sInd;
% Distances
d11 = sqrt((Xs(:,1)-Xs(:,2))' * (Xs(:,1)-Xs(:,2)));
d21 = sqrt((Xs(:,2)-Xs(:,3))' * (Xs(:,2)-Xs(:,3)));
d12 = sqrt((Xf(:,1)-Xf(:,2))' * (Xf(:,1)-Xf(:,2)));
d22 = sqrt((Xf(:,2)-Xf(:,3))' * (Xf(:,2)-Xf(:,3)));
% Change in distances
d1dot1 = zeros(Cm1,Cm2);
d2dot1 = zeros(Cm1,Cm2);
d1dot2 = zeros(Cm1,Cm2);
d2dot2 = zeros(Cm1,Cm2);

UssM1 = zeros([size(Xs),Cm1,Cm2]);
UssM2 = zeros([size(Xs),Cm1,Cm2]);
tic
fprintf([repmat('.',1,Cm1) '\n\n']);
parfor i = 1:Cm1
    fprintf('\b=\n');
    for j = 1:Cm2
        [Uss,Uus,err] = construct_motion(Xs,Xs,CB1(:,[i,sInd+j]),conn,0,0);
        UssM1(:,:,i,j) = Uss;
        d1dot1(i,j) = (Xs(:,1)-Xs(:,2))' * (Uss(:,1)-Uss(:,2));
        d2dot1(i,j) = (Xs(:,2)-Xs(:,3))' * (Uss(:,2)-Uss(:,3));
    end
end
fprintf([repmat('.',1,Cm1) '\n\n']);
parfor i = 1:Cm1
    fprintf('\b=\n');
    for j = 1:Cm2
        [Uss,Uus,err] = construct_motion(Xf,Xf,CB2(:,[i,sInd+j]),conn,0,0);
        UssM2(:,:,i,j) = Uss;
        d1dot2(i,j) = (Xf(:,1)-Xf(:,2))' * (Uss(:,1)-Uss(:,2));
        d2dot2(i,j) = (Xf(:,2)-Xf(:,3))' * (Uss(:,2)-Uss(:,3));
    end
end
toc
d1dot1 = d1dot1/d11;
d2dot1 = d2dot1/d21;
d1dot2 = d1dot2/d12;
d2dot2 = d2dot2/d22;


% Plot
dRat1 = d2dot1./d1dot1;
dRat2 = d2dot2./d1dot2;

%%
figure(2); clf;
subplot(1,2,1);
colormap(parula(2^11));
imagesc(tanh(abs(dRat1 .* dRat2)));
hold on;
contour(abs(dRat1.*dRat2),[1,1],'color','r','linewidth',1);
% contour(abs(1./dRat2),[1,1],'color','c','linewidth',1);
% contour(abs(dRat1 .* dRat2),[1,1],'color','c','linewidth',1);
hold off;
dRat12 = (dRat1 .* dRat2);
subplot(1,2,2);
imagesc(abs(dRat1.*dRat2)>1 & dRat1 > 0 & dRat2 < 0);


%% Create and simulate network
% niInd2 = 161;
% niInd1 = 202;
% niInd2 = 437;
% niInd1 = 156;
niInd2 = 244;
niInd1 = 300;
Xu = CB1(:,[niInd1 niInd2+sInd]);
[XCf,fCf] = sim_motion(Xs,Xu,conn,.005,2000,[Xs Xu],0);
[XCb,fCb] = sim_motion(Xs,Xu,conn,.005,2000,-[Xs Xu],0);
XC = cat(3,flip(XCf,3),XCb);

d1C = squeeze(sqrt(sum((XC(:,1,:) - XC(:,2,:)).^2)));
d2C = squeeze(sqrt(sum((XC(:,2,:) - XC(:,3,:)).^2)));

figure(3); clf;
visualize_network(Xs,Xu,conn);

figure(4); clf;
plot(d1C,d2C,'linewidth',1);
hold on;
plot(d2C,d1C,'linewidth',1);
plot(d11,d21,'x','linewidth',5);
plot(d12,d22,'x','linewidth',5);
line([min([d1C;d2C])  max([d1C;d2C])], [min([d1C;d2C])  max([d1C;d2C])],...
     'color', 'r');
hold off;
dSlo = diff(d2C)./diff(d1C);
disp(dSlo(abs(d1C-d2C)< .005));
disp(dRat12(niInd1,niInd2));


%% Tesselate
% Remap point
Xu2 = [W v0] * P * [Xu; ones(1,size(Xu,2))];
Xu2 = Xu2([1:d]+d,:);

s = sqrt(3);
Xsc = [-s/2 -nxSh  s/2 s-nxSh ;...
       -1/2  nySh -1/2   nySh];
Xuc = [Xu Xu2.*[1;-1]+[s/2-nxSh;-0.5 + nySh]];
connc = [1 5; 2 5; 3 5; 1 6; 2 6; 3 6; 2 7; 3 7; 4 7; 2 8; 3 8; 4 8];

[Xsa, Xua ,conna] = tesselate_network_old(Xsc,Xuc,connc,[s;0],[10,1]);
x0 = zeros(2,max(conna(:))); x0(1,2) = 1;
[XCa,fCa] = sim_motion(Xsa,Xua,conna,.002,3000,x0,0);

% Distances
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
disp(max(fCa));
disp('done');


%%
h = figure(6); clf;
fname = 'chaos_net5.gif';

set(gcf, 'Renderer', 'opengl'); 
% Combined Network Distances
D1 = sqrt(squeeze(sum(diff(XCa,1,2).^2)));
D1 = D1(1:size(Xsa,2),:);
plot(d1,d2,'k-','linewidth',1);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
% Cobweb 3
lwL = 1;
arrHL = 5;
arrHW = 5;
caF = .8;
cTr3 = [015 082 186]/255;
cSr  = [255 100 100]/255;
xMax = max(max(XCa(1,:,:))); xMin = min(min(XCa(1,:,:)));
yMax = max(max(XCa(2,:,:))); yMin = min(min(XCa(2,:,:)));
nSk = 40;
for i = 1:nSk:size(D1,2)
%     cla;
    subplot(1,4,1);
    plot(d1,d2,'k-','linewidth',1);
    hold on;
    plot([.2 1.3],[.2 1.3], '--', 'color', [200 200 200]/255);
    pInd3 = i;
    dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:); dP = dP(1:end-2);
    dPa = dP(1:end-1); dPb = [1;dP(3:end)];
    line(dP(1:end-1),[1;dP(3:end)],'color',cTr3,'linewidth',lwL);
%     for i = 1:length(dP)-2
%         ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr3);
%         set(ah,'parent',gca);
%         set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
%     end
    hold off;
    axis([.8 2.5 .8 2.5]);
    
    subplot(1,4,[2:4]); cla;
    visualize_network(XCa(:,1:max(conna(:,1)),i),...
                      XCa(:,max(conna(:,1))+1:max(conna(:,2)),i),conna,2);
    axis([xMin xMax yMin yMax]);
    drawnow;
    
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     if i == 1 
%         imwrite(imind,cm,fname,'gif', 'Loopcount',inf, 'DelayTime', 0.03); 
%     else 
%         imwrite(imind,cm,fname,'gif','WriteMode','append', 'DelayTime', 0.03); 
%     end 
end

% d2 = maps(Xs,Xu);
% f = matlabFunction(d2);

