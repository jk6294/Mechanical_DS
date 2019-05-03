% Figure 3: Designing Folding Sequence
%% Prepare Space
clear; clc;

% Subplot Indices
NRow = 20;
NCol = 14;
nBal = {[1 1 1], [1 1 1], [1 1 1], [1 1 1 1]};
cellL = [0 cumsum(cellfun('length',nBal))];
cellM = cell(sum(cellL));
for i = 1:length(nBal)
    s = nBal{i}; 
    su = cumsum(s/sum(s)*NCol);
    sd = 
    for j = 1:length(s)
        cellM(j+cellL(i)) = 
    end
end

figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nRow)-1
    for j = 1:length(nCol)
        cP = figM((nR(i)+1):(nR(i+1)-1), (nC(j)+1):(nC(j+1)-1));
        cellM(i+(j-1)*(length(nRow)-1)) = {cP(:)};
    end
end

% Figure Axis Bounds
axM = [0 9 0 3];

fig = figure(3); clf;
s = sqrt(3);
pSc = 0.6;


%% a: Slope 1 Module
subplot(NRow,NCol,cellM{1}); cla;

% Single Module
Xs1 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
dXs1 = [ 0.0  0.0 -0.0;...
         0.5 -0.5  0.5];
lSh = .6;
Xu1 = [[- lSh; 1] [ lSh; 1]];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Flipped Module
Xs1f = [Xs1(1,:); -Xs1(2,:)+.5];
dXs1f = [dXs1(1,:); -dXs1(2,:)];
Xu1f = [Xu1(1,:); -Xu1(2,:)+.5];

% Combined Module
Xs1c = [Xs1 Xs1f+[s/2;0]];
Xu1c = [Xu1 Xu1f+[s/2;0]];
conn1c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
[Xs1a,Xu1a,conn1a] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[3;1]);
dXs1G = zeros(size(Xs1a)); 
dXs1G(2,1:2:size(Xs1a,2)) = .5;
dXs1G(2,2:2:size(Xs1a,2)) = -.5;

% Show Module
construct_motion(Xs1, dXs1, Xu1, conn1, 1, 1, pSc);
construct_motion(Xs1a+[2.5;0],dXs1G,Xu1a+[2.5;0],conn1a,1,1, pSc);

axis(axM - [1 1 1.3 1.3]);


%% b: Slope 1 Propagate Combined Motion
subplot(NRow,NCol,cellM{2}); cla;
[Xs1a,Xu1a,conn1a] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[8;1]);
[XC,fC] = sim_motion(Xs1a,Xu1a,conn1a,.05,60,[Xs1a Xu1a],0);
visualize_network(XC(:,1:size(Xs1a,2),end),...
                  XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),end),conn1a, pSc);
axis(axM+[2.5 2.5 -1.3 -1.3]);


%% c: Slope -2 Module
subplot(NRow,NCol,cellM{3}); cla;

% Single Module
Xs2 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
dXs2 = [ 0.0  0.0  0.0;...
         -0.2  0.0  -0.3]*2;
Xu2 = [ 0.30  0.8;...
       -0.50  1.3];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Flipped Module
Xs2f = [Xs2(1,:); -Xs2(2,:)+.5];
dXs2f = [dXs2(1,:); -dXs2(2,:)];
Xu2f = [Xu2(1,:); -Xu2(2,:)+.5];
[Uss,Uus,err] = construct_motion(Xs2f+[6;0], dXs2f, Xu2f+[6;0], conn2,0,0);

% Combined Module
Xs2c = [Xs2 Xs2f+[s/2;0]];
Xu2c = [Xu2 Xu2f+[s/2;0]];
conn2c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
[Xs2a,Xu2a,conn2a] = tesselate_network_old(Xs2c,Xu2c,conn2c,[s;0],[3;1]);
dXs2G = zeros(size(Xs2a)); 
dXs2G(2,end) = 1;
dXs2G(2,end-6) = dXs2G(2,end)*(Uss(2,1)/Uss(2,3))^3;

% Show Module
construct_motion(Xs2, dXs2, Xu2, conn2, 1, 1, pSc);
construct_motion(Xs2a+[2.5;0],dXs2G,Xu2a+[2.5;0],conn2a,1,1, pSc);
axis(axM - [1 1 1.3 1.3]);


%% d: Slope -2 Propagate Combined Motion
subplot(NRow,NCol,cellM{4}); cla;
[Xs2a,Xu2a,conn2a] = tesselate_network_old(Xs2c,Xu2c,conn2c,[s;0],[8;1]);
[XC,fC] = sim_motion(Xs2a,Xu2a,conn2a,.05, 18, [Xs2a Xu2a],0);
visualize_network(XC(:,1:size(Xs2a,2),end),...
                  XC(:,[1:size(Xu2a,2)]+size(Xs2a,2),end),conn2a, pSc);
axis(axM+[2.5 2.5 -1.3 -1.3]);


%% e: Slope 0 Module
subplot(NRow,NCol,cellM{5}); cla;

% Single Module
Xs3 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
dXs3 = [ 0.0  0.0 -0.0;...
         0.5  0.0 -0.0];
lSh = 0.3;
Xu3 = [[-lSh; -.5] [.25*s; .25]+[-1;s]*lSh/2];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Single Module Flipped
Xs3f = [Xs3(1,:); -Xs3(2,:)+.5];
dXs3f = [dXs3(1,:); -dXs3(2,:)];
Xu3f = [Xu3(1,:); -Xu3(2,:)+.5];

% Combined Module
Xs3c = [Xs3 Xs3f+[s/2;0]];
Xu3c = [Xu3 Xu3f+[s/2;0]];
conn3c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3c,Xu3c,conn3c,[s;0],[3;1]);
dXs3G = zeros(size(Xs3a)); dXs3G(:,1) =  dXs3(:,1);

% Visualize
construct_motion(Xs3/1.2, dXs3, Xu3/1.2, conn3, 1, 1, pSc);
construct_motion(Xs3a+[2.5;0],dXs3G,Xu3a+[2.5;0],conn3a,1,1, pSc);
axis(axM - [1 1 1.3 1.3]);


%% f: Slope 0 Propagate Combined Motion
subplot(NRow,NCol,cellM{6}); cla;
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3c,Xu3c,conn3c,[s;0],[8;1]);
[XC,fC] = sim_motion(Xs3a,Xu3a,conn3a,.05,170,-[Xs3a Xu3a],0);
visualize_network(XC(:,1:size(Xs3a,2),end),...
                  XC(:,[1:size(Xu3a,2)]+size(Xs3a,2),end),conn3a, pSc);
axis(axM+[2 2 -1.3 -1.3]);


%% g: Place Modules
subplot(NRow,NCol,cellM{7}); cla;

% Construct Unit
Xs3 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
Xs3f = [Xs3(1,:); -Xs3(2,:)+.5];
lSh = .3;
Xu3 = [[-lSh; -.5] [.25*s; .25]+[-1;s]*lSh/2];

% lSh = -0.3;
% lSh2 = 0.000;
% Xu3 = [Xs3(:,1)-lSh*[s;-3]+lSh2*[3;s] Xs3(:,2)+lSh*[s;-3]+lSh2*[3;s]];

Xu3 = [-Xu3(1,:); Xu3(2,:)];
Xu3f = [Xu3(1,:); -Xu3(2,:)+.5];
Xs3c = [Xs3 Xs3f+[s/2;0]];
Xu3c = [Xu3 Xu3f+[s/2;0]];
conn1c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
nR = 5;
[Xs3a1,Xu3a1,conn3a1] = tesselate_network_old(Xs3c,Xu3c,conn1c,[s;0],[nR;1]);
% Remove Offset
Xu3a1 = Xu3a1-Xs3a1(:,1);
Xs3a1 = Xs3a1-Xs3a1(:,1);

% Modified Unit: Extend Length 1
Xs3a2 = Xs3a1+[nR*s;0];
Xu3a2 = Xu3a1+[nR*s;0];

% Modified Unit: Branch 1
Xs3a3 = [Xs3a1(1,:); -Xs3a1(2,:)]; 
Xu3a3 = [Xu3a1(1,:); -Xu3a1(2,:)]; 
Rz = rotz(60); Rz = Rz(1:2,1:2);
Xs3a3 = Rz*Xs3a3 + [(nR-.5)*s;1.5];
Xu3a3 = Rz*Xu3a3 + [(nR-.5)*s;1.5];

% Modified Unit: Branch 2
Xs3a4 = Xs3a1+[1.5*nR*s;1.5*nR];
Xu3a4 = Xu3a1+[1.5*nR*s;1.5*nR];

% Visualize
visualize_network(Xs3a1,Xu3a1,conn3a1, pSc/2);
visualize_network(Xs3a2+[1;0],Xu3a2+[1;0],conn3a1, pSc/2);
visualize_network(Xs3a3+[0;1],Xu3a3+[0;1],conn3a1, pSc/2);
axis([0 9 -.3 3.2]*3.2);


%% Combine and Simulate Modules
Xs3a = [Xs3a1 Xs3a2 Xs3a3];
Xu3a = [Xu3a1 Xu3a2 Xu3a3];
conn3a = conn3a1;
conn3a = [[conn3a(:,1), conn3a(:,2)+max(conn3a1(:,1))];...
          [conn3a1(:,1)+max(conn3a(:,1)), conn3a1(:,2)+max(conn3a(:,2))]];
conn3a = [[conn3a(:,1), conn3a(:,2)+max(conn3a1(:,1))];...
          [conn3a1(:,1)+max(conn3a(:,1)), conn3a1(:,2)+max(conn3a(:,2))]];
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3a,Xu3a,conn3a,[1;1],[1;1]);

XsDot = [zeros(size(Xs3a)) zeros(size(Xu3a))];
XsDot(2,size(Xs3a,2)) = -1;

[XC,fC] = sim_motion(Xs3a,Xu3a,conn3a,.01,3200,XsDot,0);

%% h: Visualize
subplot(NRow,NCol,cellM{6}); cla;
pInd = 600;
visualize_network(XC(:,1:size(Xs3a,2),pInd),...
                  XC(:,[1:size(Xu3a,2)]+size(Xs3a,2),pInd),conn3a, pSc/2);
axis([0 9 -.3 3.2]*3.2);


%% i: More Collapsed
subplot(NRow,NCol,cellM{9}); cla;
pInd = 3190;
XCD = -(XC(:,:,pInd)-XC(:,:,pInd-1));
XCD = 2*XCD / sqrt(sum(diag(XCD*XCD')));
R = rigidity(XC(:,:,pInd),conn3a);
[ua,sa,va] = svds(R,2,'smallest');
va = -2*va/sqrt(va'*va);
% hold on;
% quiver(XC(1,:,pInd),XC(2,:,pInd),va(1:size(XC,2),2)',...
%        va([1:size(XC,2)]+size(XC,2),2)',0,'linewidth',1.2,'color',[0,.8,.8]);
% quiver(XC(1,:,pInd),XC(2,:,pInd),va(1:size(XC,2),2)',...
%        va([1:size(XC,2)]+size(XC,2),2)',0,'linewidth',.3,'color',[1 1 1]);
% hold off;
% construct_motion(XC(:,1:size(Xs3a,2),pInd),XCD(:,1:size(Xs3a,2)),...
%                  XC(:,[1:size(Xu3a,2)]+size(Xs3a,2),pInd),conn3a,1,1, pSc/2);
visualize_network(XC(:,1:size(Xs3a,2),pInd),XC(:,[1:size(Xu3a,2)]+size(Xs3a,2),pInd),conn3a,pSc/2);
axis([0 9 -.3 3.2]*2.5);


%% Size and Save Figure
fName = 'figure3';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters'; 
fig.PaperPosition = [-2.1 -0.3 18.4 5.9];
fig.PaperSize = [14 5.5];
saveas(fig, ['Figures/' fName], 'pdf');


%% Animation
% fig = figure(4); clf;
% fName = 'animation.gif';
% dT = 0.03;
% % dT = 0.1;
% nXSh = 1.4;
% nYSh = 2;
% nSV = 1;
% Xs1a = Xs2a;
% Xu1a = Xu2a;
% conn1a = conn2a;
% XdotV = XC(:,1:size(Xs1a,2),1); XdotV = XdotV-mean(XdotV,2);
% 
% D1 = sqrt(squeeze(sum(diff(XC(:,1:size(Xs1a,2),:),1,2).^2)));
% D1 = D1(1:2:end,:); 
% 
% for i = 1:20:size(XC,3)
%     cla;
%     dP = [D1(:,i)';D1(:,i)']; dP = dP(:);
%     dPa = dP(1:end-1); dPb = [1;dP(3:end)];
%     plot(d1,d2,'k-');
%     hold on;
% %     [Us, Uu, ~] = construct_motion(XC(:,1:size(Xs1a,2),i),...
% %                                    XdotV,...
% %                                    XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),...
% %                                    conn1a, 0, 0);
% %     Us = 0.3*Us/sqrt(sum(sum(Us.^2))); 
% %     Uu = 0.3*Uu/sqrt(sum(sum(Uu.^2)));
% %     XdotV = Us;
% %     quiver(XC(1,:,i)/10 + nXSh, XC(2,:,i)/10 + nYSh,[Us(1,:) Uu(1,:)], [Us(2,:) Uu(2,:)],0);
%     plot([0 4],[0 4],'LineStyle','--','color',[0 0 0 .5]);
%     line(dP(1:end-1),[1;dP(3:end)]);
%     visualize_network(XC(:,1:size(Xs1a,2),i)/10 + [nXSh;nYSh],...
%                       XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i)/10 + [nXSh;nYSh],conn1a,.5);
%     axis([0.6 2.7 0.6 2.7]);
% %     construct_motion(XC(:,1:size(Xsa,2),i), XC(:,1:size(Xsa,2),i+1)-XC(:,1:size(Xsa,2),i), XC(:,[1:size(Xua,2)]+size(Xsa,2),i), conn, 20, 20);
% %     axis([min(min(XC(1,:,:))) max(max(XC(1,:,:))) min(min(XC(2,:,:))) max(max(XC(2,:,:)))]);
%     set(gca,'visible',0);
%     drawnow;
% 
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
% end


%% Animate Network Alone
% fig = figure(5); clf;
% fName = 'animation_net.gif';
% dT = 0.03;
% % dT = 0.1;
% nXSh = 1.6;
% nYSh = 2;
% nSV = 1;
% Xs1a = Xs3a;
% Xu1a = Xu3a;
% conn1a = conn3a;
% XdotV = XC(:,1:size(Xs1a,2),1); XdotV = XdotV-mean(XdotV,2);
% 
% for i = 1:20:size(XC,3)
%     cla;
% %     R = rigidity(XC(:,:,i),conn1a);
% %     [u,sa,va] = svds(R,1,'smallest');
%     visualize_network(XC(:,1:size(Xs1a,2),i),...
%                       XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),conn1a,1);
%     hold on;
% %     quiver(XC(1,:,i),XC(2,:,i),-va(1:size(XC,2))',-va([1:size(XC,2)]+size(XC,2))',0);
%     hold off;
%     axis([min(min(XC(1,:,:))) max(max(XC(1,:,:))) min(min(XC(2,:,:))) max(max(XC(2,:,:)))]);
%     set(gca,'visible',0);
%     drawnow;
% 
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
% end


