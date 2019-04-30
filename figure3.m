% Figure 3: Designing Folding Sequence
%% Prepare Space
clear; clc;

% Subplot Indices
nRow = [6 6 6];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
nCol = [5 5 5 5];   NCol = sum(nCol);   nC = [0 cumsum(nCol)];
figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nRow)
    for j = 1:length(nCol)
        cP = figM((nR(i)+1):(nR(i+1)-1), (nC(j)+1):(nC(j+1)-1));
        cellM(i+(j-1)*length(nRow)) = {cP(:)};
    end
end

fig = figure(3); clf;
s = sqrt(3);


%% a: Single Network: Amplify
subplot(NRow,NCol,cellM{1}); cla;
Xs1 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
dXs1 = [ 0.0  0.0 -0.0;...
         0.8  0.0 -0.0];
lSh = .2;
Xu1 = [[-lSh; -.5] [.25*sqrt(3); .25]+[-1;s]*lSh/2];

Xs1f = [Xs1(1,:); -Xs1(2,:)+.5];
dXs1f = [dXs1(1,:); -dXs1(2,:)];
Xu1f = [Xu1(1,:); -Xu1(2,:)+.5];

conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
construct_motion(Xs1, dXs1, Xu1, conn1, 1, 1);
construct_motion(Xs1f+[2;0], dXs1f, Xu1f+[2;0], conn1, 1, 1);
construct_motion(Xs1+[4;0], dXs1, Xu1+[4;0], conn1, 1, 1);
construct_motion(Xs1f+[6;0], dXs1f, Xu1f+[6;0], conn1, 1, 1);
axis([-1 8 -1.5 2]);

   
%% b: Combined Network
subplot(NRow,NCol,cellM{2}); cla;
% Left ==> Right
Xs1c = [Xs1 Xs1f+[s/2;0]];
Xu1c = [Xu1 Xu1f+[s/2;0]];
conn1c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
[Xs1a,Xu1a,conn1a] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[4;1]);
dXs1G = zeros(size(Xs1a)); dXs1G(2,1) =  .8;
construct_motion(Xs1a,dXs1G,Xu1a,conn1a,1,1);
axis([-1 8 -1.5 2]);


%% c: Propagate Combined Motion
subplot(NRow,NCol,cellM{3}); cla;
[Xs1a,Xu1a,conn1a] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[8;1]);
[XC,fC] = sim_motion(Xs1a,Xu1a,conn1a,.05,226,-[Xs1a Xu1a],0);
visualize_network(XC(:,1:size(Xs1a,2),end),...
                  XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),end),conn1a);
axis([-1 8 -1.5 2]+[3.5 3.5 0 0]);


%% d: Deployable Slope
subplot(NRow,NCol,cellM{4}); cla;
Xs2 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
dXs2 = [ 0.0  0.0 -0.0;...
         0.0 -0.8  0.0];
lSh = .6;
Xu2 = [[- lSh; 1] [ lSh; 1]];

Xs2f = [Xs2(1,:); -Xs2(2,:)+.5];
dXs2f = [dXs2(1,:); -dXs2(2,:)];
Xu2f = [Xu2(1,:); -Xu2(2,:)+.5];

conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
construct_motion(Xs2, dXs2, Xu2, conn2, 1, 1);
construct_motion(Xs2f+[2;0], dXs2f, Xu2f+[2;0], conn2, 1, 1);
construct_motion(Xs2+[4;0], dXs2, Xu2+[4;0], conn2, 1, 1);
construct_motion(Xs2f+[6;0], dXs2f, Xu2f+[6;0], conn2, 1, 1);
axis([-1 8 -1.5 2]);


%% e: Combined Network
subplot(NRow,NCol,cellM{5}); cla;
% Left ==> Right
Xs2c = [Xs2 Xs2f+[s/2;0]];
Xu2c = [Xu2 Xu2f+[s/2;0]];
conn2c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
[Xs2a,Xu2a,conn2a] = tesselate_network_old(Xs2c,Xu2c,conn2c,[s;0],[4;1]);
dXs2G = zeros(size(Xs2a)); 
dXs2G(2,1:2:size(Xs2a,2)) = .5;
dXs2G(2,2:2:size(Xs2a,2)) = -.5;
construct_motion(Xs2a,dXs2G,Xu2a,conn2a,1,1);
axis([-1 8 -1.5 2]);


%% f: Propagate Combined Motion
subplot(NRow,NCol,cellM{6}); cla;
[Xs2a,Xu2a,conn2a] = tesselate_network_old(Xs2c,Xu2c,conn2c,[s;0],[8;1]);
[XC,fC] = sim_motion(Xs2a,Xu2a,conn2a,.05,60,[Xs2a Xu2a],0);
visualize_network(XC(:,1:size(Xs2a,2),end),...
                  XC(:,[1:size(Xu2a,2)]+size(Xs2a,2),end),conn2a);
axis([-1 8 -1.5 2]+[3 3 0 0]);





% %% Animation
% fig = figure(4); clf;
% fName = 'animation.gif';
% dT = 0.03;
% % dT = 0.1;
% nXSh = 2.4;
% nYSh = 2;
% nSV = 1;
% XdotV = XC(:,1:size(Xs1a,2),1); XdotV = XdotV-mean(XdotV,2);
% 
% for i = 1:40:size(XC,3)
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
%     visualize_network(XC(:,1:size(Xs1a,2),i)/20 + [nXSh;nYSh],...
%                       XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i)/20 + [nXSh;nYSh],conn1a,.5);
%     axis([1.4 3.53 1.4 3.5]);
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
% 
% 
%% Animate Network Alone
fig = figure(5); clf;
fName = 'animation_net.gif';
dT = 0.03;
% dT = 0.1;
nXSh = 1.6;
nYSh = 2;
nSV = 1;
Xs1a = Xs2a;
Xu1a = Xu2a;
conn1a = conn2a;
XdotV = XC(:,1:size(Xs1a,2),1); XdotV = XdotV-mean(XdotV,2);

for i = 1:5:size(XC,3)
    cla;
%     dP = [D1(:,i)';D1(:,i)']; dP = dP(:);
%     dPa = dP(1:end-1); dPb = [1;dP(3:end)];
%     [Us, Uu, ~] = construct_motion(XC(:,1:size(Xs1a,2),i),...
%                                    XdotV,...
%                                    XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),...
%                                    conn1a, 0, 0);
%     Us = 0.3*Us/sqrt(sum(sum(Us.^2))); 
%     Uu = 0.3*Uu/sqrt(sum(sum(Uu.^2)));
%     XdotV = Us;
%     quiver(XC(1,:,i)/10 + nXSh, XC(2,:,i)/10 + nYSh,[Us(1,:) Uu(1,:)], [Us(2,:) Uu(2,:)],0);
    visualize_network(XC(:,1:size(Xs1a,2),i),...
                      XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),conn1a,1);
    axis([min(min(XC(1,:,:))) max(max(XC(1,:,:))) min(min(XC(2,:,:))) max(max(XC(2,:,:)))]);
%     construct_motion(XC(:,1:size(Xsa,2),i), XC(:,1:size(Xsa,2),i+1)-XC(:,1:size(Xsa,2),i), XC(:,[1:size(Xua,2)]+size(Xsa,2),i), conn, 20, 20);
%     axis([min(min(XC(1,:,:))) max(max(XC(1,:,:))) min(min(XC(2,:,:))) max(max(XC(2,:,:)))]);
    set(gca,'visible',0);
    drawnow;

    % Capture the plot as an image
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if nSV == 1
      imwrite(imind,cm,fName,'gif', 'Loopcount',inf,'DelayTime',dT);
      nSV = 0;
    else
      imwrite(imind,cm,fName,'gif','WriteMode','append','DelayTime',dT);
    end
end


