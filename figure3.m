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


% a: Single Network: Amplify
subplot(NRow,NCol,cellM{1});
Xs10 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
Xs1T = [-s/2  0    s/2;...
        -1.0  1.5 -1.0];
dXs1 = [ 0.0  0.0 -0.0;...
        -0.4  0.0 -0.0];
visualize_conic(Xs10,dXs1,[-2 2;-2 2],[400;400],0,1,1);
visualize_conic_finite(Xs10,Xs1T,[-2 2;-2 2],[400;400],5,1,1);


%% b: Simulate Single
subplot(NRow,NCol,cellM{2});
Xu1a = [0.3; 1.5];
Xu1b = [1.3; 1.5];
Xu1c = [-1.5; -1];
Xu1d = [1.0;-1];

% Xu1e = [-0.43; 1.75];
% Xu1f = [1.308; -1.25];
Xu1e = [-0.25; 1.50];
Xu1f = [1.5; -1.00];
Xu1g = [-0.29; 1.5];
Xu1h = [1.16; -1.0];

Xu1 = [Xu1g Xu1h];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
LVal = sqrt(sum((Xs10(:,conn1(:,1))-Xu1(:,conn1(:,2)-3)).^2));
[XC,fC] = sim_motion(Xs10,Xu1,conn1,.01,2000,-[Xs10 Xu1],0);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
plot(d1,d2);
visualize_network(Xs10,Xu1,conn1);
hold on;
plot([0 5],[0 5],'--');
hold off;

   
%% d: Combined Network
% subplot(NRow,NCol,cellM{3});
cla;
% Left ==> Right
Xs1c = [[Xs10(1,1:2);Xs10(2,1:2)]  [Xs10(1,1:2);Xs10(2,1:2)]+[s;0]];
Xu1c = [[Xu1(1,:); Xu1(2,:)] [Xu1(1,:)+s/2;-Xu1(2,:)+.5]];
conn1c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs1a,Xu1a,conn1a] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[5;1]);
% Right ==> Left
% [Xs1a2,Xu1a2,conn1a2] = tesselate_network_old(Xs1c2,Xu1c2,conn1c,[s;0],[2;1]);
Xs1a2 = -Xs1a(:,1:8) + [-s/2; .5];
Xu1a2 = -Xu1a(:,1:14) + [-s/2; .5];
Xs1a2 = Xs1a2 + [(size(Xs1a2,2)-2)*s/2;0];
Xu1a2 = Xu1a2 + [(size(Xs1a2,2)-2)*s/2;0];
Xs1a = Xs1a - [(size(Xs1a,2)-2)*s/2;0];
Xu1a = Xu1a - [(size(Xs1a,2)-2)*s/2;0];
conn1a2 = conn1a(find(conn1a(:,1)<=size(Xs1a2,2) &...
                      conn1a(:,2)<=(size(Xs1a,2)+size(Xu1a2,2))),:)
conn1a2(:,2) = conn1a2(:,2)-(min(conn1a2(:,2))-max(conn1a2(:,1)))+1;
% visualize_network(Xs1a,Xu1a,conn1a);
% visualize_network(Xs1a2,Xu1a2,conn1a2);
mV = max(conn1a);
mV2 = max(conn1a2);
conn1a = [conn1a+[0 mV2(1)]; conn1a2+[mV(1) mV(2)]];
[Xs1a,Xu1a,conn1a] = tesselate_network_old([Xs1a Xs1a2],[Xu1a Xu1a2],conn1a,[0;0],[1;1]);
R = rigidity([Xs1a,Xu1a],conn1a);
A = null(R');
visualize_network(Xs1a,Xu1a,conn1a(abs(A)<1e-5,:));


%% e: Propagate Combined Motion
[XC,fC] = sim_motion(Xs1a,Xu1a,conn1a,.01,5000,-[Xs1a Xu1a],0);
D1 = sqrt(squeeze(sum(diff(XC,1,2).^2)));
D1 = D1(1:size(Xs1a,2)-1,:);


%% Calculate Rigidity
figure(5); clf;
sV = zeros(1,size(XC,3));
sV = [];
VP = [];
VP2 = [];
for i = 1:5:size(XC,3)
    cla;
    R = rigidity(XC(:,:,i),conn1a);
%     histogram(svd(R),1000);
    [U sVP V] = svds(R,2,'smallest');
    sV = [sV sVP(1,1)];
    V = -5*V/sqrt(V'*V);
    VP = [VP V(:,1)];
    VP2 = [VP2 V(:,2)];
%     quiver(squeeze(XC(1,:,i)),squeeze(XC(2,:,i)),V(1:length(V)/2,1)',V([1:length(V)/2]+length(V)/2,1)',0);
    hold on;
    visualize_network(XC(:,1:size(Xs1a,2),i),...
                      XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),conn1a,.5);
    hold off;
    disp(sVP);
    drawnow;
end




%% Animation
fig = figure(4); clf;
fName = 'animation.gif';
dT = 0.03;
% dT = 0.1;
nXSh = 2.4;
nYSh = 2;
nSV = 1;
XdotV = XC(:,1:size(Xs1a,2),1); XdotV = XdotV-mean(XdotV,2);

for i = 1:40:size(XC,3)
    cla;
    dP = [D1(:,i)';D1(:,i)']; dP = dP(:);
    dPa = dP(1:end-1); dPb = [1;dP(3:end)];
    plot(d1,d2,'k-');
    hold on;
%     [Us, Uu, ~] = construct_motion(XC(:,1:size(Xs1a,2),i),...
%                                    XdotV,...
%                                    XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),...
%                                    conn1a, 0, 0);
%     Us = 0.3*Us/sqrt(sum(sum(Us.^2))); 
%     Uu = 0.3*Uu/sqrt(sum(sum(Uu.^2)));
%     XdotV = Us;
%     quiver(XC(1,:,i)/10 + nXSh, XC(2,:,i)/10 + nYSh,[Us(1,:) Uu(1,:)], [Us(2,:) Uu(2,:)],0);
    plot([0 4],[0 4],'LineStyle','--','color',[0 0 0 .5]);
    line(dP(1:end-1),[1;dP(3:end)]);
    visualize_network(XC(:,1:size(Xs1a,2),i)/20 + [nXSh;nYSh],...
                      XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i)/20 + [nXSh;nYSh],conn1a,.5);
    axis([1.4 3.53 1.4 3.5]);
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


%% Animate Network Alone
fig = figure(5); clf;
fName = 'animation_net.gif';
dT = 0.03;
% dT = 0.1;
nXSh = 1.6;
nYSh = 2;
nSV = 1;
XdotV = XC(:,1:size(Xs1a,2),1); XdotV = XdotV-mean(XdotV,2);

for i = 1:40:size(XC,3)
    cla;
    dP = [D1(:,i)';D1(:,i)']; dP = dP(:);
    dPa = dP(1:end-1); dPb = [1;dP(3:end)];
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


