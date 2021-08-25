% Figure 5: Dynamical Systems
%% Prepare Space
clear; clc;
params_fig;
set(groot,'defaulttextinterpreter','latex');


%% Figure dimensions
fig = figure(5); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [17.8  12.0];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  6.00 17.80  6.00]
        [ 0.00  0.00 17.80  6.00]];
    
% Fontsize
FS = 10;
% Distance visualization parameters
lw = .5;
gr = 0.9;
dgr = .6;
% Scaling parameters
scc = .25;

% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.18;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');


%% Superstability simulations
% Define Unit
s = sqrt(3);
Xs = [-s/2 0 s/2;...
      -1/2 1 -1/2];
Xu = [-s/2 -s/2; sqrt(s^2-s^2/4) -sqrt(s^2-s^2/4)];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Simulation
[XC,fC] = sim_motion10(Xs,Xu,conn,.01,190,[Xs,Xu],0);
D = squeeze(sqrt(sum(diff(XC(:,[1 2 3],:),1,2).^2)));
% Correct rotation
for i = 1:size(XC,3)
    Rz = rotz(atan2d(XC(2,3,i)-XC(2,1,i),XC(1,3,i)-XC(1,1,i)));
    XC(:,:,i) = Rz(1:2,1:2)'*XC(:,:,i);
end

% Combine
nR = 10;
[Xsa,Xua,conna] = network_chain_x(Xs(:,1:2), Xu, ones(1,nR));


%% Combine 2
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

% Combine
[Xsa4,Xua4,conna4] = tesselate_network_old([Xsa2 Xsa3], [Xua2 Xua3],...
                     [conna2(:,1) conna2(:,2)+max(conna3(:,1));...
                      conna3(:,1)+max(conna2(:,1)) conna3(:,2)+max(conna2(:,2))], [0;0], [1;1]);

% Equilibrium bond lengths
Leq = sqrt(sum((Xsa4(:,conna4(:,1)) - Xua4(:,conna4(:,2)-size(Xsa4,2))).^2));

% Simulate
X041 = zeros(2,size(Xsa4,2)+size(Xua4,2));
X042 = X041;
x0mag = .01;
X041(1,1) = -x0mag;
X042(:,11) = [-1/2;sqrt(3)/2]*x0mag;
disp('Simulating lower branch');
[XCa4a, dXCa4a, ta, KEa, PEa] = sim_dynamics([Xsa4 Xua4]', X041', conna4,...
                             ones(size(X041,2),1), ones(size(conna4,1),1), 2600, Leq);
disp('Simulating upper branch');
[XCa4b, dXCa4b, tb, KEb, PEb] = sim_dynamics([Xsa4 Xua4]', X042', conna4,...
                             ones(size(X041,2),1), ones(size(conna4,1),1), 2600, Leq);
XCa4a = permute(XCa4a,[2 1 3]);
XCa4b = permute(XCa4b,[2 1 3]);
                        

%% a: bottom arm
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Axes
axsh = [.07;.1];
scax = .8;
sct = sRat(pInd)*.88;
PEasc = PEa / (x0mag^2/2) * scax;
tasc = ta / max(ta) * sct;
hold on;
CVal = lines(4);
% Axes
plot(tasc+axsh(1),PEasc+axsh(2))
plot([0 0 max(tasc)]+axsh(1), [scax 0 0]+axsh(2), 'k-');
tTck = 0:.2:1;
line(ones(2,1).*(tTck*sct+axsh(1)), [0;0.02].*ones(1,length(tTck))+axsh(2),'color','k');
line([0;0.02].*ones(1,2)+axsh(1), ones(2,1).*[0 scax].*ones(1,2)+axsh(2),'color','k');
for i = 1:length(tTck)
    text(tTck(i)*sct+axsh(1),axsh(2)-.05,num2str(tTck(i)),NVTexth{:});
end
text(axsh(1)-.02,0+axsh(2),'0',NVTextr{:});
text(axsh(1)-.02,scax+axsh(2),'1',NVTextr{:});
text(0,scax/2+axsh(2),'Potential/Total Energy', NVTexth{:},'rotation',90);
text(sct+axsh(1)+.08,axsh(2),'time',NVTexth{:});
hold off;
% Draw networks
smpt = [.2 .56 .92 .99];
[~,smpt] = min(abs(smpt - tasc/sct),[],1);
scn = .04;
shn = [0.10  0.00 -0.15  0.38;...
       0.50  0.50  0.60  0.22];
thv = [0 -2 -4 -4];
for i = 1:length(smpt)
    hold on;
    plot(tasc(smpt(i))+axsh(1),PEasc(smpt(i))+axsh(2),'o',...
         'markersize',4,'linewidth',4,'color',CVal(i,:));
    hold off;
    Rz = rotz(thv(i)); Rz = Rz(1:2,1:2);
    visualize_network(Rz*XCa4a(:,1:25,smpt(i))*scn+shn(:,i)+[tasc(smpt(i));0],...
                      Rz*XCa4a(:,26:end,smpt(i))*scn+shn(:,i)+[tasc(smpt(i));0],...
                      conna4,'ucolor',CVal(i,:),'ms',3);
end
texta = '\textbf{a} \hspace{1.5mm}Potential energy along an extension of the bottom arm';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
                        

%% b: top arm
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Axes
axsh = [.07;.04];
sct = sRat(pInd)*.88;
PEbsc = PEb / (x0mag^2/2) * scax;
tbsc = tb / max(tb) * sct;
hold on;
CVal = lines(4);
% Axes
plot(tbsc+axsh(1),PEbsc+axsh(2))
plot([0 0 max(tbsc)]+axsh(1), [scax 0 0]+axsh(2), 'k-');
tTck = 0:.2:1;
line(ones(2,1).*(tTck*sct+axsh(1)), [0;0.02].*ones(1,length(tTck))+axsh(2),'color','k');
line([0;0.02].*ones(1,2)+axsh(1), ones(2,1).*[0 scax].*ones(1,2)+axsh(2),'color','k');
for i = 1:length(tTck)
    text(tTck(i)*sct+axsh(1),axsh(2)-.05,num2str(tTck(i)),NVTexth{:});
end
text(axsh(1)-.02,0+axsh(2),'0',NVTextr{:});
text(axsh(1)-.02,scax+axsh(2),'1',NVTextr{:});
text(0,scax/2+axsh(2),'Potential/Total Energy', NVTexth{:},'rotation',90);
text(sct+axsh(1)+.08,axsh(2),'time',NVTexth{:});
hold off;
% Draw networks
smpt = [.2 .56 .92 .99];
[~,smpt] = min(abs(smpt - tbsc/sct),[],1);
scn = .04;
shn = [0.10  0.00 -0.15  0.38;...
       0.40  0.42  0.49  0.25];
thv = [1 3 4 2];
for i = 1:length(smpt)
    hold on;
    plot(tbsc(smpt(i))+axsh(1),PEbsc(smpt(i))+axsh(2),'o',...
         'markersize',4,'linewidth',4,'color',CVal(i,:));
    hold off;
    Rz = rotz(thv(i)); Rz = Rz(1:2,1:2);
    visualize_network(Rz*XCa4b(:,1:25,smpt(i))*scn+shn(:,i)+[tbsc(smpt(i));0],...
                      Rz*XCa4b(:,26:end,smpt(i))*scn+shn(:,i)+[tbsc(smpt(i));0],...
                      conna4,'ucolor',CVal(i,:),'ms',3);
end

textb = '\textbf{b} \hspace{1.5mm}Potential energy along an extension of the top arm';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);


%% Save
fName = 'suppfig_superstable_dynamics1';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');


