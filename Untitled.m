clear; clc;

nT = 50;
% % Gradient along red/blue, sweep red
% RLin = linspace(0,1,nT);
% GLin = linspace(0,1,nT);
% BLin = linspace(0,1,nT);
% [RLin,GLin,BLin] = meshgrid(RLin,GLin,BLin);
% RLin = RLin(:)';
% GLin = GLin(:)';
% BLin = BLin(:)';
% CP = [RLin; GLin; BLin]';

% HLin = linspace(0,1,400);
% SLin = linspace(0,1,100);
% VLin = linspace(1,1,1);
% [HLin,SLin,VLin] = meshgrid(HLin,SLin,VLin);
% HLin = HLin(:)';
% SLin = SLin(:)';
% VLin = VLin(:)';
% CP = hsv2rgb([HLin; SLin; VLin]');
% RLin = SLin.*cos(HLin*2*pi);
% GLin = SLin.*sin(HLin*2*pi);
% BLin = VLin;

C1a = [231 178 072]/255;
C1b = [140 181 063]/255;
C1c = [047 086 151]/255;
CP1 = interp1([0 1],[C1;C2;C3],linspace(0,1,nT));

C1 = [227 036 055]/255;
C2 = [031 172 204]/255;
DL = [1 0];
DLv = linspace(0,1,1000);
CP = interp1(DL,[C1;C2],DLv);




figure(5); clf;
scatter(DLv,ones(1,length(DLv)),500,CP,'filled','marker','s');
% scatter3(RLin,GLin,BLin,100,CP,'filled');