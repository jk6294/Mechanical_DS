load quadrifolium.mat;

% Corrections
% Rotate and flip
R = rotz(90); R = R(1:2,1:2);
for i = 1:size(XCcc,3)
    XCcc(:,:,i) = R*XCcc(:,:,i);
end
XCcc(2,:,:) = -XCcc(2,:,:);

% Compress
XCcc = cat(3,XCcc(:,:,1:2:end-50),...
             XCcc(:,:,end-49:end),...
             repmat(XCcc(:,:,end),[1 1 100]));

% Animate
figure(1); clf;
animate_network(XCcc, conncc, 'figwidth',40, 'nframe',size(XCcc,3),...
                'nu',diff(max(conncc)),'scolor',[070  190 255]/255,...
                'ucolor',[50  90 200]/255,'save',1,'lwidth',.7,...
                'msize',5);