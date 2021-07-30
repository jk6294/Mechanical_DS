function [Xsa,Xua,conna,C_UNa] = network_chain_x(Xsp,XuC,XuL,C_UN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create a network chain along the x-direction, left to right.
% Specifically for 2D networks with 3 specified node modules.
% Requires at least 4 modules
% 
% Inputs:
% Xsp:      2  x 2          matrix of specified node coordinates.
% XuC:      2  x 2 x k      matrix of unspecified node coordinates
% XuL:      1  x L          vector of unspecified node modules to use: 1-k
% C_UN:     k  x 4          optional matrix of unspecified node colors
% 
% Outputs:
% Xsa:      2  x L+2        matrix of chain specified node coordinates
% Xua:      2  x 2L         matrix of chain unspecified node coordinates
% conna:    6L x 2          matrix of connections
% C_UNa:    2L x 4          matrix of unspecified node colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Preliminary
k = size(XuC,3);        % Number of unique modules
L = size(XuL,2);        % Total chain length

% If unspecified node colors not provided, use default
if(nargin == 3)
    C_UN = repmat([100 100 255]/255,k,1);
end

% Calculate module shifts
Xsh = Xsp(1,2)-Xsp(1,1);
Ysh = Xsp(2,2)+Xsp(2,1);
XuC = cat(3,XuC, [XuC(1,:,:); -XuC(2,:,:)+Ysh]);


%% Place Unspecified Nodes
Xua = zeros(2,2*L);
C_UNa = zeros(2*L,3);
for i = 1:L
    % Node Positions
    Xua(:,[1,2]+2*(i-1)) = XuC(:,:,XuL(i)+k*mod(i-1,2))+[Xsh;0]*(i-1);
    % Node Colors
    C_UNa([1,2]+2*(i-1),:) = repmat(C_UN(XuL(i),:),2,1);
end


%% Place Specified Nodes
% X Coordinates
XsL1 = Xsh*ones(1,L+2);
% Y Coordinates
if(mod(L,2)==0)
    XsL2 = repmat(Xsp(2,:),1,(L+2)/2);
else
    XsL2 = [repmat(Xsp(2,:),1,(L+1)/2) Xsp(2,1)];
end
Xsa = [cumsum(XsL1)-Xsh+Xsp(1,1);XsL2];


%% Place Connections
connEs = repmat(3:L,6,1);
if(L>1)
    connEs = [1;1;2;2;2;2;connEs(:);repmat(L+1,4,1);repmat(L+2,2,1)];
    connEu = repmat((1:6)',1,L+2)+[2:2:2*L+4]-6;
    connEu = connEu(:); connEu = connEu(connEu>0 & connEu<=2*L)+L+2;
else
    connEs = [1;1;2;2;3;3];
    connEu = [4;5;4;5;4;5];
end
conna = [connEs connEu];

end
