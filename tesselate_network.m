function [XsO, connO] = tesselate_network(Xs, conn, XTr, NTr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for efficiently tesselating networks in cardinal directions.
%
% Inputs
% Xs:       d x ns matrix of specified node coordinates
% conn:     s x 2 matrix of edges from node i to node j
% XTr:      d x 1 vector of coordinate translations per coordinate
% NTr:      d x 1 vector of number of replications per coordinate
%
% Outputs
% XsO:      d x NsO matrix of NsO tesselated specified node positions
% connO:    sO x 2 matrix of connections from tesslated nodes i to j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Parameters
d = size(Xs,1);             % Dimension of embedding
N = size(Xs,2);             % Number of specified nodes
s = size(conn,1);           % Number of edges

% Matrices of transformations
XsOProx = repmat(Xs, [1, prod(NTr)]);
connOProx = repmat(conn, [prod(NTr), 1]);

% Replicate
if(d==2)
    % Extend Coordinates and Connection Indices
    for i = 1:NTr(1)
        for j = 1:NTr(2)
            XsOProx(:,[1:N]+((i-1)*NTr(2)+(j-1))*N) = Xs + XTr.*[(i-1); (j-1)];
            connOProx([1:s]+((i-1)*NTr(2)+(j-1))*s,:) = conn + [N N]*(j-1 + (i-1)*NTr(2));
        end
    end
elseif(d==3)
    % Extend Coordinates and Connection Indices
    for i = 1:NTr(1)
        for j = 1:NTr(2)
            for k = 1:NTr(3)
                XsOProx(:,[1:N]+((i-1)*NTr(2)*NTr(3)+(j-1)*NTr(3)+(k-1))*N) = Xs + XTr.*[(i-1); (j-1); (k-1)];
                connOProx([1:s] +((i-1)*NTr(2)*NTr(3)+(j-1)*NTr(3)+(k-1))*s,:) = conn + [N nu]*(j-1 + (i-1)*NTr(2));
            end
        end
    end
end
% Filter Redundant Connections
[C,IA,IC] = uniquetol(XsOProx',1e-5,'ByRows',1);
[B,I] = sort(IA);
[~,J] = sort(I);
XsO = XsOProx(:,B);
connO = [J((IC(connOProx(:,1)))), J((IC(connOProx(:,2))))];
% connO(:,2) = connO(:,2) + max(connO(:,1));



end