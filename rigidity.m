function [R] = rigidity(Xs,conn)
% Function to construct rigidity matrix
%
% Input: Xs, conn
%
% Output: R
%

% Get Parameters
d = size(Xs,1);
n = size(Xs,2);
I = zeros(size(conn,1), n);
for i = 1:size(conn,1)
    I(i,conn(i,1)) =  1;
    I(i,conn(i,2)) = -1;
end
if(d==2)
    R = [(Xs(1,conn(:,1)) - Xs(1,conn(:,2)))' .* I,...
         (Xs(2,conn(:,1)) - Xs(2,conn(:,2)))' .* I];
elseif(d==3)
    R = [(Xs(1,conn(:,1)) - Xs(1,conn(:,2)))' .* I,...
         (Xs(2,conn(:,1)) - Xs(2,conn(:,2)))' .* I,...
         (Xs(3,conn(:,1)) - Xs(3,conn(:,2)))' .* I];
end