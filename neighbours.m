%%% Neighbourhood function for Moore neighbours (King's graph) %%%
%
% A (M by N)
%
%east       : idx + M
%south east : idx + M + 1
%south      : idx + 1
%south west : idx - M + 1
%west       : idx - M
%north west : idx - M - 1
%north      : idx - 1
%north east : idx + M -1
%
%Input: 
%   idx : Linear index
%   M   : number of rows of (padded) matrix (nrows(A)+2)
%   A   : Padded matrix of which to extract neighbours
%Output: Vector of linear indices corresponding to Moore neighbours of
%input index

function [NB, indices] = neighbours(idx, A)
M = size(A,1);
indices = [idx-1 idx+M-1 idx+M idx+M+1 idx+1 idx-M+1 idx-M idx-M-1];
NB = nonzeros(A(indices));
end

