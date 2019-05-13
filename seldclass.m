%SELDCLASS Select class subset from a square dissimilarity dataset
%
%   [DN,J] = SELDCLASS(D,C)
%   [DN,J] = D*SELDCLASS([],C)
%
% INPUT
%   A   NxN Dissimilarity Dataset
%   C   Indices of classes
%
% OUTPUT
%   DN  Subset of the dataset D
%   J   Indices of the selected objects
%
% DESCRIPTION
% The classes listed in C (indices) are extracted for the square
% dissimilarity matrix D by both, their rows (objects) as well as their
% columns (features). The indices in C refer to the label lis of D, to be
% retrieved by GETLABLIST(D).
%
% SEE ALSO
% DATASETS, SELCLASS, GETLABLIST, ISSQUARE, MAKESQUARE

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org, and
% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [D,J] = seldclass(D,n)

if nargin < 2, n = []; end
if nargin == 0 | isempty(D)
  D = prmapping(mfilename,'fixed',n); 
  return; 
end
issquare(D);

if isempty(n)
  return
end

J = [];
c = getsize(D,3);
if (any(n > c))
  error('Not that many classes')
end

for j=1:length(n)
  J = [J; findnlab(D,n(j))];
end
D = remclass(D(J,J));
return;
