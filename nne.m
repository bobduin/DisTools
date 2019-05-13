%NNE Leave-one-out Nearest Neighbor Error on a Dissimilarity Matrix
%
%   [E,LAB] = NNE(D)
%   {E,LAB] = D*NNE
%
% INPUT
%   D   NxN symmetric dissimilarity dataset
%
% OUTPUT
%   E   Leave-one-out error
%   LAB Nearest neighbor labels
%
% DESCRIPTION
% Estimates the leave-one-out error of the 1-nearest neighbor rule
% on the square dissimilairy data D. D should be a dataset containing a
% labeled square dissimilarity matrix of a set of objects related to
% itself. Diagonal of D is assumed to be zero.
% It is not needed that D is symmetric. So NNE(D) ~= NNE(D') is possible.
%
% SEE ALSO
% DATASETS, NNERROR1, NNERROR2, TESTKD

% Copyright: Robert P.W. Duin, r.p.w.duin@prtools.org and
% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [e,NNlab] = nne(D)

if nargin == 0 || isempty(D)
  e = prmapping(mfilename,'fixed');
  e = setname(e,'LOO_NN_Error');
else
  isdataset(D);
  [m,n] = size(D);
  if m ~= n,
    error('Distance matrix should be square. Use testkd.');
  end

  lab = getlab(D);
  [nlab,lablist] = renumlab(lab);
  D(1:m+1:end)   = inf;
  [d,M] = min(D');
  e = mean(nlab(M) ~= nlab);
  if islabtype(D,'crisp')
    NNlab = lablist(nlab(M),:);
  else
    labs = gettargets(D);
    NNlab = labs(M,:);
  end
end
return;
