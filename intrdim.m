%INTRDIM Estimate Intrinsic dimension from dissimilarity data
%
%   K = INTRDIM(D)
%   K = D*INTRDIM
%
% INPUT
%   D	  NxN Dissimilarity dataset
%
% OUTPUT
%   K	  Estimated intrinsic dimension
%
% DESCRIPTION
% It is assumed that the data is generated from a normal distribution in
% K dimensions and that the Euclidean distance measure has been used to
% compute D. 
% Do NOT supply squared distances! They will definitely generate the wrong
% result.
%
% EXAMPLE
% A = GENDATB(100);
% D = SQRT(DISTM(A));
% K = INTRDIM(D);
%
% A = GENDATD(100,10);
% D = SQRT(DISTM(A*SCALEM(A,'VARIANCE')));
% K = INTRDIM(D);

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function k = intrdim(d)

if nargin < 1 || isempty(d)
  k = prmapping(mfilename,'fixed');
else
  d = +d.^2;
  [m,n] = size(d);

  % Assume selfreflecting dissim matrix
  if m==n & all(diag(+d) == 0)  
    d(1:n+1:n*n) = [];
  end
  % assume Chi square distribution and estimate dof
  k = round(2*(mean(d(:)).^2)/var(d(:))); 
  if k==0
    warning(['A zero intrinsic dimensionality has been found.' newline ...
      'Did you supply a proper, non-squared dissimilarity matrix?'])
  end
end