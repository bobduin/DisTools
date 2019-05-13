%NMF Compute Non-Metricity Fraction of a dissimilarity matrix
%
%   [F,C] = NMF(D,N)
%   [F,C] = D*NMF
%
% INPUT
%   D   - Dissimilarity matrix or dataset
%   N   - Scalar integer for a fast estimation of F from N random
%         triangles. Default: use all.
%
% OUTPUT
%   F   - Non-metricity fraction
%   C   - Constant to make F=0, if added to the off-diagonal elements of D.
%
% DESCRIPTION
% This routine counts the number of triangle inequality violations D to
% estimate the NEF of D. If D is asymmetric, it is made symmetric by 
% averaging. F is the fraction of violating inequalities. Hence, for F = 0, 
% the triangle inequality holds for all triplets of the objects that
% construct the distance matrix D. 
%
% Additionally, a constant C is determined, which when added to the
% off-diagonal elements of D will take care that F = 0:
%
% If N is given, a random NxN subset of D is taken to speed up and F is an
% estimate. N=100 yields already a reasonably fast and accurate estimate.

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [f,c] = nmf(d,nmax)

if nargin < 2, nmax = []; end
if nargin < 1 || isempty(d)
  f = prmapping(mfilename,'fixed',nmax);
  return
end

d = +d;
[n,m] = size(d);
if isempty(nmax)
	nmax = n; % check all
end
	
if ~issym(d,1e-12)
  prwarning(1,'Distance matrix should be symmetric. It is made so by averaging.')
  d = 0.5*(d+d');
end

if any(d(1:n+1:end) ~= 0)
  prwarning(1,'Non-zero diagonal elements found. They are put to zero')
  d(1:n+1:end) = 0; 
end

c = 0;
ismetric = 1;

f = 0;         % Number of disobeyed triangle inequalities
if n > nmax     % take subset to speed up.
	R = randperm(n);
	R = R(1:nmax);
	d = d(R,R);
	n = nmax;
end
	
prwaitbar(n,'check triangle inequalities');
N = 0;
for j=1:n-1
  prwaitbar(n,j);
  r1 = d(:,j);
  r2 = d(j,j+1:n);
  % M checks the triangle ineqaulities; all elements of M are
  % positive or zero if this holds.
  M  = repmat(r1,1,n-j) + d(:,j+1:n) - repmat(r2,n,1);

  % The elemnets of M should ideally be compared to 0. However,
  % due to numerical inaccuracies, a small negative number should
  % be used. Experiments indicate that -1e-13 works well.
  % Otherwise, one gets wrong results.
  tol = 1e-13;
  mm  = sum(M(:) < -tol);

	N = N + (n-j)*n;
  f = f + mm;
  ismetric = ismetric & (mm == 0);
  if ~ismetric,
    cc = max(max(M));
    if cc > c,
      c = cc;
    end
  end
end
f = f/N;
prwaitbar(0);
return;
