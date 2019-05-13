%NEF Compute Negative Eigen Fraction from PE dataset or mapping
%
%   [F,R] = NEF(W)
%   [F,R] = NEF(A)
%   [F,R] = W*NEF
%   [F,R] = A*NEF
%
% INPUT
%   A   - Pseudo-Euclidean dataset
%   W   - Pseudo-Euclidean mapping
%
% OUTPUT
%   F   - Negative eigen-fraction
%   R   - Negative eigen-ratio
%
% DESCRIPTION
% The negative eigen-fraction is a ratio expressing the sum of 
% magnitudes of all negative eigenvalues divided by the sum of 
% magnitudes of all eigenvalues. The necessary numbers are stored
% in W = PE_EM(D) or in A = D*W, in which D is a symmetric dissimilarity
% matrix.
%
% The negative eigen-fraction is the ratio of the magnitude of the lowest
% negative eigenvalue and the largest positive eigenvalue.
%
% SEE ALSO
% DATASETS, MAPPINGS, PE_EM 

function [f,r] = nef(input)

  if nargin < 1 || isempty(input)
    f = prmapping(mfilename,'combiner');
  else
    eval = getspectrum(input);
    sig  = getsig(input);
    f = 1-sum(abs(eval(1:sig(1))))/sum(abs(eval));
    r = max(abs(eval(find(eval< 0)))) / max(eval(find(eval > 0)));
  end


