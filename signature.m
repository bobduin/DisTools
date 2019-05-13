%SIGNATURE Retrieve the signature from a pseudo-Euclidean dataset or mapping
%
%   SIG = SIGNATURE(W)
%   SIG = SIGNATURE(A)
%   SIG = W*SIGNATURE
%   SIG = A*SIGNATURE
%   [P,Q] = SIGNATURE(W)
%   .....
%
% INPUT
%   W    PE mapping, W = PE_EM(D), if D is a dissimilarity matrix
%   A    Dataset, vectors in PE space, A = D*W
%
% OUTPUT 
%   SIG  Signature, 2-element vector [P,Q] with numbers of dimensions of
%        the positive and the negative subspaces of the PE space.
%   P    Number of dimensions of the positive space.
%   Q    Number of dimensions of the negative space.
%
% SEE ALSO
% DATASETS, MAPPINGS, PE_EM, PE_DISTM

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [sig,sig2] = signature(a)

  if nargin < 1 || isempty(a)
    sig = prmapping(mfilename,'combiner');
  else
    if isdataset(a)
      sig = getuser(a,'pe_signature');
      if isempty(sig)
        sig = [size(a,2) 0];
      end
    elseif ismapping(a)
      if ispe_em(a)
        a = pe_em(a);
      end
      if ispe_em(a)
        sig = getdata(a,'sig');
      else
        sig = [size(a,2) 0];
      end
    else % doubles
      sig = [size(a,2) 0];
    end
    if nargout == 2
      sig2 = sig(2);
      sig = sig(1);
    end
  end
return
  
  
  