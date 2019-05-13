%PE_DISTM Square Pseudo-Euclidean (PE) Distance Between Two Datasets
%
%   D = PE_DISTM(A,B)
%   D = A*PE_DISTM([],B)
%   D = A*(B*PE_DISTM)
%   D = PE_DISTM(A)
%
% INPUT
%   A   NxK Matrix or dataset, Euclidean or pseudo-Euclidean
%   B   MxK Matrix or dataset, Euclidean or pseudo-Euclidean
%
% OUTPUT
%   D   NxM Euclidean distance dataset or matrix
%
% DESCRIPTION
% Computation of the square pseudo-Euclidean distance matrix D between two sets
% of vectors, A and B. The pseudo-Euclidean distance with the signature SIG
% (e.g. SIG = [10 5]) between vectors X and Y is computed as an indefinite
% 'Euclidean' distance:
%     D(X,Y) = (X-Y)'*J*(X-Y),
% where J is a diagonal matrix with 1's, followed by -1's.
% J = diag ([ONES(SIG(1),1);  -ONES(SIG(2),1)]);
%
% In a PE dataset the signature is stored in the user field, see
% SETSIG. This signature is derived from A. It is not stored in D as D
% does not contain vectors in a PE space.
%
% D is a dataset with the labels defined by the labels of A and feature labels 
% defined by the labels of B. 
%
% REMARKS
% Note that square pseudo-Euclidean distances can be negative. Zero
% distances can be negative due to small computational inaccuracies. The
% call 
%    D = PE_DISTM(A)
% forces symmetry and a zero diagonal.
% If A nd B are both Euclidean then D = PE_DISTM(A,B) calls D = DISTM(A,B)
%
% SEE ALSO
% PRDATASET, SETSIG, DISTM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function D = pe_distm(A,B,sym)

if nargin < 3, sym = 0; end % sym=1 forces symmetry

if nargin == 0;
  D = prmapping(mfilename,'untrained','train');
elseif nargin == 1
  D = feval(mfilename,A,A,1); % force symmetry
elseif ismapping(B)
  D = feval(mfilename,[],A);
elseif isempty(A)
  D = prmapping(mfilename,'fixed',B);
elseif isstr(B)
  D = prmapping(mfilename,'fixed',A);
else
  if ispe_dataset(A)
    siga = getsig(A);
    if ispe_dataset(B)
      sigb = getsig(B);
      if any(siga ~=sigb)
        error('Pseudo-Euclidean datasets should have the same signature')
      end
    end
    D = pedistm(A,B,siga);   % local computation, see below
  else
    if ispe_dataset(B)
      D = pedistm(A,B,getsig(B)); % local computation, see below
    else
      D = distm(A,B);
      sym = 0; % not needed here
    end
  end
  if sym % force symmetry
    D = (D+(+D)')/2;
    D = D - diag(diag(+D));
  end
end


function D = pedistm(A,B,sig)

  isdataset(A);

  a    = +A;
  b    = +B;
  [ra,ca] = size(a);
  [rb,cb] = size(b);

  if ca ~= cb,
    error ('The datasets should have the same number of features.');
  end


  J = [ones(1,sig(1))  -ones(1,sig(2))];
  D = - 2 .* a * diag(J) * b';
  D = D + ones(ra,1) * (J*(b'.*b'));
  D = D + (J * (a'.*a'))' * ones(1,rb);

  % Set object labels and feature labels
  D = setdata(A,D,getlab(B));
  D.name = 'Square PE dismat';
  if ~isempty(A.name)
    D.name = [D.name ' for ' A.name];
  end

return
  