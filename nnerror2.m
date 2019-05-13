%NNERROR2 Exact expected NN error from a dissimilarity matrix (2)
%
%   E = NNERROR2(D,M)
%   E = D*NNERROR2([],M)
%
% INPUT
%   D 	NxN dissimilarity dataset
%   M 	Vector with desired number of objects per class to be selected
%
% OUTPUT
%   E 	Expected NN errror
%
% DESCRIPTION
% Exact computation of the expected 1-nearest neighbor error. D should be a
% dataset containing a labeled square dissimilarity matrix of a set of
% objects related to itself. It is not needed that D is symmetric.
% So NNERROR2(D) ~= NNERROR2(D') is possible.
%
%   E = NNERROR2(D)
%   E = D*NNERROR2
%
% In this case a set of training set sizes is used to produce a full
% learning curve. E can be plotted by PLOTE.
%
% There is a similar routine NNERROR1 which is based on the random
% selection of objects, instead of class-wise.
%
% SEE ALSO
% DATASETS, NNERROR1, NNE, TESTKD

% Copyright: R.P.W. Duin, r.duin@ieee.org
% and Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester

function e = nnerror2(D,n)

if nargin < 2, n = []; end
if nargin < 1 || isempty(D)
  e = prmapping(mfilename,'fixed',n);
  e = setname(e,'Exp_NN_error2');
else
  isdataset(D);
	[m,m2] = size(D);
  if m ~= m2,
    error('Distance matrix should be square. Use testkd.');
  end

	if isempty(n)
		% find full curve, but gain some speed
    mc = min(classsizes(D)); % size of smallest class
		L = [1:20 22:2:40 45:5:60 70:10:100 120:20:300 350:50:1000 1100:100:10000];
		L = [L(find(L<mc-1)) mc-1];
		f = zeros(1,mc-1);
		prwaitbar(max(L),'Compute Learning Curve')
		for i=1:length(L)
			prwaitbar(max(L),L(i));
			f(L(i)) = feval(mfilename,D,L(i));
			if (i > 1) & (L(i)-L(i-1) > 1)
				for n=L(i-1):L(i)
					f(n) = f(L(i-1)) + (f(L(i))-f(L(i-1)))*(n-L(i-1))/(L(i)-L(i-1));
				end
			end
		end
		prwaitbar(0)
    
    e.error = f;
    e.xvalues = [1:length(e.error)];
    e.title = 'Learning curve 1-NN rule';
    e.xlabel = 'Size training set';
    e.ylabel = 'Expected classification error';
    e.plot = 'semilogx';
    
  else
    
    isdataset(D);
    %discheck(D);
    nlab     = getnlab(D);
    [m,mm,c] = getsize(D);

    % Number of objects per class
    nc = classsizes(D);

    % Compute for all training set sizes
    if nargin < 2
      n = min(nc)-1;
    end

    if length(n) > 1
      e = zeros(1,length(n));
      s = sprintf('Running over %i training set sizes: ',length(n));
      prwaitbar(length(n),s)
      for j=1:length(n)
        prwaitbar(length(n),j,[s int2str(j)]);
        e(j) = nnerror(D,n(j));
      end
      prwaitbar(0);
    else
      if n >= min(nc)
        error('Requested size of the training set is too large.')
      end
      % Call for the given sample size	
      [D,I] = sort(D);
      I     = reshape(nlab(I),m,m);   % Order objects according to their distances
      ee    = zeros(1,m);

      % Loop over all classes
      prwaitbar(c*c,'Computing exact 1-NN LOO errors') 
      for j = 1:c
        % Find probabilities Q that objects of other classes are not selected
        Q = ones(m,m);
        for i = 1:c
          prwaitbar(c*c,(j-1)*c+i,'Computing exact 1-NN LOO errors')
          if i~=j
            [p,q] = nnprob(nc(i),n); 
            q = [1,q];
            C = cumsum(I==i)+1;
            Q = Q.*reshape(q(C),m,m);
          end
        end

        % Find probabilitues P for objects of this class to be the first
        [p,q] = nnprob(nc(j),n); 
        p = [0,p];
        C = cumsum(I==j)+1;
        P = reshape(p(C),m,m);

        % Now estimate the prob EC it is really the NN
        J     = find(I==j);
        EC    = zeros(m,m);
        EC(J) = P(J).*Q(J);

        % Determine its error contribution
        L     = find(nlab==j);
        ee(L) = 1-sum(EC(2:m,L))./(1-EC(1,L));	% Correct for the training size
      end
      prwaitbar(0);

      % Average for the final result
      e = abs(mean(ee));
    end
  end
end
return

%NNPROB Probability of selection as the nearest neighbor
%
% [P,Q] = NNPROB(M,K)
%
% If K objects are selected out of M, then P(i) is the probability
% that the i-th object is the nearest neigbor and Q(i) is the probability
% that this object is not selected.

function [p,q] = nnprob(m,k)
p = zeros(1,m);
q = zeros(1,m);
q(1) = (m-k)/m;
p(1) = k/m;
for i=2:(m-k+1)
	q(i) = q(i-1)*(m-k-i+1)/(m-i+1);
	p(i) = q(i-1)*k/(m-i+1);
end

