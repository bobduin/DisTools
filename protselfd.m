%PROTSELFD Forward Prototype Selection for Dissimilarity Matrices
%
%   [W,E,KOPT] = PROTSELFD(D,CRIT,K)
%    W  = D*PROTSELFD([],CRIT,K)
%
% INPUT
%   D     Dataset, dissimilarity matrix, not necessarily square.
%   K     Integer, desired number of prototypes. Default: select all, just
%         rank the columns of D
%   CRIT  'super' supervised selection using 1NN error on prototypes.
%         'nn'    same as 'super'
%         'loo' - supervised selection using leave-one-out 1NN error estimation.
%         'maxdist' - unsupervised selection minimizing the maximum
%                     distance to the nearest prototype.
%         'meandist' - unsupervised selection minimizing the average
%                      distance to the nearest prototype.
%
% OUTPUT
%   W     Selection mapping ('feature selection'). +W: Prototype indices.
%   E     Error stimate as a function of number of selected prototypes
%         (for supervised selection only reliable for prototype sizes >= class size)
%   KOPT  Estimate for best size in avoiding peaking 
%         (supervised selection only)
%
% DESCRIPTION
% This procedure for optimizing the representation set of a dissimilarity 
% matrix is based on a greedy, forward selection of prototypes. This
% implies that for larger values of K always the selected features always
% start with the ones found for smaller K.
%
% In case of a supervised selection, D should be a labeled dataset with 
% prototype labels stored as feature labels. The 1NN error to the nearest
% prototype based on the given dissimilarities is used as a criterion. 
% In case of a leave-one-out error estimation it is assumed that the first
% objects in D (rows) correspond with the objects used for the 
% representation set (columns).
% 
% For K=1 just a single prototype has to be returned, but computing the
% 1NN error is not possible as all objects are assigned to the same class.
% In that case the centre object of the largest class will be returned.
%
% Note that the search continues untill K prototypes are found. This might
% be larger than desired due to peaking (overtraining). Therefor an
% estimate for the optimal number of prototype is returned in KOPT. 
%
% The prototype selection may be applied by C = B*W(:,1:KSEL), in which B
% is a dissimilarity matrix based on the same representation set as A (e.g.
% A itself) and C is a resulting dissimilarity matrix in which the KSEL
% (e.g. KOPT) best prototypes are selected.
%
% In case of unsupervised selection the maximum or the mean distances to
% the nearest prototype are minimized. These criteria are the same as used
% in the k-centers and k-mediods cluster procedures.
%
% A random selection just selects at random a set of columns of D.
%
% REFERENCE
% E. Pekalska, R.P.W. Duin, and P. Paclik, Prototype selection for
% dissimilarity-based classification, Pattern Recognition, 
% vol. 39, no. 2, 2006, 189-208.
%
% SEE ALSO
% DATASETS, MAPIINGS, KNNDC, DISEX_PROTSELFD

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [R,e,D] = protselfd(varargin)


		 
  varargin = shiftargin(varargin,{'char','prmapping'});
  argin = setdefaults(varargin,[],'super',[]);
  if mapping_task(argin,'definition')
    R = define_mapping(argin,'untrained','Forward Prototype Sel');
    return
  end
    
  [D,crit,ksel] = deal(argin{:});
  k = size(D,2);
  if isempty(ksel), ksel = k; end
  switch lower(crit)
    case {'loo','super',''}
      [R,e,D,J,nlab,clab] = protselfd_super_init(D,crit,ksel);
    case {'maxdist','meandist'}
      [R,e] = protselfd_unsuper(D,crit,ksel);
    case {'random'}
      N = randperm(size(D,2));
      N = N(1:ksel);
      R = featsel(k,N);
    otherwise
      error('Unknown selection criterion')
  end
      
return  


function [R,e,D,J,nlab,clab] = protselfd_super_init(D,crit,ksel)
% this routine takes care of the initialisation of supervised selection

  isdataset(D);
  [m,k,c] = getsize(D);
  if isempty(ksel), ksel = k; end
  if strcmp(crit,'loo') | strcmp(crit,'LOO')
    if k > m
      error('More rows than columns expected for dissimilarity matrix')
    end
    issquare(D(1:k,:));
    D(1:k,:) = D(1:k,:) + 1e100*eye(k); % get rid of diagonal for LOO
  end
  
  %Initialise by the centre of the largest class
  cc = classsizes(D);
  [cmax,n] = max(cc); % n is the largest class
  lablist = getlablist(D);
  nlab = getnlab(D);
  clab = renumlab(getfeatlab(D),lablist);
  R = find(nlab == n);
  C = find(clab == n);
  dd = +D(R,C);
  [dmin,rmin] = sort(dd,1); % find one but most remote object
  [dmin,cmin] = min(dmin(end-1,:)); % find prototype for which this is minimum
  R = C(cmin);
  
  e = zeros(1,ksel);
  [nlab,clab] = renumlab(getlabels(D),getfeatlab(D));
  [dd,J] = min(+D(:,R),[],2);
  e(1) = sum(clab(R(J)) ~= nlab);
  
	if ksel > 1
		% this will be a deep recursive call !!!
		prwaitbar(ksel,'Forward prototype selection')
		[R,e,D,J,nlab,clab] = protselfd_super(D,ksel,R,J,e,nlab,clab);
		prwaitbar(0);
	end
  e = e(1:length(+R))/m;
  R = featsel(k,R);
  
  % Find optimal number of prototypes in avoiding peaking
  
  Jopt = find(e==min(e));
  D = floor((Jopt(end)+Jopt(1))/2);
  
return
  
function [R,e,D,J,nlab,clab] = protselfd_super(D,ksel,R,J,e,nlab,clab)

  [m,k,c] = getsize(D);
  d = +D;
  S = [1:k];   % all candidates
  S(R) = [];   % exclude ones we have
  emin = inf;
  dmin = inf;
  r = length(R);
  prwaitbar(ksel,r);
  for j=S      % run over all candidates left
		% the following tricky statements finds the nearest neighobor indices n
		% for all objects to their nearest prototype (n=1) or the candidate 
		% prototype (n=2). In ds the minimum distances are stored and used for
		% solving ties later.
    [ds,n] = min([d(m*(R(J')'-1)+[1:m]'),d(:,j)],[],2);
		% the labels of the nearest prototypes and the candidates
    cclab = [clab(R(J)') repmat(clab(j),m,1)];
		% compute the nearest neighbor error using the computed n
    ee = sum(cclab(m*(n-1)+[1:m]') ~= nlab);
    de = sum(ds);
		% if better, use it
    if ee < emin || ((ee == emin) && (de < dmin))
      emin = ee;
      jmin = j;
      JJ = [J repmat(r+1,m,1)];
      Jmin = JJ(m*(n-1)+[1:m]');
      Rmin = [R jmin];
      dmin = de;
    end
  end

  if emin <= e(r) || 1 % we even continue if emin increases due to peaking
    e(r+1) = emin;
    R = Rmin;
    if (r+1) < ksel
			[R,e,D,J,nlab,clab] = protselfd_super(D,ksel,R,Jmin,e,nlab,clab);
    end
  end
  
return

%PROTSELFD_UNSUPER Forward prototype selection
%
%               N = PROTSELFD_UNSUPER(D,CRIT,P)
%
% INPUT
%   D     Square dissimilarity matrix, zeros on diagonal
%   P     Number of prototypes to be selected
%   CRIT  'maxdist' or 'meandist'
%
% OUTPUT
%   N     Indices of selected prototypes
%
% DESCRIPTION
% Sort objects given by square dissim matrix D using a greedy approach
% such that the maximum NN distance from all objects (prototypes)
% to the first K: max(min(D(:,N(1:K),[],2)) is minimized.
%
% This routines tries to sample the objects such that they are evenly
% spaced judged from their dissimilarities. This may be used as
% initialisation in KCENTRES. It works reasonably, but not very good.
%
% SEE ALSO
% KCENTRES, KMEDIODS

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [N,e] = protselfd_unsuper(d,crit,p)

d = +d;
[m,k] = size(d);
if nargin < 3 || isempty(p), p = k; end
if nargin < 2 || isempty(crit), crit = 'max'; end
L = 1:k;
N = zeros(1,p);
switch crit
  case 'maxdist'
    [dummy,n] = min(max(d));    % this is the first (central) prototype
  case 'meandist'
    [dummy,n] = min(mean(d));   % this is the first (central) prototype
end
e = d(:,n);                 % store here the distances to the nearest prototype (dNNP)
f = min(d,repmat(e,1,k));   % replace distances that are larger than dNNP by dNNP
N(1) = n;                   % ranking of selected prototypes
L(n) = [];                  % candidate prototypes (all not yet selected objects)

for j=2:p                   % extend prototype set
  switch crit               % select the next prototype out of candidates in L
    case 'maxdist'
      [dummy,n] = min(max(f(:,L))); 
    case 'meandist'
      [dummy,n] = min(mean(f(:,L))); 
  end
  e = min([d(:,L(n)) e],[],2);   % update dNNP
  f = min(d,repmat(e,1,k));      % update replacement of distances that are larger
                                 % than dNNP by dNNP
  N(j) = L(n);                   % update list of selected prototypes
  L(n) = [];                     % update list of candidate prototypes
end
N = featsel(k,N);

