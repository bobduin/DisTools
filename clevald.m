%CLEVALD Classifier evaluation (learning curve) for dissimilarity data
%
%   E = CLEVALD(D,CLASSF,TRAINSIZES,REPSIZE,NREPS,T,TESTFUN)
%
% INPUT
%   D          Square dissimilarity dataset
%   CLASSF     Classifiers to be evaluated (cell array)
%   TRAINSIZE  Vector of training set sizes, used to generate subsets of D
%              (default [2,3,5,7,10,15,20,30,50,70,100]). TRAINSIZE is per
%              class unless D has no priors set or has soft labels.
%   REPSIZE    Representation set size per class (>=1), or fraction (<1)
%              (default total training set). If REPSIZE > TRAINSIZE
%              testobjects are included for representation.
%   NREPS      Number of repetitions (default 1)
%   T          Test dataset (default [], use remaining samples in A)
%   TESTFUN    Mapping,evaluation function (default classification error)
%
% OUTPUT
%   E          Error structure (see PLOTE) containing training and test
%              errors
%
% DESCRIPTION
% Generates at random, for all class sizes defined in TRAINSIZES, training
% sets out of the dissimilarity dataset D. The representation set is either
% equal to the training set (REPSIZE = []), or a fraction of it (REPSIZE <1)
% or a random subset of it of a given size (REPSIZE>1). This set is used
% for training the untrained classifiers CLASSF. The resulting trained
% classifiers are tested on the training objects and on the left-over test 
% objects, or, if supplied, the testset T. This procedure is then repeated 
% NREPS times. The default test routine is classification error estimation 
% by TESTC([],'crisp'). 
%
% The returned structure E contains several fields for annotating the plot
% produced by PLOTE. They may be changed by the users. Removal of the field
% 'apperror' (RMFIELD(E,'apperror')) suppresses the draw of the error
% curves for the training set.
%
% Training set generation is done "with replacement" and such that for each
% run the larger training sets include the smaller ones and that for all
% classifiers the same training sets are used.
%
% This function uses the RAND random generator and thereby reproduces 
% if its seed is reset (see RAND).
% If CLASSF uses RANDN, its seed should be reset as well.
%
% SEE ALSO
% MAPPINGS, DATASETS, CLEVAL, TESTC, PLOTE

% R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function e = clevald(a,classf,learnsizes,repsize,nreps,t,testfun)

	prtrace(mfilename);

  if (nargin < 7) | isempty(testfun)
    testfun = testc([],'crisp');
  end;
  if (nargin < 6) 
    t = [];
  end;
  if (nargin < 5) | isempty(nreps);
    nreps = 1;
  end;
  if (nargin < 4) 
    repsize = [];
  end
  if (nargin < 3) | isempty(learnsizes);
    learnsizes = [2,3,5,7,10,15,20,30,50,70,100];
  end;
  if ~iscell(classf), classf = {classf}; end

  % Assert that all is right.
  isdataset(a); issquare(a); ismapping(classf{1}); 
  if (~isempty(t)), isdataset(t); end
  
  if numel(repsize) > 1
    error('Size of the representation set should be a scalar')
  end

  % Remove requested class sizes that are larger than the size of the
  % smallest class.

	[m,k,c] = getsize(a);
	if ~isempty(a,'prior') & islabtype(a,'crisp')
		classs = true;
		mc = classsizes(a);
		toolarge = find(learnsizes >= min(mc));
		if (~isempty(toolarge))
			prwarning(2,['training set class sizes ' num2str(learnsizes(toolarge)) ...
									 ' larger than the minimal class size; removed them']);
			learnsizes(toolarge) = [];
		end
	else
		if islabtype(a,'crisp') & isempty(a,'prior')
			prwarning(1,['No priors found in dataset, class frequencies are used.' ...
			newline '            Training set sizes hold for entire dataset']);
		end
		classs = false;
		toolarge = find(learnsizes >= m);
		if (~isempty(toolarge))
			prwarning(2,['training set sizes ' num2str(learnsizes(toolarge)) ...
									 ' larger than number of objects; removed them']);
			learnsizes(toolarge) = [];
		end
	end
  learnsizes = learnsizes(:)';

  % Fill the error structure.

  nw = length(classf(:));
  datname = getname(a);

  e.n       = nreps;
  e.error   = zeros(nw,length(learnsizes));
  e.std     = zeros(nw,length(learnsizes));
  e.apperror   = zeros(nw,length(learnsizes));
  e.appstd     = zeros(nw,length(learnsizes));
  e.xvalues = learnsizes(:)';
	if classs
		e.xlabel   = 'Training set size per class';
	else
		e.xlabel   = 'Training set size';
	end
  e.names   = [];
  if (nreps > 1)
    e.ylabel= ['Averaged error (' num2str(nreps) ' experiments)'];
  elseif (nreps == 1)
    e.ylabel = 'Error';
  else
    error('Number of repetitions NREPS should be >= 1.');
  end;
  if (~isempty(datname))
    if isempty(repsize)
      e.title = [datname ', Rep. Set = Train Set'];
    elseif repsize < 1
      e.title = [datname ', Rep. size = ' num2str(repsize) ' Train size'];
    else
      e.title = [datname ', Rep. size = ' num2str(repsize) ' per class'];
    end
  end
  if (learnsizes(end)/learnsizes(1) > 20)
    e.plot = 'semilogx';        % If range too large, use a log-plot for X.
  end

  % Report progress.
	
	s1 = sprintf('cleval: %i classifiers: ',nw);
	prwaitbar(nw,s1);

  % Store the seed, to reset the random generator later for different
  % classifiers.

	seed = rand('state');

  % Loop over all classifiers (with index WI).

  for wi = 1:nw
		
    if (~isuntrained(classf{wi}))
      error('Classifiers should be untrained.')
    end
    name = getname(classf{wi});
    e.names = char(e.names,name);
    prwaitbar(nw,wi,[s1 name]);

    % E1 will contain the error estimates.

    e1 = zeros(nreps,length(learnsizes));
    e0 = zeros(nreps,length(learnsizes));

    % Take care that classifiers use same training set.

    rand('state',seed); seed2 = seed;

		% For NREPS repetitions...
		
		s2 = sprintf('cleval: %i repetitions: ',nreps);
		prwaitbar(nreps,s2);

		for i = 1:nreps
	
			prwaitbar(nreps,i,[s2 int2str(i)]);
      % Store the randomly permuted indices of samples of class CI to use in
      % this training set in JR(CI,:).

			if classs

				JR = zeros(c,max(learnsizes));
			
				for ci = 1:c

					JC = findnlab(a,ci);

					% Necessary for reproducable training sets: set the seed and store
					% it after generation, so that next time we will use the previous one.
					rand('state',seed2);

					JD = JC(randperm(mc(ci)));
					JR(ci,:) = JD(1:max(learnsizes))';
					seed2 = rand('state'); 
				end
				
			elseif islabtype(a,'crisp')
				
				rand('state',seed2); % get seed for reproducable training sets
				% generate indices for the entire dataset taking care that in
				% the first 2c objects we have 2 objects for every class
				[a1,a2,I1,I2] = gendat(a,2*ones(1,c));
				JD = randperm(m-2*c);
				JR = [I1;I2(JD)];
				seed2 = rand('state'); % save seed for reproducable training sets
				
			else  % soft labels
				
				rand('state',seed2); % get seed for reproducable training sets
				JR = randperm(m);
				seed2 = rand('state'); % save seed for reproducable training sets
				
			end

			li = 0;										% Index of training set.
			
			nlearns = length(learnsizes);
			s3 = sprintf('cleval: %i sizes: ',nlearns);
			prwaitbar(nreps,s3);
			
			for j = 1:nlearns 
				
				nj = learnsizes(j);
				
				prwaitbar(nlearns,j,[s3 int2str(j) ' (' int2str(nj) ')']);
				li = li + 1;

        % J will contain the indices for this training set.

        J = [];
        R = [];
				
				if classs
					for ci = 1:c
						J = [J;JR(ci,1:nj)'];
						if isempty(repsize)
							R = [R JR(ci,1:nj)];
						elseif repsize < 1
							R = [R JR(ci,1:ceil(repsize*nj))];
						else
							%R = [R JR(ci,1:min(nj,repsize))];
							R = [R JR(ci,1:repsize)];
						end
					end;
				else
					J = JR(1:nj);
					if isempty(repsize)
						R = J;
					elseif repsize < 1
						R = JR(1:ceil(repsize*nj));
					else
						%R = JR(1:min(nj,repsize));
						R = JR(1:repsize);
					end
				end;
					
				trainset = a(J,R);
				trainset = setprior(trainset,getprior(trainset,0));
				w = trainset*classf{wi}; 					% Use right classifier.
				e0(i,li) = trainset*w*testfun;
				if (isempty(t))
  				Jt = ones(m,1);
					Jt(J) = zeros(size(J));
					Jt = find(Jt); 								% Don't use training set for testing.
					testset = a(Jt,R);
					testset = setprior(testset,getprior(testset,0));
					e1(i,li) = testset*w*testfun;
				else
					testset = t(:,R);
					testset = setprior(testset,getprior(testset,0));
					e1(i,li) = testset*w*testfun;
				end

			end
			prwaitbar(0);

		end
		prwaitbar(0);

    % Calculate average error and standard deviation for this classifier
    % (or set the latter to zero if there's been just 1 repetition).

		e.error(wi,:) = mean(e1,1);
		e.apperror(wi,:) = mean(e0,1);
		if (nreps == 1)
			e.std(wi,:) = zeros(1,size(e.std,2));
			e.appstd(wi,:) = zeros(1,size(e.appstd,2));
		else
			e.std(wi,:) = std(e1)/sqrt(nreps);
			e.appstd(wi,:) = std(e0)/sqrt(nreps);
		end
	end
	prwaitbar(0);

	% The first element is the empty string [], remove it.
	e.names(1,:) = [];

return
