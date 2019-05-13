%PE_SVC SVC for PE spaces
%
%   W = PE_SVC(A,C)
%   W = A*PE_SVC([],C)
%  
% INPUT
%   A	      Pseudo-Euclidean dataset
%   C       Trade_off parameter in the support vector classifier.
%           Default C = 1;
%  
% OUTPUT
%   W       Mapping: Support Vector Classifier
%
% DESCRIPTION  
% Computation of the linear SVC classifier for the Pseudo-Euclidean 
% dataset A.  Note that testsets should be defined in the same PE space 
% as A.
%
% Warning: class prior probabilities in A are neglected.
%
% EXAMPLE
% trainset = gendatm;
% testset = gendatm;
% Dtrain = trainset*proxm(trainset,'m',1);
% Dtest  = testset*proxm(testset,'m',1);
% w = pe_em(Dtrain);
% Xtrain = Dtrain*w;
% Xtest = Dtest*w;
% v = pe_svc(Xtrain);
% Xtest*v*testc
%
% SEE ALSO
% MAPPINGS, DATASETS, PE_EM

% R.P.W. Duin, r.duin@37steps.com

function w = pe_svc(a,c)

  if nargin < 2, c = 1; end

	if nargin == 0 || isempty(a)
		w = prmapping(mfilename,'untrained',{c});
		w = setname(w,'PE SVC');
		
	elseif ~ismapping(c) % training
		
		if ~ispe_dataset(a)
			prwarning(1,'Dataset is Euclidean')
			w = svc(a,[],c);
		else
			ktrain = pe_kernelm(a,a);
			v = svc(ktrain,0);
			w = prmapping(mfilename,'trained',{v,a},getlablist(a),size(a,2),getsize(a,3));
		end
		
	else % execution, testset is in a, trained mapping is in c
		
		w = getdata(c,1);
		trainset = getdata(c,2);
		ktest = pe_kernelm(a,trainset);
		w = ktest*w;
    
  end
		
return