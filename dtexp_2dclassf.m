%DTEXP_2DCLASSF  DisTools 2D classf example

delfigs                   % remove all figures
a = gendatb;              % generate a 2d dataset, 2 classes : 50 50
fs2ds = proxm([],'m',2);  % define the feature space to dissimilarity
                          % space mapping based on Minkowsky-2 distances
                          
% some 2D feature based classifiers
scatterd(a);              % scatterplot of the dataset
w1 = a*nmc;               % nearest mean classifier in feature space
w2 = a*knnc;              % kNN classifier in feature spae
w3 = a*libsvc([],proxm([],'p',3),1); % SVM on polynomial-3 kernel
plotc({w1,w2,w3})         % plot all classifiers

% Fisher in disspace [1 1]
w4 = gendat(a,[1 1])*fs2ds; % define disspace on 2 prototypes
w4 = w4*fisherc;            % give it fisherc as classifier
w4 = a*setname(w4,'fdsc1'); % give it a name and train it

% Fisher in disspace [2 2]
w5 = gendat(a,[2 2])*fs2ds; % define disspace on 4 prototypes
w5 = w5*fisherc;            % give it fisherc as classifier
w5 = a*setname(w5,'fdsc2'); % give it a name and train it

% Fisher in disspace on 10% of dataset as prototypes
w6 = gendat(a,0.1)*fs2ds;   % define disspace on 10% of datasets
w6 = w6*fisherc;            % give it fisherc as classifier
w6 = a*setname(w6,'fdsc-0.1'); % give it a name and train it

% Fisher in PCA-0,99 reduced disspace 
w7 = fs2ds*pcam([],0.99);     % define disspace on 10% of datasets
w7 = w7*fisherc;             % give it fisherc as classifier
w7 = a*setname(w7,'fdsc-pca'); % give it a name and train it

% create another figure and plot these classifiers
figure;
scatterd(a);
plotc({w4,w5,w6,w7})

% have a look at the 2D dissimilarity based classifier
parsc(w7)                   