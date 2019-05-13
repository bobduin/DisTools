%ASYMMETRY Estimate asymmetry of dissimilarity matrix
%
%   ASYM = ASYMMETRY(D)
%   ASYM = D*ASYMMETRY
%
% The asymmetry is defined as the average over 2*|D-D'|./(|D|+|D'|).

function asym = asymmetry(d)

  if nargin < 1 | isempty(d)
    asym = prmapping(mfilename,'fixed');
    return
  end

  d = +d;
  e = abs(d-d')./(abs(d)+abs(d')+1e-100); %avoid division by 0
  L = find(~isnan(e));
  asym = mean(mean(e(L)));

return
