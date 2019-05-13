%MAKESYM Make a dissimilarity matrix symmetric
%
%   [B,C] = MAKESYM(A,PROC)
%   [B,C] = A*MAKESYM(PROC)
%
% INPUT
%   A     Square dataset or matrix
%   PROC  String with procedure: 'average' (default), 'min', 'max',
%         'lower' (use lower triangle part) or 'upper' (user upper part).
%
% OUTPUT
%   B   Symmetric dataset or matrix computed as (A+A')/2
%   C   Asymmetric remaining part, dataset or matrix computed as (A-A')/2
%
% DESCRIPTION
% B is a symmetric matrix obtained combining A and A' by PROC.
%
%SEE ALSO
%DATASETS, ISSYM, ISSQUARE

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% EWI Faculty, Delft University of Technology and
% School of Computer Science, University of Manchester


function [d1,d2] = makesym(varargin)

varargin = shiftargin(varargin,'char');
argin = setdefaults(varargin,[],'average');
if mapping_task(argin,'definition')
  d1 = define_mapping(argin,'fixed');
  return
end

[d,proc] = deal(argin{:});
issquare(d);
switch lower(proc)
  case 'average'
    d1 = 0.5 * (d + d');
  case 'min'
    d1 = min(d,d');
  case 'max'
    d1 = max(d,d');
  case 'lower'
    e = tril(+d);
    d1 = setdata(d,e+e'-diag(diag(e)));
  case 'upper'
    e = triu(+d);
    d1 = setdata(d,e+e'-diag(diag(e)));
  otherwise
    error('Unkknown symmetrization procedure');
end
if nargout == 2,
  d2 = 0.5 * (d - d');
end
return
