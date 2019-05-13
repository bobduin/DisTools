%DTEXP_CLASSF  DisTools 2D classf example

if (exist('d') ~= 1)
  d = (gendatb([200 200])*mink2).^4;
end

fprintf('Size: %4.0f %4.0f\n',size(d));
fprintf('Classizes: ');
fprintf('%4.0f ',classsizes(d));
fprintf('\n')
w = d*pe_em;
fprintf('Signature: %4.0f %4.0f\n',signature(w));
fprintf('Negative Eigen-Fraction: %5.3f\n',w*nef);
fprintf('Non-Metric Fraction: %5.3f\n',d*nmf([],100));
fprintf('Asymetry: %5.3f\n',d*asymmetry);
fprintf('LOO NN error %5.3f\n',d*nne)

delfigs
figure; scatterd(d*pcam(d,2));
figure; scatterd(d*pe_em(d,[2 0]));
figure; scatterd(d*pe_em(d,[1 1]));
%options.inspect = 0;
options.st = 0;
figure; scatterd(d*mds(d,[],options));
figure; imagesc(+d);
figure; plotspectrum(w)
showfigs