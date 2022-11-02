function [L, fitdHCP, baseVec] = fitHeadHolder(params, HCP, PAP, normVec, sizeInf)
% Fit the Head Holder
x = params(1);
y = params(2);
theta = mod(params(3), 2*pi);

P1 = orthPos([x; y], PAP, normVec);
baseVec1 = orthVec(theta, normVec, 1);
baseVec2 = cross(baseVec1, normVec);
baseVec = [baseVec1, baseVec2];

W = [0, 1, 1, 0;
     0, 0, 1, 1];

fitdHCP = P1 + ([sizeInf.HH.lenH, sizeInf.HH.widH].*baseVec)*W;

L = sum((fitdHCP - HCP).^2, 'all', 'omitnan'); % Residual Sum of Squares



end

