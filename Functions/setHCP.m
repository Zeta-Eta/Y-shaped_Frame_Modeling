function HCP = setHCP(initHCP, PAP, normVec)
% set Holes' Central Points
% normVec need to be unit vectors

HCP = NaN(size(initHCP));
for n = 1:size(normVec, 2)
    HCP(:, :, n) = projectP2PonP(initHCP(:, :, n), normVec(:, n), PAP(:, n), normVec(:, n));
%     HCP(:, :, n) = initHCP(:, :, n) + (normVec(:, n)'*(PAP(:, n) - initHCP(:, :, n))).*normVec(:, n);
end

% len = (normVec'*(PAP - initHCP))./norm(normVec, 2);
% but notice that: norm(normVec, 2) == 1 (normVec need to be unit vectors);

end

