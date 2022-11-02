function [bevelVrtx, normVec4Bevel, drctVec4Bevel, drctVec4Level] = getBevel(baseP, cntrP, normVec, sizeInf, sink)
% get Bevel's Vertices
baseVec1 = cntrP - baseP;
baseVec1 = baseVec1./norm(baseVec1, 2);
baseVec2 = cross(baseVec1, normVec);
drctVec4Bevel = cos(sizeInf.theta).*baseVec1 + sin(sizeInf.theta).*normVec;
normVec4Bevel = cross(baseVec2, drctVec4Bevel);

bevelVrtx = NaN(3, 4);

bevelVrtx(:, 2:3) = baseP - sink.*normVec4Bevel - sizeInf.len.*drctVec4Bevel + [-0.5, 0.5].*sizeInf.wid.*baseVec2;
bevelVrtx(:, [1, 4]) = projectP2PonP(bevelVrtx(:, 2:3), drctVec4Bevel, ...
    [0, 0, sizeInf.hig - sink]', [0, 0, 1]');

tempVec = bevelVrtx(:, 4) - bevelVrtx(:, 1);
tempVecLen = norm(tempVec, 2);
tempSin = sizeInf.wid./tempVecLen;
tempCos = sqrt(1 - tempSin.^2);
tempUnitVec1 = tempVec./tempVecLen;
tempUnitVec2 = cross([0, 0, 1]', tempUnitVec1);

drctVec4Level = tempCos.*tempUnitVec1 + tempSin.*tempUnitVec2;


end

