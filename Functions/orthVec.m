function vecXYZ = orthVec(theta, normVec, normON)
% orthogonal vector
vecXY = [cos(theta); sin(theta)];
vecXYZ = [vecXY; sum(-vecXY.*normVec(1:2, :), 1, 'omitnan')./normVec(3, :)];

if normON
    vecXYZ = vecXYZ./vecnorm(vecXYZ, 2, 1);
end

end

