function posXYZ = orthPos(posXY, PAP, normVec)
% orthogonal position
posXYZ = [posXY; sum([PAP.*normVec; -posXY.*normVec(1:2, :)], 1, 'omitnan')./normVec(3, :)];

% x = posXY(1, :);
% y = posXY(2, :);
% z = sum([PAP.*normVec; -[x; y].*normVec(1:2, :)])./normVec(3, :);
% posXYZ = [x; y; z];

end