function holeVrtx = getHole(inputP, drctVec, inputPonP, normVec, baseVec, r, n)
% get Hole's Vertices
theta = linspace(0, 2*pi*(1 - 1/n), n);

Vrtx = [];
for i = 1:size(inputP, 2)
    Vrtx = [Vrtx, inputP(:, i) + r.*(baseVec*[sin(theta); cos(theta)])]; % *** CCW
end

holeVrtx = projectP2PonP(Vrtx, drctVec, inputPonP, normVec);


end

