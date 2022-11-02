function outputPonP = projectP2PonP(inputP, drctVec, inputPonP, normVec)
% project Points to Points on a Plane
outputPonP = inputP + (normVec'*(inputPonP - inputP)).*drctVec./(normVec'*drctVec);

end