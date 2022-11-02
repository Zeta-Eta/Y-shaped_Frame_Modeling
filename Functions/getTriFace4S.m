function TriFace = getTriFace4S(Vrtx, sidFace, edgeConON, figON, normalizationON)
% get Triangulation Face for Side
sidFaceNum = size(sidFace, 1);
TriFace = cell(sidFaceNum, 1);
if normalizationON
    for s = 1:sidFaceNum
        tempFace = sidFace(s, :);
        TriFace{s} = tempFace(triangulation4orthProj({Vrtx}, edgeConON, figON));
    end
else
    for s = 1:sidFaceNum
        tempFace = sidFace(s, :);
        TriFace{s} = tempFace(triangulation4orthProj({Vrtx(:, tempFace)}, edgeConON, figON));
    end
end


end

