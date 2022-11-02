function TriFace = getTriFace4TB(FaceID, topON, Vrtx, outFace, holFace, edgeConON, figON)
% get Triangulation Face for Top and Bottom
tempID = find(holFace.ID == FaceID);
tempIDnum = length(tempID);
tempCell = cell(1, 1 + tempIDnum);
origFace = [];
if topON
    tempFace = outFace;
else
    tempFace = flip(outFace, 2);
end
tempCell{1} = Vrtx(:, tempFace);
origFace = [origFace, tempFace];

for i = 1:tempIDnum
    if topON
        tempFace = flip(holFace.top(tempID(i), 1:end-1), 2);
    else
        tempFace = holFace.bot(tempID(i), 1:end-1);
    end
    tempCell{i + 1} = Vrtx(:, tempFace);
    origFace = [origFace, tempFace];
end

TriFace = origFace(triangulation4orthProj(tempCell, edgeConON, figON));

end

