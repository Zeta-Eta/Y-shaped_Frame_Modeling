function sidFace = getSideFace(botFace, topFace, obj)
% get Side Face
if strcmp(obj, 'Rectangle')
    sidFace = [reshape(botFace(:, 1:end-1)', [], 1), ...
        reshape(botFace(:, 2:end)', [], 1), ...
        reshape(topFace(:, 2:end)', [], 1), ...
        reshape(topFace(:, 1:end-1)', [], 1)];
elseif strcmp(obj, 'Corner')
    sidFace = reshape(topFace', [], 1);
    sidFace = repmat(sidFace, 1, 3) + (1:-1:-1);
    sidFace = [sidFace, reshape(botFace', [], 1)];
elseif strcmp(obj, 'Corner4BluntEdge')
    sidFace = [repmat(reshape(topFace', [], 1), 1, 3) + (1:-1:-1), ...
        repmat(reshape(botFace', [], 1), 1, 3) + (-1:1)];
end

end