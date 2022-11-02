function levelVrtx = getLevel(P, drctVec)
% get Level's Vertices
len = drctVec(1:2, :)\diff(P(1:2, :), 1, 2);
levelVrtx = P(:, 1) + len(1).*drctVec(:, 1);
% or use: 
% levelVrtx = P(:, 2) - len(2).*drctVec(:, 2);

end

