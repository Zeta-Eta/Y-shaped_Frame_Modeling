function P = getCorner(P0, v1, v2, len)
% get Corner
cos = v1'*v2;
v = v2 - v1.*cos; % could also use: v = cross(cross(v1, v2), v1)
v = v./norm(v, 2);

tan = sqrt((1 - cos)./(1 + cos)); % Half-angle formula
P = P0 + v.*tan.*len;

end

