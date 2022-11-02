function GM = GMof3Dobj(Vrtx, Face)
%% Geometric Measurements of a 3D object
% Faces should be divided into triangles and the numbering of them needs to
% follow the right-hand or left-hand rule. 
% ↑ But these two rules cannot be mixed. 

A = 0; % Area
V = 0; % Volume
CofA = zeros(3, 1); % Centroid of Area
CofV = zeros(3, 1); % Centroid of Volume

for n = 1:size(Face, 1)
    v1 = Vrtx(:, Face(n, 1));
    v2 = Vrtx(:, Face(n, 2));
    v3 = Vrtx(:, Face(n, 3));
    
    tempA = norm(cross(v2 - v1, v3 - v1), 2)./2;
    tempV = v1'*cross(v2, v3)./6;
    A = A + tempA;
    V = V + tempV;
    CofA = CofA + mean([v1, v2, v3], 2).*tempA;
    CofV = CofV + mean([v1, v2, v3, zeros(3, 1)], 2).*tempV;
end
CofA = CofA./A;
CofV = CofV./V;

GM.A = A;
GM.V = V;
GM.CofA = CofA;
GM.CofV = CofV;


end

