function P_2D = orthProj(P_3D)
% orthographic projection 3D to 2D
% 1st Point is convex
normVec = cross(P_3D(:, 1) - P_3D(:, end), P_3D(:, 2) - P_3D(:, 1));
normVec = normVec./norm(normVec, 2);

axisVec = [0, 0, 1]';

if all(normVec == axisVec)
    R = eye(3);
elseif all(normVec == -axisVec)
    R = rotx(180);
elseif normVec(3) >= 0
    rotAxisVec = normVec + axisVec;
    rotAxisVec = rotAxisVec./norm(rotAxisVec, 2);
    
    R = rotMatrix(pi, rotAxisVec);
    R = rotz(180)*R;
elseif normVec(3) < 0
    rotAxisVec = normVec - axisVec;
    rotAxisVec = rotAxisVec./norm(rotAxisVec, 2);
    
    R = rotMatrix(pi, rotAxisVec);
    R = rotx(180)*rotz(180)*R;
end

P_2D = R(1:2, :)*P_3D;

end

