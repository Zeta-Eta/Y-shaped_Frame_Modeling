function R = rotMatrix(theta, k)
% Rotation Matrix
% through an angle 'theta' (counter clock wise)
% about the rotation axis 'k' (a unit vector: k = [k_x; k_y; k_z])

% Cross-product Matrix
K = [  0  , -k(3),  k(2);
      k(3),   0  , -k(1);
     -k(2),  k(1),   0  ];

% Rodrigues' Rotation Formula
R = eye(3) + sin(theta).*K + (1 - cos(theta)).*(K^2);

end

