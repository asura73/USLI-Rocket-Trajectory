function T = rotmat(axis, angle)

%angle in radians


% return 3d matrix to rotate a column vector
% by angle. if axis not betweeb [1 3]
% return identity matrix

switch axis
    case 1
        T = [1, 0, 0; ...
            0, cos(angle), sin(angle); ...
            0, -sin(angle), cos(angle)];
    case 2
        T = [cos(angle), 0, -sin(angle); ...
            0, 1, 0; ...
            sin(angle), 0, cos(angle)];
    case 3
        T = [cos(angle), sin(angle), 0; ...
            -sin(angle), cos(angle), 0; ...
            0, 0, 1];
    otherwise
        T = eye(3);
end


end