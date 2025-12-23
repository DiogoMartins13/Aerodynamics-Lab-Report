function [U] = vortexInfluence(TargetPoint, P1_Array, P2_Array)
% Computes the velocity induced by many points over one target point
% TargetPoint: vector 1x3
% P1_Array, P2_Array: matrices Nx3 (list of coordinates)

    % Vectors for all the panels at the same time
    r1 = TargetPoint - P1_Array; 
    r2 = TargetPoint - P2_Array;
    r0 = P2_Array - P1_Array;

    % Cross product along the rows
    CP = cross(r1, r2, 2); 
    CP_sq = sum(CP.^2, 2);

    % Norm of r1 and r2
    r1_norm = sqrt(sum(r1.^2, 2));
    r2_norm = sqrt(sum(r2.^2, 2));

    % Scalar r0
    term1 = r1 ./ r1_norm;
    term2 = r2 ./ r2_norm;
    dot_term = dot(r0, term1 - term2, 2);

    % Manage singularities (if the point is on the vortex)
    toll = 1e-10;
    CP_sq(CP_sq < toll) = toll;

    % Final velocity
    coeff = (1 ./ (4*pi)) .* (dot_term ./ CP_sq);
    U = coeff .* CP; 
end