function [U] = uVortex(center, extreme_1, extreme_2, L2G_TransfMatrix, G2L_TransfMatrix)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:    center: coordinates where the velocity is computed
%           extreme_1 & _2: panel extrema
%           L2G_TransfMatrix: 2x2 matrix, local -> global coordinates
%           G2L_TransfMatrix: 2x2 matrix, global -> local coordinates
% Output:   U : x and y velocity components.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Global to local coordinates
center = G2L_TransfMatrix * center;
extreme_1 = G2L_TransfMatrix * extreme_1;
extreme_2 = G2L_TransfMatrix * extreme_2;

%% Local u and v
r1 = center - extreme_1;
theta_1 = atan2(r1(2), r1(1));
r2 = center - extreme_2;
theta_2 = atan2(r2(2), r2(1));

% self induction:
if (abs(theta_1)<10^(-8) && abs(theta_2)>3); theta_1=0; theta_2=pi; end
if (abs(theta_2)<10^(-8) && abs(theta_1)>3); theta_2=0; theta_1=-pi; end

% velocity component
u = theta_2 - theta_1;
u = u / (2*pi);
v = (0.5/pi) * log(norm(r2)/norm(r1));


%% Local to global coordinates

U = L2G_TransfMatrix * [u;v];

if abs(U(1)) < 10^(-12)
    U(1) = 0;
end
if abs(U(2)) < 10^(-12)
    U(2) = 0;
end


