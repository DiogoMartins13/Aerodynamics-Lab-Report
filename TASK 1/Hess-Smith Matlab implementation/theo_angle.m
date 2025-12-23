%% Hess Smith method

clc
close all
clear

addpath mat_functions
MAX=[];
%% Input
AoAr=-3:0.5:3;
U_inf5 = 1;      % Far-field velocity [m/s]
for AoA=-3:0.5:3
U_inf_x = U_inf5 * cos(deg2rad(AoA));
U_inf_y = U_inf5 * sin(deg2rad(AoA));
U_inf=[];
U_inf = [U_inf_x; U_inf_y];

Chord = 1;
NPanels = 101;

LE_X_Position = 0;
LE_Y_Position = 0;

% Create profile (with xfoil)

[x,y]=createProfile('6412',NPanels,Chord);

geo.x=x;
geo.y=y;

% figure
% plot(x,y,'o-')
% axis equal

% Create discretization & initialization

[centers,normals,tangent,extrema_1,extrema_2,alpha,lengths,L2G_TransfMatrix,G2L_TransfMatrix] = CreatePanels(geo);
        
NCols = sum(NPanels) + 1;
NRows = NCols;
A = zeros(NRows,NCols);     % system coefficients
B = zeros(NRows,1);         % known terms

% Fill A

for i = 1:NPanels
    local_center = centers(i, :)';
    local_normal = normals(i, :)';

    for j = 1:NPanels
        local_extreme_1 = extrema_1(j, :)';
        local_extreme_2 = extrema_2(j, :)';

        local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
        local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));

        A(i, j) = dot(uSource(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), local_normal);

        A(i, sum(NPanels)+1) = A(i, sum(NPanels)+1) + dot(uVortex(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), local_normal);
    end

end

% Create a_v, c_s, and c_v vectors

first_centers = centers(1, :)';
first_tangent = tangent(1, :)';

last_centers = centers(end, :)';
last_tangent = tangent(end, :)';

last_a = 0;
for j = 1:NPanels
    local_extreme_1 = extrema_1(j, :)';
    local_extreme_2 = extrema_2(j, :)';
    local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
    local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(uSource(first_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), first_tangent);
    last_a = last_a + dot(uVortex(first_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), first_tangent);

    a = a + dot(uSource(last_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), last_tangent);
    last_a = last_a + dot(uVortex(last_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), last_tangent);

    A(sum(NPanels) + 1, j) = a;
end

A(sum(NPanels) + 1, sum(NPanels) + 1) = last_a;

% Create B, the known terms in the system 

for j = 1:NPanels
    local_normal = normals(j, :)';
    B(j) = - dot(U_inf, local_normal);
end

first_tangent = tangent(1, :)';
last_tangent = tangent(end, :)';
B(sum(NPanels) + 1) = - dot(U_inf, (first_tangent + last_tangent));

% Solve the linear system

solution = linsolve(A,B);

% Compute velocity and Cp

q = solution(1:NPanels);  % source intensity
gamma = solution(NPanels + 1);  % vortex intensity
vel_vect = zeros(NPanels, 2);
Cp = zeros(NPanels,1);
U_inf_abs = sqrt(U_inf(1)^2 + U_inf(2)^2); 

for i = 1:NPanels
    local_center = centers(i, :)';
    V_i = U_inf; 
    for j = 1:NPanels
        local_extreme_1 = extrema_1(j, :)';
        local_extreme_2 = extrema_2(j, :)';
        local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
        local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));
        
        u_s = uSource(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix);
        V_i = V_i + q(j) .* u_s;
        
        u_v = uVortex(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix);
        V_i = V_i + gamma .* u_v;
    end
    
    vel_vect(i, 1) = V_i(1); % Vx
    vel_vect(i, 2) = V_i(2); % Vy
end

vel_vect_tan = zeros(NPanels, 1);

for i = 1:NPanels
    V_i = vel_vect(i, :)';
    t_i = tangent(i, :)';
    V_t = dot(V_i, t_i);
    vel_vect_tan(i) = V_t;
    Cp(i) = 1 - (vel_vect_tan(i) ./ U_inf_abs).^2;
end


% Cp distribution

[~, le_idx] = min(geo.x);  % leading edge index
idx_low_plot = 1:le_idx; 
idx_up_plot = (le_idx + 1):NPanels; 
x_low = flipud(centers(idx_low_plot, 1));
cp_lower = flipud(Cp(idx_low_plot));
x_up = centers(idx_up_plot, 1);
cp_upper = Cp(idx_up_plot);
% figure
% hold on
% plot(x_up / Chord, -cp_upper, 'b-', 'DisplayName', 'Suction side')
% plot(x_low / Chord, -cp_lower, 'r-', 'DisplayName', 'Pressure side')
% xlabel('x/c')
% ylabel('-C_p')
% title('Pressure Coefficient (Cp)')
% grid on
% legend('show') 


% CL and Cm

n_U_inf = [-sin(deg2rad(AoA)); cos(deg2rad(AoA))];
Cl_contributions = zeros(NPanels, 1);

for i = 1:NPanels
    n_i = normals(i, :)';
    dot_n_L = dot(n_i, n_U_inf);
    Cl_contributions(i) = -Cp(i) * (lengths(i) / Chord) * dot_n_L;
end

Cl = sum(Cl_contributions);
% fprintf('Lift Coefficient (Cl): %.4f\n', Cl);

Cm_contributions = zeros(NPanels, 1);

for i = 1:NPanels
    x_i = centers(i, 1);
    y_i = centers(i, 2);
    n_xi = normals(i, 1);
    n_yi = normals(i, 2);
    cross_term = x_i * n_yi - y_i * n_xi;
    Cm_contributions(i) = Cp(i) * (lengths(i) / (Chord^2)) * cross_term;
end

Cm = sum(Cm_contributions);
% fprintf('Leading Edge Moment Coefficient (Cm_LE): %.4f\n', Cm);
MAX=[MAX,max(-Cp)];
end
%%
[min,index]=min(MAX)
AoAr(index)