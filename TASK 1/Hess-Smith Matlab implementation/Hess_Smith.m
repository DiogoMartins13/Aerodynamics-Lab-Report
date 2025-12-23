%% Hess Smith method
% Diogo Martins
% Gabriele Pagnoni
% Riccardo Rossetti
% Tommaso Rossi
% Lorenzo Rota

%% Clear data
clc
close all
clear

%% Get Xfoil CP Data
addpath mat_functions

data1 = readtable('cp_NACA0008_a1_v2.dat');
xfoil_a1 = table2array(data1);
cp_xfoil_a1 = xfoil_a1(:,3);
Cp_top_xfoil_a1 = xfoil_a1(1:51,3);
Cp_bottom_xfoil_a1 = xfoil_a1(53:102,3);
x_top_xfoil_a1 = xfoil_a1(1:51,1);
x_bottom_xfoil_a1 = xfoil_a1(53:102,1);

data2 = readtable('cp_NACA0008_a2_v2.dat');
xfoil_a2 = table2array(data2);
cp_xfoil_a2 = xfoil_a2(:,2);
Cp_top_xfoil_a2 = xfoil_a2(1:52,2);
Cp_bottom_xfoil_a2 = xfoil_a2(54:102,2);
x_top_xfoil_a2 = xfoil_a2(1:52,1);
x_bottom_xfoil_a2 = xfoil_a2(54:102,1);

%% Input

U_inf = 1;      % Far-field velocity [m/s]
AoA = 1;        % Angle of attack
U_inf_x = U_inf * cos(deg2rad(AoA));
U_inf_y = U_inf * sin(deg2rad(AoA));
U_inf = [U_inf_x; U_inf_y];

Chord = 1;
NPanels = 101;

LE_X_Position = 0;
LE_Y_Position = 0;

fprintf('Angle of Attack = %.4f \n \n', AoA);
if AoA ~= 1 && AoA ~= 2
    fprintf('Comparision between XFOIL data and Hess-Smith only for AoA = 1 and AoA = 2! \n \n')
end

%% Create profile (with xfoil)

[x,y]=CreateProfile('0008',NPanels,Chord);

geo.x=x;
geo.y=y;

figure
plot(x,y,'o-')
axis equal
title ('NACA 0008')

%% Create discretization & initialization

[centers,normals,tangent,extrema_1,extrema_2,alpha,lengths,L2G_TransfMatrix,G2L_TransfMatrix] = CreatePanels(geo);
        
NCols = sum(NPanels) + 1;
NRows = NCols;
A = zeros(NRows,NCols);     % system coefficients
B = zeros(NRows,1);         % known terms

%% Fill matrix A

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

%% Create a_v, c_s, and c_v vectors

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

%% Create vector B, the known terms in the system 

for j = 1:NPanels
    local_normal = normals(j, :)';
    B(j) = - dot(U_inf, local_normal);
end

first_tangent = tangent(1, :)';
last_tangent = tangent(end, :)';
B(sum(NPanels) + 1) = - dot(U_inf, (first_tangent + last_tangent));

%% Solve the linear system

solution = linsolve(A,B);

%% Compute velocity and Cp

q = solution(1:NPanels);        % source intensity
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
    V_i = vel_vect(i, :);
    t_i = tangent(i, :);
    V_t = dot(V_i, t_i);
    vel_vect_tan(i) = V_t;
    Cp(i) = 1 - (vel_vect_tan(i) ./ U_inf_abs).^2;
end

%% Cp distribution

[~, Cp_max_idx] = max(Cp);
Cp_bottom = flipud(Cp(1:Cp_max_idx));
Cp_top = Cp((Cp_max_idx+1):end);
x_bottom = flipud(centers((1:Cp_max_idx), 1));
x_top = centers(((Cp_max_idx+1):end), 1);
figure
hold on
plot(x_top / Chord, -Cp_top, 'b-', 'DisplayName', 'Suction side')
plot(x_bottom / Chord, -Cp_bottom, 'r-', 'DisplayName', 'Pressure side')
xlabel('x/c')
ylabel('-C_p')
title('Pressure Coefficient (Cp) - Hess Smith')
grid on
legend('show') 

%% CL and Cm

n_U_inf = [-sin(deg2rad(AoA)); cos(deg2rad(AoA))];
Cl_i = zeros(NPanels, 1);

for i = 1:NPanels
   n_i = normals(i, :)';
   Cl_i(i) = -Cp(i) * (lengths(i) / Chord) * normals(i,2);
end

Cl = sum(Cl_i);
fprintf('Lift Coefficient (Cl): %.4f\n', Cl);

% Kutta Joukowski

Circ = gamma * sum(lengths(1:NPanels));
Cl_KJ = 2 * Circ./U_inf_abs;

fprintf('Lift Coefficient with Kutta Joukowski theorem (Cl_KJ): %.4f\n \n', Cl_KJ);

Cm_i = zeros(NPanels, 1);
Cm_i_ac = zeros(NPanels, 1);
x_ac = Chord/4;

for i = 1:NPanels
   Cm_i(i) = Cp(i) * (lengths(i) / (Chord^2)) * normals(i,2) * centers(i,1);
   Cm_i_ac(i) = Cp(i) * (lengths(i) / (Chord^2)) * normals(i,2) * (centers(i,1)-x_ac);

   Cm_xfoil_a1(i) = cp_xfoil_a1(i) * (lengths(i) / (Chord^2)) * normals(i,2) * centers(i,1);
   Cm_xfoil_a1_ac(i) = cp_xfoil_a1(i) * (lengths(i) / (Chord^2)) * normals(i,2) * (centers(i,1)-x_ac);

   Cm_xfoil_a2(i) = cp_xfoil_a2(i) * (lengths(i) / (Chord^2)) * normals(i,2) * centers(i,1);
   Cm_xfoil_a2_ac(i) = cp_xfoil_a2(i) * (lengths(i) / (Chord^2)) * normals(i,2) * (centers(i,1)-x_ac);
end

Cm_LE = sum(Cm_i);
fprintf('Leading Edge Moment Coefficient (Cm_LE): %.8f\n', Cm_LE);
Cm_AC = sum(Cm_i_ac);
fprintf('Aerodynamic Center Moment Coefficient (Cm_AC): %.8f\n \n', Cm_AC);

Cm_LE_xfoil_a1 = sum(Cm_xfoil_a1);
fprintf('Leading Edge Moment Coefficient (Cm_LE) [XFOIL AoA = 1]: %.8f\n', Cm_LE);
Cm_AC_xfoil_a1 = sum(Cm_xfoil_a1_ac);
fprintf('Aerodynamic Center Moment Coefficient (Cm_AC) [XFOIL AoA = 1]: %.8f\n \n', Cm_AC);

Cm_LE_xfoil_a2 = sum(Cm_xfoil_a2);
fprintf('Leading Edge Moment Coefficient (Cm_LE) [XFOIL AoA = 2]: %.8f\n', Cm_LE);
Cm_AC_xfoil_a2 = sum(Cm_xfoil_a2_ac);
fprintf('Aerodynamic Center Moment Coefficient (Cm_AC) [XFOIL AoA = 2]: %.8f\n', Cm_AC);

%% Plots

figure()
hold on
if AoA == 1
   plot(x_top_xfoil_a1 / Chord, (-Cp_top)-(-Cp_top_xfoil_a1), 'g--')
    plot(x_bottom_xfoil_a1 / Chord, (-Cp_bottom)-(-Cp_bottom_xfoil_a1), 'y--')
end
if AoA == 2
   plot(x_top_xfoil_a2 / Chord, (-Cp_top)-(-Cp_top_xfoil_a2), 'g--')
    plot(x_bottom_xfoil_a2 / Chord, (-Cp_bottom)-(-Cp_bottom_xfoil_a2), 'y--')
end
xlabel('x/c')
ylabel('-C_p')
title('Difference between Hess-Smith and Xfoil')
legend('Top Hess-Smith - Xfoil','Bottom Hess-Smith - Xfoil')
hold off

figure()
plot(x, 2*y, 'LineWidth', 1.5, 'Color', 'k');
hold on
plot(x_top / Chord, -Cp_top, 'b-')
plot(x_bottom / Chord, -Cp_bottom, 'r-')
if AoA == 1
   plot(x_top_xfoil_a1 / Chord, -Cp_top_xfoil_a1, 'g--')
    plot(x_bottom_xfoil_a1 / Chord, -Cp_bottom_xfoil_a1,'LineStyle', '--', 'Color', '#FFA500', 'LineWidth', 1)
end
if AoA == 2
   plot(x_top_xfoil_a2 / Chord, -Cp_top_xfoil_a2, 'g--')
    plot(x_bottom_xfoil_a2 / Chord, -Cp_bottom_xfoil_a2,'LineStyle', '--', 'Color', '#FFA500', 'LineWidth', 1)
end
xlabel('x/c')
ylabel('-C_p')
title('NACA 0008')
legend ('Airfoil', 'Suction Side Hess-Smith','Pressure Side Hess-Smith', 'Suction Side Xfoil', 'Pressure Side Xfoil')
grid on

% annotation('textbox', [0.73, 0.75, 0.25, 0.15], ...
%     'String', sprintf('\\alpha = %.0f°\nC_L = %.4f\nC_m_A_C = %.4f', AoA, Cl, Cm_AC), ...
%     'FitBoxToText', 'on', ...
%     'BackgroundColor', 'white', ...
%     'EdgeColor', 'black', ...
%     'Margin', 5, ...
%     'FontSize', 12);

hold off