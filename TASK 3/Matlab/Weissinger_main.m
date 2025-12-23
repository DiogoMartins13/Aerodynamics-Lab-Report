close all
clear 
clc

%% Information for the final plots

% To create the final plots where the informations about the two planes are
% overlapped:
% 1. run it using 'Cessna' function;
% 2. comment the clear in line 2, and the 'Cessna' function (lines 22 and 23);
% 3. remove the comment from 'Piper' function;
% 4. run it again.

%% Data

U_Inf_Mag = 67;
alpha_deg = 2; % Angle of attack (having RotationAngle_Y fixed (= 0))
beta = 0;
U_Inf = [cosd(alpha_deg)*cosd(beta), sind(beta), sind(alpha_deg)] .* U_Inf_Mag;
rho = 1.225;

[config] = Cessna172();
AircraftName = 'Cessna';
% [config] = PiperPA28180Cherokee();
% AircraftName = 'Piper';

%% Preliminary computations

% Computing the span
config.SemiSpan = config.Span./2; % half wing 0 to b/2
% Computing the surface per body
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2); 
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle); 
% Computing the tip chord
config.TipChord = config.RootChord .* config.TaperRatio;
% Computing the MAC by definition
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio)); 

%% Create the geometry structure

ControlPoints = cell(config.NBodies, 1);    
InducedPoints = cell(config.NBodies, 1);
Normals = cell(config.NBodies, 1);
InfiniteVortices = cell(config.NBodies, 1);
Vortices = cell(config.NBodies, 1);
internalMesh = cell(config.NBodies, 1);
WingExtremes = cell(config.NBodies, 1);

for iBody = 1:config.NBodies
   [ControlPoints{iBody}, InducedPoints{iBody}, Normals{iBody}, InfiniteVortices{iBody}, Vortices{iBody}, internalMesh{iBody}, WingExtremes{iBody}] = createStructure(config, iBody);
end
  
%% Vectorization & Normal correction

totalPanels = 0;
for i = 1:config.NBodies
    totalPanels = totalPanels + config.ChordwiseDiscr(i) * 2 * config.SemiSpanwiseDiscr(i);
end

% Memory allocation
CP_All = zeros(totalPanels, 3); 
Normals_All = zeros(totalPanels, 3); 
IndP_All = zeros(totalPanels, 3); 

V_InfRoot_P1 = zeros(totalPanels, 3);
V_InfRoot_P2 = zeros(totalPanels, 3);
V_Bound_P1 = zeros(totalPanels, 3);
V_Bound_P2 = zeros(totalPanels, 3);
V_InfTip_P1 = zeros(totalPanels, 3);
V_InfTip_P2 = zeros(totalPanels, 3);

Map_Body = zeros(totalPanels, 1);
Map_Col = zeros(totalPanels, 1); 
Map_Row = zeros(totalPanels, 1); 

idx = 0;

for iBody = 1:config.NBodies
   for C = 1:config.ChordwiseDiscr(iBody)
       for S = 1:2*config.SemiSpanwiseDiscr(iBody)
           idx = idx + 1;
           Map_Body(idx) = iBody;
           Map_Row(idx)  = C;
           Map_Col(idx)  = S;
           CP_All(idx, :)   = ControlPoints{iBody}{C,S}.Coords;
           IndP_All(idx, :) = InducedPoints{iBody}{C,S}.Coords;
           
           % If the normal points downward we change the sign
           current_normal = Normals{iBody}{C,S}.Coords;
           if current_normal(3) < 0
               current_normal = -current_normal;
               Normals{iBody}{C,S}.Coords = current_normal; 
           end
           Normals_All(idx, :) = current_normal;
           
           % Extraction of vortex geometry
           V_InfRoot_P1(idx, :) = InfiniteVortices{iBody}{C,S}.Root.toInfty;
           V_InfRoot_P2(idx, :) = InfiniteVortices{iBody}{C,S}.Root.onWing;
           
           V_Bound_P1(idx, :)   = Vortices{iBody}{C,S}.Root;
           V_Bound_P2(idx, :)   = Vortices{iBody}{C,S}.Tip;
           
           V_InfTip_P1(idx, :)  = InfiniteVortices{iBody}{C,S}.Tip.onWing;
           V_InfTip_P2(idx, :)  = InfiniteVortices{iBody}{C,S}.Tip.toInfty;
       end
   end
end

%% Construction of the matrix

matrixA = zeros(totalPanels, totalPanels);
for i = 1:totalPanels
    Target = CP_All(i, :);       
    Normal = Normals_All(i, :);  
    
    % Influence of all panels (j) on panel (i)
    U_Bound = vortexInfluence(Target, V_Bound_P1, V_Bound_P2);
    U_InfRoot = vortexInfluence(Target, V_InfRoot_P1, V_InfRoot_P2);
    U_InfTip  = vortexInfluence(Target, V_InfTip_P1, V_InfTip_P2);
    U_Total = U_Bound + U_InfRoot + U_InfTip; 
    
    % Projection on the normal 
    matrixA(i, :) = U_Total * Normal'; 
end

%% Known term and solution

knownTerm = - (Normals_All * U_Inf');
Solution = linsolve(matrixA, knownTerm);

Gamma = cell(config.NBodies, 1);

for iBody = 1:config.NBodies
   Gamma{iBody} = zeros(config.ChordwiseDiscr(iBody), config.SemiSpanwiseDiscr(iBody)*2);
end

for k = 1:totalPanels
    b = Map_Body(k);
    r = Map_Row(k);
    c = Map_Col(k);
    Gamma{b}(r, c) = Solution(k);
end

%% Visualization

% PLOT 2D
figure('Name', 'Gamma Distribution');

if config.NBodies == 3
    % Body 1 and 2 MUST have the same config.ChordwiseDiscr
    N_semi_span_2 = config.SemiSpanwiseDiscr(2);
    Gamma_Tip_Left = Gamma{2}(:, 1:N_semi_span_2);
    Gamma_Tip_Right = Gamma{2}(:, N_semi_span_2+1:end);
    Gamma_Root = Gamma{1};
    Gamma_Wing_Total = [Gamma_Tip_Left, Gamma_Root, Gamma_Tip_Right];

    % Wing
    subplot(2,1,1);
    pcolor(Gamma_Wing_Total);
    shading interp; axis tight; colorbar;
    title('Main Wing');
    xlabel('Spanwise Index'); ylabel('Chordwise Index');

    % Tail
    subplot(2,1,2);
    pcolor(Gamma{3});
    shading interp; axis tight; colorbar;
    title('Tail');
    xlabel('Spanwise Index'); ylabel('Chordwise Index');

else
    for i = 1:config.NBodies
        subplot(config.NBodies, 1, i);
        pcolor(Gamma{i}); shading interp; axis tight; colorbar;
        title(['Body ', num2str(i), ' Gamma Distribution']);
        xlabel('Spanwise Index'); ylabel('Chordwise Index');
    end
end


%% Visualization 3D

figure('Name', '3D Geometry');
hold on; grid on; axis equal; 
view([-30, 30]); 
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title([AircraftName, ' wings view and normals direction']);
show_normals = true;  % 'false' to hide the normals
normal_scale = 0.5;
all_gamma = [];

for k=1:config.NBodies 
    all_gamma = [all_gamma; Gamma{k}(:)];
end

clim([min(all_gamma), max(all_gamma)]); 
colormap jet;

for iBody = 1:config.NBodies
    Mesh = internalMesh{iBody};
    [Nr, Nc] = size(Mesh);
    for i = 1:Nr
        for j = 1:Nc
            P1 = Mesh{i,j}.TERoot; P2 = Mesh{i,j}.TEtip;
            P3 = Mesh{i,j}.LEtip;  P4 = Mesh{i,j}.LERoot;
            
            X = [P1(1) P2(1) P3(1) P4(1)];
            Y = [P1(2) P2(2) P3(2) P4(2)];
            Z = [P1(3) P2(3) P3(3) P4(3)];
            
            patch(X, Y, Z, Gamma{iBody}(i,j), 'EdgeColor', 'none', 'FaceAlpha', 0.9);
        end
    end
    
    if show_normals
        % Only the first line of panels along the span
        num_arrows = Nc;
        QX = zeros(1, num_arrows); QY = zeros(1, num_arrows); QZ = zeros(1, num_arrows);
        QU = zeros(1, num_arrows); QV = zeros(1, num_arrows); QW = zeros(1, num_arrows);
        
        for j = 1:Nc
            P1 = Mesh{1,j}.TERoot; P2 = Mesh{1,j}.TEtip;
            P3 = Mesh{1,j}.LEtip;  P4 = Mesh{1,j}.LERoot;
            Center = (P1+P2+P3+P4)/4;
            
            NormVec = Normals{iBody}{1,j}.Coords;
            
            QX(j) = Center(1); QY(j) = Center(2); QZ(j) = Center(3);
            QU(j) = NormVec(1); QV(j) = NormVec(2); QW(j) = NormVec(3);
        end
        
        quiver3(QX, QY, QZ, QU * normal_scale, QV * normal_scale, QW * normal_scale, ...
            0, 'r', 'LineWidth', 1.5, 'AutoScale', 'off');
    end
end

colorbar;


%% Compute the 2D and 3D Lift

L_Total = 0;
Lift_Per_Body = zeros(config.NBodies, 1);

for iBody = 1:config.NBodies
    Nspan = 2 * config.SemiSpanwiseDiscr(iBody);
    dy = config.Span(iBody) / Nspan;

    Gamma_Section = sum(Gamma{iBody}, 1); 
    L_2D_dist = rho * U_Inf_Mag * Gamma_Section * cosd(config.DihedralAngle(iBody));

    L_Body = sum(L_2D_dist * dy);
    L_Total = L_Total + L_Body;
    Lift_Per_Body(iBody) = L_Body;
end

if config.NBodies == 3
    L_MainWing = Lift_Per_Body(1) + Lift_Per_Body(2);
    L_Tail = Lift_Per_Body(3);
    S_ref = config.Surface(1) + config.Surface(2);
    CL_Total = L_Total / (0.5 * rho * U_Inf_Mag^2 * S_ref);

    fprintf('--------------------------------------\n');
    fprintf('Lift:\n');
    fprintf('  > MAIN WING:  %.2f N\n', L_MainWing);
    fprintf('  > TAIL:       %.2f N\n', L_Tail);
    fprintf('--------------------------------------\n');
    fprintf('Total Lift L = %.2f N\n', L_Total);
    fprintf('Total Lift Coefficient C_L = %.4f \n', CL_Total);
    fprintf('--------------------------------------\n');
else
    L_MainWing = Lift_Per_Body(1); 
    L_Tail = Lift_Per_Body(2);
    S_ref = config.Surface(1);
    CL_Total = L_Total / (0.5 * rho * U_Inf_Mag^2 * S_ref);

    fprintf('--------------------------------------\n');
    fprintf('Lift:\n');
    fprintf('  > MAIN WING:  %.2f N\n', L_MainWing);
    fprintf('  > TAIL:       %.2f N\n', L_Tail);
    fprintf('--------------------------------------\n');
    fprintf('Total Lift L = %.2f N\n', L_Total);
    fprintf('Total Lift Coefficient C_L = %.4f \n', CL_Total);
    fprintf('--------------------------------------\n');
end


%% Compute 2D and 3D induced drag 

Di_Total = 0; 
alpha_ind_deg = cell(config.NBodies, 1);
Di_Distribution = cell(config.NBodies, 1);
Drag_Per_Body = zeros(config.NBodies, 1); 

totalPanels = 0;
for i = 1:config.NBodies
    totalPanels = totalPanels + config.ChordwiseDiscr(i) * 2 * config.SemiSpanwiseDiscr(i);
end

All_Wake_Root_P1 = zeros(totalPanels, 3);
All_Wake_Root_P2 = zeros(totalPanels, 3);
All_Wake_Tip_P1 = zeros(totalPanels, 3);
All_Wake_Tip_P2 = zeros(totalPanels, 3);
All_Gamma = zeros(totalPanels, 1);

idx = 0;
for jBody = 1:config.NBodies
    for C = 1:config.ChordwiseDiscr(jBody)
        for S = 1:2*config.SemiSpanwiseDiscr(jBody)
            idx = idx + 1;
            
            % Gamma of this panel
            All_Gamma(idx) = Gamma{jBody}(C, S);
            
            % Wake's coordinates (root)
            All_Wake_Root_P1(idx, :) = InfiniteVortices{jBody}{C,S}.Root.toInfty;
            All_Wake_Root_P2(idx, :) = InfiniteVortices{jBody}{C,S}.Root.onWing;
            
            % Wake's coordinates (tip)
            All_Wake_Tip_P1(idx, :)  = InfiniteVortices{jBody}{C,S}.Tip.onWing;
            All_Wake_Tip_P2(idx, :)  = InfiniteVortices{jBody}{C,S}.Tip.toInfty;
        end
    end
end

% Computing drag on the strip

for iBody = 1:config.NBodies

    Nspan = 2 * config.SemiSpanwiseDiscr(iBody);
    
    Di_strip_dist = zeros(1, Nspan);
    alpha_ind_local = zeros(1, Nspan);

    for S = 1:Nspan
        
        % 2D lift (summing Gamma along the chord) 
        Gamma_Strip_Total = sum(Gamma{iBody}(:, S));
        
        % Induced control point at quarter of the chord
        P_LeadingEdge = internalMesh{iBody}{1, S}.LERoot;   
        P_TrailingEdge = internalMesh{iBody}{end, S}.TERoot; 
        Chord_Vector = P_TrailingEdge - P_LeadingEdge;
        ICP = P_LeadingEdge + 0.25 * Chord_Vector; 
        
        Normal = Normals{iBody}{1,S}.Coords;

        % Induced velocity
        % Root vortices influence
        U_Root = vortexInfluence(ICP, All_Wake_Root_P1, All_Wake_Root_P2);
        
        % Tip vortices influence
        U_Tip = vortexInfluence(ICP, All_Wake_Tip_P1, All_Wake_Tip_P2);
        
        Vind_total = sum(U_Root .* All_Gamma, 1) + sum(U_Tip .* All_Gamma, 1);
        
        w_ind = dot(Vind_total, Normal);
        alpha_i = atan(-w_ind / U_Inf_Mag); 
        
        P_Root_Strip = Vortices{iBody}{1,S}.Root;
        P_Tip_Strip = Vortices{iBody}{1,S}.Tip;
        delta_b = norm(P_Tip_Strip - P_Root_Strip); 
        
        L_Strip_Force = rho * U_Inf_Mag * Gamma_Strip_Total * delta_b * cosd(config.DihedralAngle(iBody));
        
        Di_strip_dist(S) = L_Strip_Force * sin(alpha_i);
        alpha_ind_local(S) = alpha_i;
    end
    
    alpha_ind_deg{iBody} = rad2deg(alpha_ind_local);
    Di_Distribution{iBody} = Di_strip_dist;

    Di_Body = sum(Di_strip_dist); 
    Di_Total = Di_Total + Di_Body;
    Drag_Per_Body(iBody) = Di_Body;
end

if config.NBodies == 3
    Di_MainWing = Drag_Per_Body(1) + Drag_Per_Body(2);
    Di_Tail = Drag_Per_Body(3);
    S_ref = config.Surface(1) + config.Surface(2);
    CDi_Total = Di_Total / (0.5 * rho * U_Inf_Mag^2 * S_ref);

    fprintf('--------------------------------------\n');
    fprintf('Induced Drag:\n');
    fprintf('  > MAIN WING:  %.4f N\n', Di_MainWing);
    fprintf('  > TAIL:       %.4f N\n', Di_Tail);
    if Di_Tail < 0
        fprintf('--------------------------------------\n');
        fprintf('This is not the induced drag of the tail, but the drag induced considering even the effect of the main wing!\n');
    end
    fprintf('--------------------------------------\n');
    fprintf("Total Induced Drag Di = %.6f N\n", Di_Total);
    fprintf("Total Induced Drag Coefficient C_Di = %.6f\n", CDi_Total);
    fprintf('--------------------------------------\n');
else
    Di_MainWing = Drag_Per_Body(1);
    Di_Tail = Drag_Per_Body(2);
    S_ref = config.Surface(1);
    CDi_Total = Di_Total / (0.5 * rho * U_Inf_Mag^2 * S_ref);

    fprintf('--------------------------------------\n');
    fprintf('Induced Drag:\n');
    fprintf('  > MAIN WING:  %.4f N\n', Di_MainWing);
    fprintf('  > TAIL:       %.4f N\n', Di_Tail);
    if Di_Tail < 0
        fprintf('--------------------------------------\n');
        fprintf('This is not the induced drag of the tail, but the drag induced considering even the effect of the main wing!\n');
    end
    fprintf('--------------------------------------\n');
    fprintf("Total Induced Drag Di = %.6f N\n", Di_Total);
    fprintf("Total Induced Drag Coefficient C_Di = %.6f\n", CDi_Total);
    fprintf('--------------------------------------\n');
end


%% Post-processing

if config.NBodies == 3
    WingBodies = [1, 2]; % Cessna
else
    WingBodies = [1];    % Piper
end

Wing_Eta = [];
Wing_Gamma = [];
Wing_Cdi_Local = [];

Min_Y = 1e9; 
Max_Y = -1e9;
Total_Wing_Area_Check = 0; 
q_inf = 0.5 * rho * U_Inf_Mag^2;

for k = 1:length(WingBodies)
    iBody = WingBodies(k);

    % Data check
    if ~exist('Di_Distribution', 'var') || length(Di_Distribution) < iBody
        error('Missing Di_Distribution. Compute first the Drag!');
    end

    Di_Force_Vec = Di_Distribution{iBody}; 
    Gamma_Vec = sum(Gamma{iBody}, 1);

    [~, N_Strips] = size(internalMesh{iBody});

    for s = 1:N_Strips
        Mesh_Strip = internalMesh{iBody}{1, s};

        y_root = (Mesh_Strip.LERoot(2) + Mesh_Strip.TERoot(2))/2;
        y_tip  = (Mesh_Strip.LEtip(2)  + Mesh_Strip.TEtip(2))/2;
        y_center = (y_root + y_tip) / 2;

        % Local chord
        c_root = norm(Mesh_Strip.TERoot - Mesh_Strip.LERoot);
        c_tip  = norm(Mesh_Strip.TEtip  - Mesh_Strip.LEtip);
        c_local = (c_root + c_tip) / 2;

        dy_local = abs(y_tip - y_root);
        Area_Strip = c_local * dy_local;

        % If the area of the strip is close to 0 we don't consider it
        if Area_Strip < 1e-6, continue; end

        Total_Wing_Area_Check = Total_Wing_Area_Check + Area_Strip;

        cdi_val = Di_Force_Vec(s) / (q_inf * Area_Strip);

        % Delete wrong values
        if abs(cdi_val) > 0.5, cdi_val = NaN; end

        Wing_Eta       = [Wing_Eta, y_center];
        Wing_Gamma     = [Wing_Gamma, Gamma_Vec(s)];
        Wing_Cdi_Local = [Wing_Cdi_Local, cdi_val];

        Min_Y = min(Min_Y, min(Mesh_Strip.LERoot(2), Mesh_Strip.LEtip(2)));
        Max_Y = max(Max_Y, max(Mesh_Strip.LERoot(2), Mesh_Strip.LEtip(2)));
    end
end

b_wing = Max_Y - Min_Y;
Wing_Eta = 2 * Wing_Eta / b_wing;

[Wing_Eta, sortIdx] = sort(Wing_Eta);
Wing_Gamma     = Wing_Gamma(sortIdx);
Wing_Cdi_Local = Wing_Cdi_Local(sortIdx);

% Symmetry imposed
if abs(beta) < 1e-5
    Wing_Cdi_Local(isnan(Wing_Cdi_Local)) = 0; 

    N_pts = length(Wing_Cdi_Local);
    Half = floor(N_pts/2);

    Avg_Cdi = (Wing_Cdi_Local(1:Half) + fliplr(Wing_Cdi_Local(end-Half+1:end))) / 2;

    Wing_Cdi_Local(1:Half) = Avg_Cdi;
    Wing_Cdi_Local(end-Half+1:end) = fliplr(Avg_Cdi);

    Avg_Gam = (Wing_Gamma(1:Half) + fliplr(Wing_Gamma(end-Half+1:end))) / 2;
    Wing_Gamma(1:Half) = Avg_Gam;
    Wing_Gamma(end-Half+1:end) = fliplr(Avg_Gam);
end

if ~exist('AircraftName', 'var'), AircraftName = ['Config ' num2str(config.NBodies)]; end

% Elliptic gamma (only for the wing)
L_Wing_Only = 0;
for k=1:length(WingBodies), L_Wing_Only = L_Wing_Only + Lift_Per_Body(WingBodies(k)); end

Gamma_0 = (4 * L_Wing_Only) / (pi * rho * U_Inf_Mag * b_wing);

Arg_Sqrt = 1 - Wing_Eta.^2;
Arg_Sqrt(Arg_Sqrt<0) = 0;
Gamma_Ell = Gamma_0 * sqrt(Arg_Sqrt);

Current_Data.Name = AircraftName;
Current_Data.Eta = Wing_Eta;
Current_Data.Gamma_W = Wing_Gamma;
Current_Data.Gamma_E = Gamma_Ell;
Current_Data.cdi_local = Wing_Cdi_Local;

if ~exist('ComparisonStore', 'var'), ComparisonStore = {}; end
Found = false;
for k = 1:length(ComparisonStore)
    if strcmp(ComparisonStore{k}.Name, AircraftName)
        ComparisonStore{k} = Current_Data; Found = true; break;
    end
end
if ~Found, ComparisonStore{end+1} = Current_Data; end

%% Plots

Colors = {'b', 'r', 'g', 'm', 'k'};

% Figure 1: Circulation
figure('Name', 'Gamma comparison', 'Color', 'w'); hold on; grid on;
title('Circulation distribution', 'FontSize', 14);
xlabel('Normalized span: \eta = 2y/b'); ylabel('Circulation: \Gamma [m^2/s]');

for k = 1:length(ComparisonStore)
    c_idx = mod(k-1, length(Colors)) + 1;
    Data = ComparisonStore{k};
    plot(Data.Eta, Data.Gamma_W, [Colors{c_idx} '-'], 'LineWidth', 2, 'DisplayName', [Data.Name ' (Sim)']);
    plot(Data.Eta, Data.Gamma_E, [Colors{c_idx} '--'], 'LineWidth', 1.5, 'DisplayName', [Data.Name ' (Ell)']);
end
legend('Location', 'South');
xlim([-1.05 1.05]);

% Figure 2: CDi
figure('Name', 'Induced Drag comparison', 'Color', 'w'); hold on; grid on;
title('C_{D_i} distribution', 'FontSize', 14);
xlabel('Normalized span: \eta = 2y/b'); ylabel('C_{D_i}');

ymax_plot = 0;
for k = 1:length(ComparisonStore)
    c_idx = mod(k-1, length(Colors)) + 1;
    Data = ComparisonStore{k};
    plot(Data.Eta, Data.cdi_local, [Colors{c_idx} '-'], 'LineWidth', 2, 'DisplayName', Data.Name);

    valid_vals = Data.cdi_local(~isnan(Data.cdi_local));
    if ~isempty(valid_vals)
        ymax_plot = max(ymax_plot, max(valid_vals)); 
    end
end

if ymax_plot > 0, ylim([0, ymax_plot * 1.2]); end
xlim([-1.05 1.05]);
legend('Location', 'North');


%% Polar curve

Total_Panels = sum(config.ChordwiseDiscr .* 2 .* config.SemiSpanwiseDiscr);
AIC_Drag = zeros(Total_Panels, Total_Panels);
Global_ICP = zeros(Total_Panels, 3);
idx = 0;

for iBody = 1:config.NBodies
    for C = 1:config.ChordwiseDiscr(iBody)
        for S = 1:2*config.SemiSpanwiseDiscr(iBody)
            idx = idx + 1;
            P_Root = Vortices{iBody}{C,S}.Root;
            P_Tip  = Vortices{iBody}{C,S}.Tip;
            Global_ICP(idx, :) = (P_Root + P_Tip) / 2;
        end
    end
end

for i = 1:Total_Panels
    Target = Global_ICP(i, :);
    Norm   = Normals_All(i, :);
    U_Root = vortexInfluence(Target, V_InfRoot_P1, V_InfRoot_P2);
    U_Tip  = vortexInfluence(Target, V_InfTip_P1, V_InfTip_P2);
    AIC_Drag(i, :) = (U_Root + U_Tip) * Norm';
end

Alpha_Range = -4:2:14;
N_Points = length(Alpha_Range);
CL_Polar = zeros(1, N_Points);
CDi_Polar = zeros(1, N_Points);
S_ref_Total = sum(config.Surface(WingBodies)); 

U_Mag_Base = norm(U_Inf);

for k = 1:N_Points
    a_deg = Alpha_Range(k);
    U_Curr = [cosd(a_deg) 0 sind(a_deg)] * U_Mag_Base;

    RHS = - (Normals_All * U_Curr');
    G_Vec = linsolve(matrixA, RHS);
    W_Vec = AIC_Drag * G_Vec;

    L_It = 0; D_It = 0; idx_p = 0;

    for iBody = 1:config.NBodies
        Dihed = config.DihedralAngle(iBody);
        for C = 1:config.ChordwiseDiscr(iBody)
            for S = 1:2*config.SemiSpanwiseDiscr(iBody)
                idx_p = idx_p + 1;
                Gam = G_Vec(idx_p);
                w_i = W_Vec(idx_p);

                P_R = Vortices{iBody}{C,S}.Root;
                P_T = Vortices{iBody}{C,S}.Tip;
                dl  = norm(P_T - P_R);

                al_i = atan(-w_i / U_Mag_Base);
                dL = rho * U_Mag_Base * Gam * dl * cosd(Dihed);
                dDi = dL * sin(al_i);

                L_It = L_It + dL;
                D_It = D_It + dDi;
            end
        end
    end
    CL_Polar(k) = L_It / (0.5 * rho * U_Mag_Base^2 * S_ref_Total);
    CDi_Polar(k) = D_It / (0.5 * rho * U_Mag_Base^2 * S_ref_Total);
end

Current_Polar.Name = AircraftName;
Current_Polar.Alpha = Alpha_Range;
Current_Polar.CL = CL_Polar;
Current_Polar.CDi = CDi_Polar;

if ~exist('PolarStore', 'var'), PolarStore = {}; end
Found = false;
for k = 1:length(PolarStore)
    if strcmp(PolarStore{k}.Name, AircraftName)
        PolarStore{k} = Current_Polar; Found = true; break;
    end
end
if ~Found, PolarStore{end+1} = Current_Polar; end

figure('Name', 'Polar Curve', 'Color', 'w'); 
subplot(1,2,1); hold on; grid on; title('C_L - \alpha'); xlabel('\alpha [deg]'); ylabel('C_L');
subplot(1,2,2); hold on; grid on; title('Lift-Drag Polar'); xlabel('C_{D_i}'); ylabel('C_L');

for k = 1:length(PolarStore)
    c_idx = mod(k-1, length(Colors)) + 1;
    Dat = PolarStore{k};
    subplot(1,2,1); plot(Dat.Alpha, Dat.CL, [Colors{c_idx} '-o'], 'LineWidth', 2, 'DisplayName', Dat.Name);
    subplot(1,2,2); plot(Dat.CDi, Dat.CL, [Colors{c_idx} '-o'], 'LineWidth', 2, 'DisplayName', Dat.Name);
end
legend('Location', 'SouthEast');

