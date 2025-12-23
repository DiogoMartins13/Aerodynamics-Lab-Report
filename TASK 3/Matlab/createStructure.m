function [ControlPoints, InducedPoints, Normals, InfiniteVortices, Vortices, internalMesh, WingExtremes2Export] = createStructure(config, iBody)

XAngle = config.DihedralAngle(iBody);
XRot = [ 1         0              0      ;
         0    cosd(XAngle)  -sind(XAngle);
         0    sind(XAngle)   cosd(XAngle) ];

YAngle = config.RotationAngle_Y(iBody);
YRot = [ cosd(YAngle)     0        sind(YAngle);
              0           1            0
        -sind(YAngle)     0        cosd(YAngle)];

% Global position of the LE
GlobalPos = [config.LEPosition_X(iBody), config.LEPosition_Y(iBody), config.LEPosition_Z(iBody)];

%% Right Semi-Wing

TipChordVal = config.RootChord(iBody) * config.TaperRatio(iBody);

Local_RootLE = [0, 0, 0];
Local_RootTE = [config.RootChord(iBody), 0, 0];

IsoscelesOffset = (config.RootChord(iBody) - TipChordVal) / 2;

% Sweep angle effect
SweepOffset = config.SemiSpan(iBody) * tand(config.SweepAngle(iBody));

TipLE_X = IsoscelesOffset + SweepOffset;

Local_TipLE  = [TipLE_X, config.SemiSpan(iBody), 0];
Local_TipTE  = [TipLE_X + TipChordVal, config.SemiSpan(iBody), 0];

% Function to rotate and traslate
TransformPoint = @(P) (YRot * XRot * P')' + GlobalPos;

WingExtremes.RootLE = TransformPoint(Local_RootLE);
WingExtremes.RootTE = TransformPoint(Local_RootTE);
WingExtremes.TipLE  = TransformPoint(Local_TipLE);
WingExtremes.TipTE  = TransformPoint(Local_TipTE);

% Saving (optional)
WingExtremes2Export = cell(3,1); 
WingExtremes2Export{1}.LE = WingExtremes.TipLE;
WingExtremes2Export{1}.TE = WingExtremes.TipTE;

%% Mesh Generation

% Spanwise discr.
LEDiscr = linspace(0, 1, config.SemiSpanwiseDiscr(iBody)+1);
RootDiscr = linspace(0, 1, config.ChordwiseDiscr(iBody)+1);

NumChord = config.ChordwiseDiscr(iBody);
NumSpanSemi = config.SemiSpanwiseDiscr(iBody);
NumSpanTot = NumSpanSemi * 2;

internalMesh = cell(NumChord, NumSpanTot);
Vortices = cell(NumChord, NumSpanTot);
ControlPoints = cell(NumChord, NumSpanTot);
InducedPoints = cell(NumChord, NumSpanTot);
Normals = cell(NumChord, NumSpanTot);
InfiniteVortices = cell(NumChord, NumSpanTot);

InfiniteVorticesLength = 50 * config.RootChord(1);

% RX semi_wing
LE_Line = @(s) WingExtremes.RootLE + s * (WingExtremes.TipLE - WingExtremes.RootLE);
TE_Line = @(s) WingExtremes.RootTE + s * (WingExtremes.TipTE - WingExtremes.RootTE);

for i = 1:NumChord
    c_start = RootDiscr(i);
    c_end   = RootDiscr(i+1);
    
    for j = 1:NumSpanSemi
        s_start = LEDiscr(j);
        s_end   = LEDiscr(j+1);
        
        % Interp. before on the span, then on the chord
        P_LE_root = LE_Line(s_start); P_TE_root = TE_Line(s_start);
        P_LE_tip  = LE_Line(s_end);   P_TE_tip  = TE_Line(s_end);
        
        P4 = P_LE_root + c_start * (P_TE_root - P_LE_root); % LE Root
        P3 = P_LE_tip  + c_start * (P_TE_tip  - P_LE_tip);  % LE Tip
        P2 = P_LE_tip  + c_end   * (P_TE_tip  - P_LE_tip);  % TE Tip
        P1 = P_LE_root + c_end   * (P_TE_root - P_LE_root); % TE Root
        
        internalMesh{i, j}.LERoot = P4;
        internalMesh{i, j}.LEtip  = P3;
        internalMesh{i, j}.TEtip  = P2;
        internalMesh{i, j}.TERoot = P1;
        
        % Bound vortex
        Vortices{i, j}.Root = P4 + 0.25 * (P1 - P4);
        Vortices{i, j}.Tip  = P3 + 0.25 * (P2 - P3);
        
        % Collocation points
        CP_Root_Geom = P4 + 0.75 * (P1 - P4);
        CP_Tip_Geom  = P3 + 0.75 * (P2 - P3);
        ControlPoints{i, j}.Coords = (CP_Root_Geom + CP_Tip_Geom) / 2;
        
        % Induced points
        InducedPoints{i, j}.Coords = (Vortices{i, j}.Root + Vortices{i, j}.Tip) / 2;
        
        % Normals
        vecA = P2 - P4; vecB = P3 - P1;
        N = cross(vecA, vecB); 
        Normals{i, j}.Coords = N / norm(N);
        
        % Infinite Vortices
        InfiniteVortices{i, j}.Root.onWing = Vortices{i, j}.Root;
        InfiniteVortices{i, j}.Root.toInfty = Vortices{i, j}.Root + [InfiniteVorticesLength, 0, 0];
        
        InfiniteVortices{i, j}.Tip.onWing = Vortices{i, j}.Tip;
        InfiniteVortices{i, j}.Tip.toInfty = Vortices{i, j}.Tip + [InfiniteVorticesLength, 0, 0];
    end
end

%% Left Semi-Wing
% Symmetry wrt XZ plane (Y -> -Y)
% The RX semi-wing columns are 1:NumSpanSemi
% The LX semi-wing columns will be NumSpanSemi+1 : NumSpanTot

for i = 1:NumChord
    for j = 1:NumSpanSemi
        idx_sym = NumSpanSemi + j;
        MeshOrig = internalMesh{i, j};
        VortOrig = Vortices{i, j};
        CPOrig   = ControlPoints{i, j};
        IndOrig  = InducedPoints{i, j};
        NormOrig = Normals{i, j};
        
        % Mesh
        internalMesh{i, idx_sym}.LERoot = [MeshOrig.LEtip(1),  -MeshOrig.LEtip(2),  MeshOrig.LEtip(3)];
        internalMesh{i, idx_sym}.LEtip  = [MeshOrig.LERoot(1), -MeshOrig.LERoot(2), MeshOrig.LERoot(3)];
        internalMesh{i, idx_sym}.TERoot = [MeshOrig.TEtip(1),  -MeshOrig.TEtip(2),  MeshOrig.TEtip(3)];
        internalMesh{i, idx_sym}.TEtip  = [MeshOrig.TERoot(1), -MeshOrig.TERoot(2), MeshOrig.TERoot(3)];
        
        % Vortices (invert Root/Tip to maintain the correct circulation)
        Vortices{i, idx_sym}.Root = [VortOrig.Tip(1),  -VortOrig.Tip(2),  VortOrig.Tip(3)];
        Vortices{i, idx_sym}.Tip  = [VortOrig.Root(1), -VortOrig.Root(2), VortOrig.Root(3)];
        
        % Control Points e Induced Points
        ControlPoints{i, idx_sym}.Coords = [CPOrig.Coords(1), -CPOrig.Coords(2), CPOrig.Coords(3)];
        InducedPoints{i, idx_sym}.Coords = [IndOrig.Coords(1), -IndOrig.Coords(2), IndOrig.Coords(3)];
        
        % Normals
        Normals{i, idx_sym}.Coords = [NormOrig.Coords(1), -NormOrig.Coords(2), NormOrig.Coords(3)];
        
        % Infinite Vortices
        InfiniteVortices{i, idx_sym}.Root.onWing  = Vortices{i, idx_sym}.Root;
        InfiniteVortices{i, idx_sym}.Root.toInfty = Vortices{i, idx_sym}.Root + [InfiniteVorticesLength, 0, 0];
        InfiniteVortices{i, idx_sym}.Tip.onWing   = Vortices{i, idx_sym}.Tip;
        InfiniteVortices{i, idx_sym}.Tip.toInfty  = Vortices{i, idx_sym}.Tip + [InfiniteVorticesLength, 0, 0];
    end
end

%% Organization
% Now we have RX (Root->Tip), LX (Tip->Root)
% We want LX (Tip->Root), DX (Root->Tip)

Full_Mesh = cell(NumChord, NumSpanTot);
Full_Vort = cell(NumChord, NumSpanTot);
Full_CP   = cell(NumChord, NumSpanTot);
Full_Ind  = cell(NumChord, NumSpanTot);
Full_Norm = cell(NumChord, NumSpanTot);
Full_InfV = cell(NumChord, NumSpanTot);

for i = 1:NumChord
    for j = 1:NumSpanSemi
        idx_old = NumSpanSemi + j;
        idx_new = NumSpanSemi - j + 1;
        
        Full_Mesh{i, idx_new} = internalMesh{i, idx_old};
        Full_Vort{i, idx_new} = Vortices{i, idx_old};
        Full_CP{i, idx_new}   = ControlPoints{i, idx_old};
        Full_Ind{i, idx_new}  = InducedPoints{i, idx_old};
        Full_Norm{i, idx_new} = Normals{i, idx_old};
        Full_InfV{i, idx_new} = InfiniteVortices{i, idx_old};
    end
    
    for j = 1:NumSpanSemi
        idx_old = j;
        idx_new = NumSpanSemi + j;
        
        Full_Mesh{i, idx_new} = internalMesh{i, idx_old};
        Full_Vort{i, idx_new} = Vortices{i, idx_old};
        Full_CP{i, idx_new}   = ControlPoints{i, idx_old};
        Full_Ind{i, idx_new}  = InducedPoints{i, idx_old};
        Full_Norm{i, idx_new} = Normals{i, idx_old};
        Full_InfV{i, idx_new} = InfiniteVortices{i, idx_old};
    end
end

internalMesh = Full_Mesh;
Vortices = Full_Vort;
ControlPoints = Full_CP;
InducedPoints = Full_Ind;
Normals = Full_Norm;
InfiniteVortices = Full_InfV;

WingExtremes2Export{1}.LE = internalMesh{1, end}.LEtip;
WingExtremes2Export{1}.TE = internalMesh{end, end}.TEtip;
WingExtremes2Export{2}.LE = internalMesh{1, NumSpanSemi+1}.LERoot;
WingExtremes2Export{2}.TE = internalMesh{end, NumSpanSemi+1}.TERoot;
WingExtremes2Export{3}.LE = internalMesh{1, 1}.LEtip;
WingExtremes2Export{3}.TE = internalMesh{end, 1}.TEtip;

end
