function [config] = Cessna172 ()

% Cessna 172 Skyhawk

config.NBodies = 3;
% Body 1 rectangular part of main wing, Body 2 trapezohidal part of main wing, Body 3 tail
config.RootChord = [1.625; 1.625; 1.4]; % [m]
config.DihedralAngle = [0; 0; 0]; % [°]
config.SweepAngle = [0; 0; 0]; % [°]
config.TaperRatio = [1; 0.672; 0.8/1.4]; % =TipChord/RootChord
config.AspectRatio = [(5.36/1.625); ((1.625+1.092)*5.64/2); (3.4^2/3.74)];
config.Span = [5.36; 5.64; 3.4];
config.LEPosition_X = [0; 0; 4.25];
config.LEPosition_Y = [0; 5.36/2; 0];
config.LEPosition_Z = [0; 0; -0.65];
Tilt_angle = -6.5;
config.RotationAngle_X = [0; 0; 0];
config.RotationAngle_Y = [0; 0; Tilt_angle];
config.RotationAngle_Z = [0; 0; 0];

% Discretization options
config.SemiSpanwiseDiscr = [15; 15; 20];
config.ChordwiseDiscr = [20; 20; 20];

end