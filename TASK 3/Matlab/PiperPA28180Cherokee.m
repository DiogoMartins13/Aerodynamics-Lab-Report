function [config] = PiperPA28180Cherokee ()

% Piper PA-28-180 Cherokee

config.NBodies = 2;
config.RootChord = [1.6; 0.762]; % [m]
config.DihedralAngle = [7; 0]; % [°]
config.SweepAngle = [0; 0]; % [°]
config.TaperRatio = [1; 1]; % =TipChord/RootChord
config.AspectRatio = [5.63; 4.1];
config.Span = [9.144; 3.048];
config.LEPosition_X = [0; 4.14];
config.LEPosition_Y = [0; 0];
config.LEPosition_Z = [0; 0.5];
Tilt_angle = -4;
config.RotationAngle_X = [0; 0];
config.RotationAngle_Y = [0; Tilt_angle];
config.RotationAngle_Z = [0; 0];

% Discretization options
config.SemiSpanwiseDiscr = [30; 20];
config.ChordwiseDiscr = [20; 20];

end