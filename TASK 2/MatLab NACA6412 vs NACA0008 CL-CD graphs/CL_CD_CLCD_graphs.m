%% Task 2

clc
close all
clear


%% CL and CD data

% Reynolds = 10^6 and N = 9, for -4º, -3º, -2º
CL_NACA0008 = [-0.5214, -0.3905, -0.2637, -0.1334, 0, 0.1334, 0.2637, 0.3905, 0.5214]
CD_NACA0008 = [0.00857, 0.00710, 0.00552, 0.00437, 0.00388, 0.00437, 0.00552, 0.00710, 0.00857]
CL_NACA6412 = [0.1464, 0.2224, 0.2998, 0.3827, 0.4789, 0.5878, 0.7459, 0.9515, 1.0387]
CD_NACA6412 = [0.00854, 0.00819, 0.00806, 0.00808, 0.00817, 0.00829, 0.00778, 0.00773, 0.00819]

alpha = linspace (-4,4,9)


%% Plots

figure()
subplot(1,3,1)
hold on
plot(alpha, CL_NACA0008)
plot(alpha, CL_NACA6412)
xlabel('Alpha (º)')
ylabel('C_l')
title('C_l Comparison')
legend ('NACA0008','NACA6412')
grid on
hold off

subplot(1,3,2)
hold on
plot(alpha, CD_NACA0008)
plot(alpha, CD_NACA6412)
xlabel('Alpha (º)')
ylabel('C_d')
title('C_d Comparison')
grid on
hold off
legend ('NACA0008', 'NACA6412')

subplot(1,3,3)
hold on
plot(alpha, CL_NACA0008./CD_NACA0008)
plot(alpha, CL_NACA6412./CD_NACA6412)
title('C_l/C_d Comparison')
xlabel('Alpha (º)')
ylabel('C_l/C_d')
grid on
hold off
legend ('NACA0008', 'NACA6412')


