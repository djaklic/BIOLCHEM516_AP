clc;
clear;

time = [0 46 69 92 115 138 161 184 207 230];
samples = table2array(load("biolchem516_AP2_partC.mat").x);

figure;
subplot(1,2,1);
scatter(time, samples');
title("Activity Assay Raw Data")
xlabel("time (s)")
ylabel("absorbance")

[P_1 S_1, fit_1, error_1, slope_ci_1, inter_ci_1] = ...
    calc_params(time,samples(:,1));
[P_2 S_2, fit_2, error_2, slope_ci_2, inter_ci_2] = ...
    calc_params(time,samples(:,2));
[P_3 S_3, fit_3, error_3, slope_ci_3, inter_ci_3] = ...
    calc_params(time,samples(:,3));
[P_4 S_4, fit_4, error_4, slope_ci_4, inter_ci_4] = ...
    calc_params(time,samples(:,4));
[P_5_15 S_5_15, fit_5_15, error_5_15, slope_ci_5_15, inter_ci_5_15] = ...
    calc_params(time,samples(:,5));
[P_5_16 S_5_16, fit_5_16, error_5_16, slope_ci_5_16, inter_ci_5_16] = ...
    calc_params(time,samples(:,6));
[P_5_17 S_5_17, fit_5_17, error_5_17, slope_ci_5_17, inter_ci_5_17] = ...
    calc_params(time,samples(:,7));

subplot(1,2,2);
hold on
errorbar(time, fit_1, error_1);
errorbar(time, fit_2, error_2);
errorbar(time, fit_3, error_3);
errorbar(time, fit_4, error_4);
errorbar(time, fit_5_15, error_5_15);
errorbar(time, fit_5_16, error_5_16);
errorbar(time, fit_5_17, error_5_17);
title("Activity Assay Fit and Error")
xlabel("time (s)")
legend("stage 1","stage 2", ...
    "stage 3","stage 4", ...
    "stage 5-15","stage 5-16","stage 5-17");
hold off

abs_velocities = [P_1(1) P_2(1) P_3(1) ...
    P_4(1) P_5_15(1) P_5_16(1) P_5_17(1)];
abs_velocities_errors = [slope_ci_1(2) - slope_ci_1(1), ...
    slope_ci_2(2) - slope_ci_2(1), ...
    slope_ci_3(2) - slope_ci_3(1), ...
    slope_ci_4(2) - slope_ci_4(1), ...
    slope_ci_5_15(2) - slope_ci_5_15(1), ...
    slope_ci_5_16(2) - slope_ci_5_16(1), ...
    slope_ci_5_17(2) - slope_ci_5_17(1)];
conc_velocities = abs_velocities/20000; 

%error = error of slope + error of ext. coeff
syms A_vel EC
C_vel = A_vel/EC;
C_vel_error_1 = PropError( ...
    C_vel,[A_vel EC],[abs_velocities(1) 20000],[abs_velocities_errors(1) 200]);
C_vel_error_2 = PropError( ...
    C_vel,[A_vel EC],[abs_velocities(2) 20000],[abs_velocities_errors(2) 200]);
C_vel_error_3 = PropError( ...
    C_vel,[A_vel EC],[abs_velocities(3) 20000],[abs_velocities_errors(3) 200]);
C_vel_error_4 = PropError( ...
    C_vel,[A_vel EC],[abs_velocities(4) 20000],[abs_velocities_errors(4) 200]);
C_vel_error_5_15 = PropError( ...
    C_vel,[A_vel EC],[abs_velocities(5) 20000],[abs_velocities_errors(5) 200]);
C_vel_error_5_16 = PropError( ...
    C_vel,[A_vel EC],[abs_velocities(6) 20000],[abs_velocities_errors(6) 200]);
C_vel_error_5_17 = PropError( ...
    C_vel,[A_vel EC],[abs_velocities(7) 20000],[abs_velocities_errors(7) 200]);

function [P S fit error slope_ci inter_ci] = calc_params(time, sample)
    [P, S] = polyfit(time, sample',1);
    [fit, error] = polyval(P,time,S);
    [ci] = polyparci(P,S);
    slope_ci = ci(:,1);
    inter_ci = ci(:,2);
end