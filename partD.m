clc;
clear;

time = [0 20 30 40 50 60 70 80 90 ...
    100 110 120 130 140 150 160 170 180 190 200];
samples = table2array(load("biolchem516_AP2_partD.mat").x);
E = 20000; %extinction coeff M^-1*cm^-1
E_mM = E/1000; %mM^-1*cm^-1

figure;
subplot(1,2,1)
scatter(time, samples');
title("Activity Assay Raw Data")
xlabel("time (s)")
ylabel("absorbance")

[P_10, S_10, fit_10, error_10, slope_ci_10, inter_ci_10] = ...
    calc_params(time, samples(:,1)');
[P_30, S_30, fit_30, error_30, slope_ci_30, inter_ci_30] = ...
    calc_params(time, samples(:,2)');
[P_50, S_50, fit_50, error_50, slope_ci_50, inter_ci_50] = ...
    calc_params(time, samples(:,3)');
[P_100, S_100, fit_100, error_100, slope_ci_100, inter_ci_100] = ...
    calc_params(time, samples(:,4)');
[P_200, S_200, fit_200, error_200, slope_ci_200, inter_ci_200] = ...
    calc_params(time, samples(:,5)');

subplot(1,2,2);
hold on
errorbar(time, fit_10, error_10);
errorbar(time, fit_30, error_30);
errorbar(time, fit_50, error_50);
errorbar(time, fit_100, error_100);
errorbar(time, fit_200, error_200);
title("Activity Assay Fit and Error")
xlabel("time (s)")
legend("10 uM","30 uM","50 uM","100 uM","200 uM");
hold off

concentrations = [0.01 0.03 0.05 0.1 0.2]; %mM
velocities = [P_10(1) P_30(1) P_50(1) P_100(1) P_200(1)]/E_mM; %mM/sec

syms abs e
vel = abs/e;
vel_error = PropError(vel,[abs e],[2.8e-3 200],[1e-04 2]);

figure;
scatter(concentrations, velocities);

concentrations_inverse = 1./concentrations(2:end);
velocities_inverse = 1./velocities(2:end);
[P_MM , S_MM, fit_MM, error_MM, slope_ci_MM, inter_ci_MM] = ...
    calc_params(concentrations_inverse, velocities_inverse);
figure;
hold on
scatter(concentrations_inverse, velocities_inverse);
plot(concentrations_inverse, fit_MM);
hold off

slope_error = slope_ci_MM(2)-P_MM(1);
slope_error_2 = P_MM(1)-slope_ci_MM(1);
intercept_error = inter_ci_MM(2)-P_MM(2);
intercept_error_2 = P_MM(2) - inter_ci_MM(1);

Vmax = 1/P_MM(2);
Vmax_error = intercept_error/(P_MM(2)^2);
km = Vmax*P_MM(1);
km_error = Vmax*slope_error + P_MM(1)*Vmax_error;

% error = slope error + y-intercept error
% 1/V = (Km/Vmax)*(1/S) + 1/Vmax
%syms m x b
%y = m*x+b;
%y_error = PropError(y,[m x b],[212159 20 426403],[0.001 0]);

%syms L g
%T = 2*pi*sqrt(L/g);
%a = PropError(T,[L g],[10 9.81],[0.001 0]);

function [P S fit error slope_ci inter_ci] = calc_params(time, sample)
    [P, S] = polyfit(time, sample,1);
    [fit, error] = polyval(P,time,S);
    [ci] = polyparci(P,S);
    slope_ci = ci(:,1);
    inter_ci = ci(:,2);
end