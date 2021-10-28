clc;
clear;

sc = [0.237 0.291 0.271 0.177 0.184 0.279 ...
    0.364 0.33 0.37 0.439 0.78 1.328];
concentrations = [0 0.01 0.025 0.05 0.1 ...
    0.15 0.2 0.3 0.4 0.5 1 2]; %mg/ml
samples = table2array(load("biolchem516_AP2_partB.mat").x);

samples_mean = mean(samples);
samples_sd = std(samples);

[P S] = polyfit(concentrations, sc, 1);
[fit_sc, error_sc] = polyval(P,concentrations,S);

[P_inv S_inv] = polyfit(sc, concentrations, 1);
[fit_sample, error_sample] = polyval(P_inv,samples_mean,S_inv);


figure;
hold on
scatter(concentrations, sc, 'b.')
plot(concentrations, fit_sc, 'g-')
errorbar(fit_sample, samples_mean, samples_sd, 'c.')
errorbar(fit_sample, samples_mean, error_sample, 'horizontal' ,'r.')
scatter(fit_sample, samples_mean, 'r*')
title("BCA Assay Standard Curve and Sample Predictions")
xlabel("concentration (mg/ml)")
ylabel("absorbance")
legend("sc raw data","sc linear fit", ...
    "sample abs. error","sample conc. error","sample conc. prediction");
xlim([-0.25 2.25])
ylim([0 1.5])
hold off