%% description
%%This is for data analysis worksheet question Q6:
% you have learnt from our first-day class how to fit a linear curve to
% data using Excel. Please refer to FITdemo.xls for details.

% Here, we introduce a matlab function polyfit to do the same work.
% Update by QY 2021

%% load data
clear all; close all;
% 
data = [ 50 22.35
25 13.82
10 6.05 
5 2.78
1 0.49
0.5 0.20
0.1 0.033];

A = data(:, 2); % data column 2: absorbance
C = data(:, 1); % data column 1: concentration

%%  plot data
figure(1)
plot(C, A, 'o')
xlabel('concentration (mg/ml)')
ylabel('absorbance')
%% fit data to a line
% plot absorbance vs concentration

figure(2)
plot(C, A, 'o')
xlabel('concentration (mg/ml)')
ylabel('absorbance')

[P S] = polyfit(C(3:end), A(3:end), 1);  %I am not fitting to the higher concentration points
%  P = polyfit(X,Y,N), here N = 1 for linear fit

hold on
plot(C, polyval(P,C), 'r-')  % plot the fit result as a red line
hold off

xlabel('concentration (mg/ml)') %label my axes
ylabel('absorbance')

%% plot and fit concentration vs. absorbance

figure(3)
plot(A, C, 'o')

[P S] = polyfit(A(3:end), C(3:end), 1);

hold on
plot(A, polyval(P,A), 'r-')
hold off
ylabel('concentration (mg/ml)')
xlabel('absorbance')

x = .75; % at an absorbance of x = 0.75, what is the estimate concentration C1
[C1,deltaC1] = polyval(P,x,S);  %% this is one way of calculating errors
% Delta C1 runs error propogation

disp(['C1 = ' num2str(C1) '+/-' num2str(deltaC1)]) %display the calculated values on the command line


ExtCoeff = 1/P(1)* 66776; % molar extinction coefficient is calculated as slope/pathlength with unit [Abs/mg/mL]/[cm]
% molecular weight of BSA (66,776g/mol)
disp(['Molar Extinction Coefficient = ' num2str(ExtCoeff)])


