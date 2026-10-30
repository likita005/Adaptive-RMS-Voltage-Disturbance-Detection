%% TESTING: NOISE ROBUSTNESS AND HARMONIC DISTORTION
% Run main.m first, then run this file.
% This script uses variables created by main.m.
%
% The code below was reconstructed from the supplied testing PDF.
% Review and execute it in MATLAB before relying on its results.

%% =========================================================
% PART 1: NOISE ROBUSTNESS TESTING
%% =========================================================

%% NOISE ROBUSTNESS TESTING
% Run main.m first. This script uses variables created by main.m:
% v, t, L, N_fixed, N_fast, N_detect, omega, change_threshold,
% refractory_samples, hold_samples, Vrms_ideal, Vrms_nominal,
% true_condition.
%
% This file was reconstructed from the supplied PDF. Review and run
% in MATLAB before using it for final results.

fprintf('\n========================================\n');

fprintf('          NOISE ROBUSTNESS TESTING\n');

fprintf('========================================\n');

rng(42);     % Reproducible noise

snr_levels = [Inf 30 20 10];
num_tests = length(snr_levels);

RMSE_fixed_noise = zeros(1,num_tests);
RMSE_adaptive_noise = zeros(1,num_tests);

MAE_fixed_noise = zeros(1,num_tests);
MAE_adaptive_noise = zeros(1,num_tests);

Accuracy_fixed_noise = zeros(1,num_tests);
Accuracy_adaptive_noise = zeros(1,num_tests);

signal_power = mean(v.^2);

for k = 1:num_tests

%% 1. ADD NOISE

if isinf(snr_levels(k))
v_test = v;
else
noise_power = signal_power / ...
(10^(snr_levels(k)/10));

noise = sqrt(noise_power)*randn(size(v));

v_test = v + noise;
end

%% 2. FIXED 200-SAMPLE RMS

fixed_test = nan(1,L);

for n = N_fixed:L
fixed_test(n) = sqrt( ...
mean(v_test(n-N_fixed+1:n).^2));
end

%% 3. FAST 100-SAMPLE RMS

fast_test = nan(1,L);

for n = N_fast:L
fast_test(n) = sqrt( ...
mean(v_test(n-N_fast+1:n).^2));
end

%% 4. 50-SAMPLE CHANGE DETECTOR

detector_test = nan(1,L);

for n = N_detect:L

idx = n-N_detect+1:n;

X = [sin(omega*t(idx))', ...
cos(omega*t(idx))'];

coefficients = X\v_test(idx)';

detector_test(n) = norm(coefficients)/sqrt(2);

end

%% 5. ADAPTIVE WINDOW SELECTION

adaptive_test = nan(1,L);

last_trigger = -inf;
fast_until = 0;

for n = N_fixed:L

difference = abs( ...
detector_test(n)-fixed_test(n));

if difference > change_threshold && ...
(n-last_trigger) > refractory_samples

last_trigger = n;
fast_until = max(fast_until, ...

n+hold_samples-1);

end

if n <= fast_until
adaptive_test(n) = fast_test(n);
else
adaptive_test(n) = fixed_test(n);
end

end

%% 6. COMMON VALID SAMPLES

valid_test = ~isnan(fixed_test) & ...
~isnan(adaptive_test);

%% 7. RMSE AND MAE

error_fixed = fixed_test(valid_test) - ...
Vrms_ideal(valid_test);

error_adaptive = adaptive_test(valid_test) - ...
Vrms_ideal(valid_test);

RMSE_fixed_noise(k) = sqrt(mean(error_fixed.^2));

RMSE_adaptive_noise(k) = sqrt(mean(error_adaptive.^2));

MAE_fixed_noise(k) = mean(abs(error_fixed));

MAE_adaptive_noise(k) = mean(abs(error_adaptive));

%% 8. CLASSIFICATION

fixed_pu_test = fixed_test/Vrms_nominal;

adaptive_pu_test = adaptive_test/Vrms_nominal;

fixed_class_test = nan(1,L);
adaptive_class_test = nan(1,L);

fixed_class_test(valid_test) = 0;
adaptive_class_test(valid_test) = 0;

fixed_class_test(valid_test & ...
fixed_pu_test < 0.1) = 3;

fixed_class_test(valid_test & ...
fixed_pu_test >= 0.1 & fixed_pu_test < 0.9) = 1;

fixed_class_test(valid_test & ...
fixed_pu_test > 1.1) = 2;

adaptive_class_test(valid_test & ...
adaptive_pu_test < 0.1) = 3;

adaptive_class_test(valid_test & ...
adaptive_pu_test >= 0.1 & adaptive_pu_test < 0.9) = 1;

adaptive_class_test(valid_test & ...
adaptive_pu_test > 1.1) = 2;

Accuracy_fixed_noise(k) = 100*mean( ...
fixed_class_test(valid_test) == ...
true_condition(valid_test));

Accuracy_adaptive_noise(k) = 100*mean( ...
adaptive_class_test(valid_test) == ...
true_condition(valid_test));

%% DISPLAY RESULTS

fprintf('\nSNR: %g dB\n',snr_levels(k));

fprintf('Fixed RMSE    : %.4f V\n', ...
RMSE_fixed_noise(k));

fprintf('Adaptive RMSE : %.4f V\n', ...
RMSE_adaptive_noise(k));

fprintf('Fixed MAE     : %.4f V\n', ...
MAE_fixed_noise(k));

fprintf('Adaptive MAE : %.4f V\n', ...
MAE_adaptive_noise(k));

fprintf('Fixed Accuracy    : %.2f %%\n', ...
Accuracy_fixed_noise(k));

fprintf('Adaptive Accuracy : %.2f %%\n', ...
Accuracy_adaptive_noise(k));

end

Fixed RMSE    : 31.4354 V
Adaptive RMSE : 23.3276 V
Fixed MAE     : 8.3472 V
Adaptive MAE : 4.1723 V

Fixed RMSE     : 31.4307 V
Adaptive RMSE : 23.3380 V
Fixed MAE      : 8.9556 V
Adaptive MAE : 4.8658 V
Fixed RMSE     : 31.7081 V
Adaptive RMSE : 23.7555 V
Fixed MAE      : 10.8859 V
Adaptive MAE : 7.0119 V
Fixed RMSE     : 37.5396 V
Adaptive RMSE : 31.9360 V
Fixed MAE      : 22.4725 V
Adaptive MAE : 19.7479 V

%% 9. RESULTS TABLE

NoiseResults = table( ...
snr_levels', ...
RMSE_fixed_noise', ...
RMSE_adaptive_noise', ...
MAE_fixed_noise', ...
MAE_adaptive_noise', ...
Accuracy_fixed_noise', ...
Accuracy_adaptive_noise', ...
'VariableNames', { ...
'SNR_dB', ...
'Fixed_RMSE', ...
'Adaptive_RMSE', ...
'Fixed_MAE', ...
'Adaptive_MAE', ...
'Fixed_Accuracy', ...
'Adaptive_Accuracy'});

fprintf('\n========================================\n');

fprintf('          FINAL NOISE RESULTS\n');

fprintf('========================================\n');

disp(NoiseResults);

Inf       31.435        23.328        8.3472        4.1723           92.44            96.378
30       31.431        23.338        8.9556        4.8658           92.45            96.398
20       31.708        23.756        10.886        7.0119          91.134            94.388
10        37.54        31.936        22.473        19.748          84.512            87.083

%% 10. RMSE COMPARISON PLOT

figure('Name','Noise RMSE Comparison','Color','w');

plot(1:num_tests,RMSE_fixed_noise,'-ob', ...
'LineWidth',1.5);

hold on;

plot(1:num_tests,RMSE_adaptive_noise,'-sr', ...
'LineWidth',1.5);

grid on;

xticks(1:num_tests);
xticklabels({'Clean','30 dB','20 dB','10 dB'});

xlabel('Noise Condition');
ylabel('RMSE (V)');

title('RMS Estimation Error Under Noise');

legend('Fixed RMS','Adaptive RMS','Location','best');

%% 11. ACCURACY COMPARISON PLOT

figure('Name','Noise Accuracy Comparison','Color','w');

plot(1:num_tests,Accuracy_fixed_noise,'-ob', ...
'LineWidth',1.5);

hold on;

plot(1:num_tests,Accuracy_adaptive_noise,'-sr', ...
'LineWidth',1.5);

grid on;

xticks(1:num_tests);
xticklabels({'Clean','30 dB','20 dB','10 dB'});

xlabel('Noise Condition');
ylabel('Classification Accuracy (%)');

title('Classification Accuracy Under Noise');

legend('Fixed RMS','Adaptive RMS','Location','best');

fprintf('\nNoise testing completed successfully.\n');


%% =========================================================
% PART 2: HARMONIC DISTORTION TESTING
%% =========================================================

%% HARMONIC DISTORTION TESTING
% Run main.m first. This script uses variables created by main.m:
% t, L, N_fixed, N_fast, N_detect, omega, Vpeak, Vrms_nominal,
% sag, swell, interruption, change_threshold, refractory_samples,
% hold_samples.
%
% This file was reconstructed from the supplied PDF. Review and run
% in MATLAB before using it for final results.

fprintf('\n========================================\n');

fprintf('       HARMONIC DISTORTION TESTING\n');

fprintf('========================================\n');

%% 1. HARMONIC PARAMETERS

h3 = 0.15;
h5 = 0.10;

h7 = 0.05;

%% 2. GENERATE DISTORTED VOLTAGE

% Fundamental and harmonic components

harmonic_wave = sin(omega*t) ...
+ h3*sin(3*omega*t) ...
+ h5*sin(5*omega*t) ...
+ h7*sin(7*omega*t);

% Apply voltage disturbance amplitudes

amplitude = ones(1,L);

amplitude(sag) = 0.50;
amplitude(swell) = 1.30;
amplitude(interruption) = 0.02;

v_harmonic = Vpeak*amplitude.*harmonic_wave;

%% 3. TOTAL NOMINAL RMS REFERENCE

% RMS includes fundamental and all harmonic components

Vrms_harmonic_nominal = Vrms_nominal*sqrt( ...
1+h3^2+h5^2+h7^2);

Vrms_harmonic_ideal = ...
Vrms_harmonic_nominal*amplitude;

%% 4. FIXED 200-SAMPLE RMS

harmonic_fixed = nan(1,L);

for n = N_fixed:L

segment = v_harmonic(n-N_fixed+1:n);

harmonic_fixed(n) = sqrt(mean(segment.^2));

end

%% 5. 100-SAMPLE RMS

harmonic_fast = nan(1,L);

for n = N_fast:L

segment = v_harmonic(n-N_fast+1:n);

harmonic_fast(n) = sqrt(mean(segment.^2));

end

%% 6. 50-SAMPLE FUNDAMENTAL CHANGE DETECTOR

harmonic_detector = nan(1,L);

for n = N_detect:L

idx = n-N_detect+1:n;

X = [sin(omega*t(idx))', ...
cos(omega*t(idx))'];

coefficients = X\v_harmonic(idx)';

harmonic_detector(n) = ...
norm(coefficients)/sqrt(2);

end

%% 7. IMPROVED ADAPTIVE WINDOW SELECTION
% Compare fundamental estimates one complete cycle apart.
% This reduces false triggers caused by periodic harmonic variation.

harmonic_adaptive = nan(1,L);
harmonic_window = nan(1,L);

last_trigger = -inf;
fast_until = 0;

for n = N_fixed:L

if n > N_fixed + N_detect

previous_estimate = harmonic_detector(n-N_fixed);

current_estimate = harmonic_detector(n);

difference = abs(current_estimate-previous_estimate);

% Detect actual fundamental magnitude changes
if difference > change_threshold && ...
(n-last_trigger) > refractory_samples

last_trigger = n;

fast_until = max(fast_until, ...
n+hold_samples-1);

end

end

% Select estimation window
if n <= fast_until

harmonic_adaptive(n) = harmonic_fast(n);

harmonic_window(n) = N_fast;

else

harmonic_adaptive(n) = harmonic_fixed(n);

harmonic_window(n) = N_fixed;

end

end

%% 8. PER-UNIT VOLTAGE

harmonic_fixed_pu = ...
harmonic_fixed/Vrms_harmonic_nominal;

harmonic_adaptive_pu = ...
harmonic_adaptive/Vrms_harmonic_nominal;

%% 9. CLASSIFICATION

harmonic_fixed_class = nan(1,L);
harmonic_adaptive_class = nan(1,L);

valid_harmonic = ...
~isnan(harmonic_fixed) & ...
~isnan(harmonic_adaptive);

harmonic_fixed_class(valid_harmonic) = 0;
harmonic_adaptive_class(valid_harmonic) = 0;

harmonic_fixed_class(valid_harmonic & ...
harmonic_fixed_pu < 0.1) = 3;

harmonic_fixed_class(valid_harmonic & ...
harmonic_fixed_pu >= 0.1 & ...

harmonic_fixed_pu < 0.9) = 1;

harmonic_fixed_class(valid_harmonic & ...
harmonic_fixed_pu > 1.1) = 2;

harmonic_adaptive_class(valid_harmonic & ...
harmonic_adaptive_pu < 0.1) = 3;

harmonic_adaptive_class(valid_harmonic & ...
harmonic_adaptive_pu >= 0.1 & ...
harmonic_adaptive_pu < 0.9) = 1;

harmonic_adaptive_class(valid_harmonic & ...
harmonic_adaptive_pu > 1.1) = 2;

%% 10. PERFORMANCE EVALUATION

fixed_harmonic_error = ...
harmonic_fixed(valid_harmonic) - ...
Vrms_harmonic_ideal(valid_harmonic);

adaptive_harmonic_error = ...
harmonic_adaptive(valid_harmonic) - ...
Vrms_harmonic_ideal(valid_harmonic);

RMSE_fixed_harmonic = ...
sqrt(mean(fixed_harmonic_error.^2));

RMSE_adaptive_harmonic = ...
sqrt(mean(adaptive_harmonic_error.^2));

MAE_fixed_harmonic = ...
mean(abs(fixed_harmonic_error));

MAE_adaptive_harmonic = ...
mean(abs(adaptive_harmonic_error));

accuracy_fixed_harmonic = 100*mean( ...
harmonic_fixed_class(valid_harmonic) == ...
true_condition(valid_harmonic));

accuracy_adaptive_harmonic = 100*mean( ...
harmonic_adaptive_class(valid_harmonic) == ...
true_condition(valid_harmonic));

%% 11. DISPLAY RESULTS

fprintf('\nHARMONIC COMPONENTS\n');

fprintf('3rd harmonic : %.1f %%\n',h3*100);

3rd harmonic : 15.0 %

fprintf('5th harmonic : %.1f %%\n',h5*100);

5th harmonic : 10.0 %

fprintf('7th harmonic : %.1f %%\n',h7*100);

7th harmonic : 5.0 %

fprintf('\nTotal nominal RMS : %.3f V\n', ...
Vrms_harmonic_nominal);

Total nominal RMS : 233.990 V

fprintf('\n========================================\n');

fprintf('HARMONIC TEST PERFORMANCE\n');

fprintf('========================================\n');

fprintf('\nFixed RMS:\n');


fprintf('RMSE           : %.4f V\n',RMSE_fixed_harmonic);


fprintf('MAE            : %.4f V\n',MAE_fixed_harmonic);


fprintf('Accuracy : %.2f %%\n',accuracy_fixed_harmonic);


fprintf('\nAdaptive RMS:\n');


fprintf('RMSE           : %.4f V\n',RMSE_adaptive_harmonic);


fprintf('MAE          : %.4f V\n',MAE_adaptive_harmonic);


fprintf('Accuracy : %.2f %%\n',accuracy_adaptive_harmonic);


%% 12. DISTORTED INPUT VOLTAGE

figure('Name','Harmonic Input Voltage','Color','w');

plot(t,v_harmonic,'b');

grid on;

xlabel('Time (s)');
ylabel('Voltage (V)');

title('Voltage Waveform with Harmonic Distortion');

%% 13. FIXED VS ADAPTIVE RMS

figure('Name','Harmonic RMS Comparison','Color','w');

plot(t,harmonic_fixed,'b','LineWidth',1.2);

hold on;

plot(t,harmonic_adaptive,'r','LineWidth',1.2);

yline(Vrms_harmonic_nominal,'--k','Nominal RMS');

grid on;

xlabel('Time (s)');
ylabel('RMS Voltage (V)');

title('Fixed vs Adaptive RMS Under Harmonic Distortion');

legend('Fixed RMS','Adaptive RMS','Nominal RMS');

%% 14. PER-UNIT RMS

figure('Name','Harmonic Per Unit RMS','Color','w');

plot(t,harmonic_adaptive_pu,'r','LineWidth',1.2);

hold on;

yline(0.1,'--k','Interruption');
yline(0.9,'--m','Sag');

yline(1.1,'--g','Swell');

grid on;

xlabel('Time (s)');
ylabel('Voltage (pu)');

title('Adaptive RMS Under Harmonic Distortion');

ylim([0 1.5]);

%% 15. WINDOW SELECTION

figure('Name','Harmonic Window Selection','Color','w');

stairs(t,harmonic_window,'LineWidth',1.2);

grid on;

xlabel('Time (s)');
ylabel('Window Length (samples)');

title('Adaptive Window Selection Under Harmonics');

yticks([100 200]);

ylim([0 220]);

fprintf('\nHarmonic distortion experiment completed.\n');
