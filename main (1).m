%% ADAPTIVE RMS-BASED VOLTAGE DISTURBANCE DETECTION
% Reconstructed from the supplied Adaptive_RMS_Project (2).pdf.
% Simulates sag, swell and interruption; compares fixed and adaptive RMS.
% Please review and run in MATLAB to validate the reconstructed source.

clc;
clear;
close all;

%% =========================================================
% ADAPTIVE RMS-BASED VOLTAGE DISTURBANCE DETECTION
% Detection of Voltage Sag, Swell and Interruption
%
% Fixed RMS window       : 200 samples
% Adaptive fast window   : 100 samples
% Change detection       : 50-sample LS estimator
%% =========================================================

%% 1. SIMULATION PARAMETERS

f = 50;                         % Supply frequency
fs = 10000;                     % Sampling frequency
T = 1;                          % Simulation duration

Vrms_nominal = 230;        % Nominal RMS voltage
Vpeak = sqrt(2)*Vrms_nominal;

t = 0:1/fs:T-1/fs;
L = length(t);

omega = 2*pi*f;

N_fixed = 200;                  % One-cycle window
N_fast = 100;                   % Half-cycle window
N_detect = 50;                  % Change detection window

fprintf('\n========================================\n');

fprintf(' ADAPTIVE RMS VOLTAGE DISTURBANCE MODEL\n');

fprintf('========================================\n');

%% 2. GENERATE VOLTAGE SIGNAL

v = Vpeak*sin(omega*t);

%% 3. INTRODUCE VOLTAGE DISTURBANCES

sag_start = 0.20;

sag_end = 0.35;

swell_start = 0.50;
swell_end = 0.65;

int_start = 0.75;
int_end = 0.85;

sag = (t >= sag_start) & (t < sag_end);
swell = (t >= swell_start) & (t < swell_end);
interruption = (t >= int_start) & (t < int_end);

v(sag) = 0.50*Vpeak*sin(omega*t(sag));

v(swell) = 1.30*Vpeak*sin(omega*t(swell));

v(interruption) = 0.02*Vpeak*sin(omega*t(interruption));

%% 4. IDEAL RMS REFERENCE

Vrms_ideal = Vrms_nominal*ones(1,L);

Vrms_ideal(sag) = 0.50*Vrms_nominal;
Vrms_ideal(swell) = 1.30*Vrms_nominal;
Vrms_ideal(interruption) = 0.02*Vrms_nominal;

%% 5. FIXED 200-SAMPLE RMS ESTIMATION

Vrms_fixed = nan(1,L);

for n = N_fixed:L

segment = v(n-N_fixed+1:n);

Vrms_fixed(n) = sqrt(mean(segment.^2));

end

%% 6. 100-SAMPLE RMS ESTIMATION

Vrms_fast = nan(1,L);

for n = N_fast:L

segment = v(n-N_fast+1:n);

Vrms_fast(n) = sqrt(mean(segment.^2));

end

%% 7. 50-SAMPLE LEAST-SQUARES CHANGE DETECTOR
% This estimator is used ONLY for change detection.
% It is not used as the final total RMS output.

Vrms_detector = nan(1,L);

for n = N_detect:L

idx = n-N_detect+1:n;

X = [sin(omega*t(idx))', ...
cos(omega*t(idx))'];

coefficients = X\v(idx)';

Vrms_detector(n) = norm(coefficients)/sqrt(2);

end

%% 8. MULTI-LEVEL ADAPTIVE WINDOW SELECTION

Vrms_adaptive = nan(1,L);

selected_window = nan(1,L);

% Initial parameters
change_threshold = 0.05*Vrms_nominal;

hold_samples = N_fast;

refractory_samples = 50;

last_trigger = -inf;

fast_until = 0;

for n = N_fixed:L

% Difference between short estimator and full-cycle RMS
difference = abs(Vrms_detector(n)-Vrms_fixed(n));

% Detect significant voltage change
if difference > change_threshold && ...
(n-last_trigger) > refractory_samples

last_trigger = n;

fast_until = max(fast_until,n+hold_samples-1);

end

% Select the RMS estimation window
if n <= fast_until

Vrms_adaptive(n) = Vrms_fast(n);

selected_window(n) = N_fast;

else

Vrms_adaptive(n) = Vrms_fixed(n);

selected_window(n) = N_fixed;

end

end

%% 9. PER-UNIT RMS VOLTAGE

Vpu_fixed = Vrms_fixed/Vrms_nominal;

Vpu_adaptive = Vrms_adaptive/Vrms_nominal;

%% 10. DISTURBANCE CLASSIFICATION
% 0 = Normal
% 1 = Sag
% 2 = Swell
% 3 = Interruption

condition_fixed = nan(1,L);
condition_adaptive = nan(1,L);

valid = ~isnan(Vrms_adaptive) & ~isnan(Vrms_fixed);

% Fixed-window classification
condition_fixed(valid) = 0;

condition_fixed(valid & Vpu_fixed < 0.1) = 3;

condition_fixed(valid & ...
Vpu_fixed >= 0.1 & Vpu_fixed < 0.9) = 1;

condition_fixed(valid & ...
Vpu_fixed > 1.1) = 2;

% Adaptive classification
condition_adaptive(valid) = 0;

condition_adaptive(valid & Vpu_adaptive < 0.1) = 3;

condition_adaptive(valid & ...
Vpu_adaptive >= 0.1 & Vpu_adaptive < 0.9) = 1;

condition_adaptive(valid & ...
Vpu_adaptive > 1.1) = 2;

%% 11. TRUE DISTURBANCE LABELS

true_condition = zeros(1,L);

true_condition(sag) = 1;
true_condition(swell) = 2;
true_condition(interruption) = 3;

%% 12. PERFORMANCE EVALUATION

% Use identical valid samples for both estimators

fixed_error = Vrms_fixed(valid)-Vrms_ideal(valid);

adaptive_error = Vrms_adaptive(valid)-Vrms_ideal(valid);

RMSE_fixed = sqrt(mean(fixed_error.^2));

RMSE_adaptive = sqrt(mean(adaptive_error.^2));

MAE_fixed = mean(abs(fixed_error));

MAE_adaptive = mean(abs(adaptive_error));

accuracy_fixed = 100*mean( ...
condition_fixed(valid)==true_condition(valid));

accuracy_adaptive = 100*mean( ...
condition_adaptive(valid)==true_condition(valid));

%% 13. DETECTION AND RECOVERY DELAY

event_names = {'Sag','Swell','Interruption'};

event_starts = [sag_start swell_start int_start];

event_ends = [sag_end swell_end int_end];

event_codes = [1 2 3];

detection_delay = nan(1,3);
recovery_delay = nan(1,3);

for k = 1:3

start_index = find(t >= event_starts(k),1);

detected = find( ...
condition_adaptive(start_index:end)==event_codes(k),1);

if ~isempty(detected)

detection_index = start_index+detected-1;

detection_delay(k) = ...
(t(detection_index)-event_starts(k))*1000;

end

recovery_index = find(t >= event_ends(k),1);

recovered = find( ...
condition_adaptive(recovery_index:end)==0,1);

if ~isempty(recovered)

actual_index = recovery_index+recovered-1;

recovery_delay(k) = ...
(t(actual_index)-event_ends(k))*1000;

end

end

%% 14. DISPLAY SIMULATION PARAMETERS

fprintf('\nSIMULATION PARAMETERS\n');

fprintf('Supply frequency            : %.1f Hz\n',f);

Supply frequency        : 50.0 Hz

fprintf('Sampling frequency          : %.0f Hz\n',fs);

Sampling frequency      : 10000 Hz

fprintf('Nominal RMS voltage         : %.1f V\n',Vrms_nominal);

Nominal RMS voltage     : 230.0 V

fprintf('Fixed RMS window              : %d samples\n',N_fixed);

Fixed RMS window       : 200 samples

fprintf('Fast RMS window               : %d samples\n',N_fast);

Fast RMS window        : 100 samples

fprintf('Change detector               : %d samples\n',N_detect);

Change detector        : 50 samples

%% 15. DISPLAY RMS PERFORMANCE

fprintf('\n========================================\n');

fprintf('RMS ESTIMATION PERFORMANCE\n');

fprintf('========================================\n');

fprintf('\nFixed RMS:\n');


fprintf('RMSE : %.4f V\n',RMSE_fixed);


fprintf('MAE       : %.4f V\n',MAE_fixed);


fprintf('\nAdaptive RMS:\n');


fprintf('RMSE : %.4f V\n',RMSE_adaptive);


fprintf('MAE       : %.4f V\n',MAE_adaptive);


%% 16. DISPLAY CLASSIFICATION PERFORMANCE

fprintf('\n========================================\n');

fprintf('CLASSIFICATION PERFORMANCE\n');

fprintf('========================================\n');

fprintf('Fixed RMS accuracy          : %.2f %%\n',accuracy_fixed);

Fixed RMS accuracy    : 92.44 %

fprintf('Adaptive RMS accuracy: %.2f %%\n',accuracy_adaptive);

Adaptive RMS accuracy: 96.38 %

%% 17. DISPLAY DETECTION AND RECOVERY DELAYS

fprintf('\n========================================\n');

fprintf('DETECTION AND RECOVERY DELAYS\n');

fprintf('========================================\n');

for k = 1:3

fprintf('\n%s:\n',event_names{k});

fprintf('Detection delay : %.2f ms\n',detection_delay(k));

fprintf('Recovery delay      : %.2f ms\n',recovery_delay(k));

end

Sag:
Swell:
Interruption:

%% 18. INPUT VOLTAGE PLOT

figure('Name','Input Voltage','Color','w');

plot(t,v,'b','LineWidth',1);

grid on;

xlabel('Time (s)');
ylabel('Voltage (V)');

title('Simulated Voltage with Sag, Swell and Interruption');

xline(sag_start,'--r','Sag start');
xline(sag_end,'--r','Sag end');

xline(swell_start,'--m','Swell start');
xline(swell_end,'--m','Swell end');

xline(int_start,'--k','Interruption start');
xline(int_end,'--k','Interruption end');

%% 19. FIXED VS ADAPTIVE RMS PLOT

figure('Name','RMS Comparison','Color','w');

plot(t,Vrms_fixed,'b','LineWidth',1.2);
hold on;

plot(t,Vrms_adaptive,'r','LineWidth',1.2);

yline(Vrms_nominal,'--k','Nominal RMS');

grid on;

xlabel('Time (s)');
ylabel('RMS Voltage (V)');

title('Fixed vs Adaptive RMS Estimation');

legend('Fixed 200-sample RMS', ...
'Adaptive RMS','Nominal RMS', ...
'Location','best');

%% 20. ADAPTIVE WINDOW SELECTION PLOT

figure('Name','Window Selection','Color','w');

stairs(t,selected_window,'LineWidth',1.2);

grid on;

xlabel('Time (s)');
ylabel('Window Length (samples)');

title('Adaptive RMS Window Selection');

yticks([100 200]);
ylim([0 220]);

%% 21. PER-UNIT RMS PLOT

figure('Name','Adaptive RMS Per Unit','Color','w');

plot(t,Vpu_adaptive,'r','LineWidth',1.2);

hold on;

yline(0.1,'--k','Interruption threshold');
yline(0.9,'--m','Sag threshold');
yline(1.1,'--g','Swell threshold');

grid on;

xlabel('Time (s)');
ylabel('Voltage (pu)');

title('Adaptive RMS Voltage in Per-Unit');

ylim([0 1.5]);

%% 22. CLASSIFICATION PLOT

figure('Name','Adaptive Classification','Color','w');

stairs(t,condition_adaptive,'LineWidth',1.2);

grid on;

xlabel('Time (s)');
ylabel('Detected Condition');

title('Adaptive RMS Disturbance Classification');

yticks([0 1 2 3]);

yticklabels({'Normal','Sag','Swell','Interruption'});

ylim([-0.5 3.5]);

%% 23. RMS ESTIMATION ERROR PLOT

figure('Name','RMS Estimation Error','Color','w');

plot(t(valid), ...
Vrms_fixed(valid)-Vrms_ideal(valid), ...
'b','LineWidth',1);

hold on;

plot(t(valid), ...
Vrms_adaptive(valid)-Vrms_ideal(valid), ...
'r','LineWidth',1);

yline(0,'--k');

grid on;

xlabel('Time (s)');
ylabel('Estimation Error (V)');

title('Fixed vs Adaptive RMS Estimation Error');

legend('Fixed RMS','Adaptive RMS');

%% 24. RESULTS SUMMARY

fprintf('\n========================================\n');

fprintf('SIMULATION COMPLETED SUCCESSFULLY\n');

fprintf('========================================\n');

fprintf('\nStable-region RMS values:\n');


stable_normal = t >= 0.05 & t < 0.15;
stable_sag = t >= 0.25 & t < 0.30;
stable_swell = t >= 0.55 & t < 0.60;
stable_interruption = t >= 0.78 & t < 0.82;

fprintf('Normal       : %.2f V\n', ...
mean(Vrms_adaptive(stable_normal),'omitnan'));


fprintf('Sag          : %.2f V\n', ...
mean(Vrms_adaptive(stable_sag),'omitnan'));


fprintf('Swell        : %.2f V\n', ...
mean(Vrms_adaptive(stable_swell),'omitnan'));


fprintf('Interruption : %.2f V\n', ...
mean(Vrms_adaptive(stable_interruption),'omitnan'));
