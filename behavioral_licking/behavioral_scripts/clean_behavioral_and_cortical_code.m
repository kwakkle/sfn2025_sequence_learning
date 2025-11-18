% Naive vs Learned w/ SEM

close all
clear all
clc

load ("FMTEL.mat")
load ("FMTNL.mat")

shadetype = 'SEM';
Time = -1:0.01:4;
sm_val = 20;

figure
hold on

% [NEW] Shaded backgrounds

y_limits = [0 0.25];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D
fill([2 3 3 2],         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Response

text(-0.7, 0.125, 'Expert', 'Color', [0.4 0.7 1], 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left')  % light blue
text(-0.7, 0.045, 'Naive', 'Color', [0.6 1 0.6], 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left')  % light green

text(0.15, 0.245, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.245, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.245, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.245, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(2.5, 0.245, 'RW', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

box off
grid on

text(-0.65, 0.21, '$$\mathbf{p^{cond} < 0.01^*}$$', 'Interpreter', 'latex', ...
    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.65, 0.19, '$$\mathbf{p^{stim} < 0.01^*}$$', 'Interpreter', 'latex', ...
    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.65, 0.17, '$$\mathbf{p^{inter} < 0.01^*}$$', 'Interpreter', 'latex', ...
    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

% rest is the same -----

stdshade_smooth3(FMTNL, Time, 0.5, [0.4 0.7 1], sm_val, 'SEM'); % Naive
stdshade_smooth3(FMTEL, Time, 0.5, [0.6 1 0.6], sm_val, 'SEM'); % Expert

title('Wild Type - Licking Across Sessions: Naive vs Expert')

xlim([-1 4])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Lick Probability')
%% Wild Type: Expert vs Reversal vs Oddball
close all
clear all
clc

load ("FMTEL.mat")
load ("FMTRL.mat")
load ("FMTOL.mat")

shadetype = 'SEM';
sm_val = 20;
Time = -1:0.01:4; % Ensure this matches your actual data resolution

figure
hold on

% Shaded stimulus backgrounds
y_limits = [0 0.25];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D
fill([2 3 3 2],         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Response

% Plot smoothed traces with SEM
stdshade_smooth3(FMTEL, Time, 0.5, 'r', sm_val, 'SEM'); % Expert
stdshade_smooth3(FMTOL, Time, 0.5, [0.90 0.70 0.10], sm_val, shadetype); % Reversal AL
stdshade_smooth3(FMTRL, Time, 0.5, [0.45 0.00 0.70], sm_val, shadetype); % Reversal RL

% Labels for each stimulus
text(2.5, 0.245, 'RW', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(-0.50, 0.128, {'Standard', '(ABCD)'}, 'Color', 'r', 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(2, 0.04, {'Oddball', '(AAAA)'}, 'Color', [0.90 0.70 0.10], 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(1.2, 0.215, {'Reversal', '(DCBA)'}, 'Color', [0.45 0.00 0.70], 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')

text(-0.8, 0.21, '$$\mathbf{p^{seq} \;\;< 0.01^*}$$', 'Interpreter', 'latex', ...
    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.8, 0.19, '$$\mathbf{p^{stim} \;< 0.01^*}$$', 'Interpreter', 'latex', ...
    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.8, 0.17, '$$\mathbf{p^{inter} > 0.05}$$', 'Interpreter', 'latex', ...
    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

xlim([-1 4])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Lick Probability')
title('Wild Type - Licking Across Sessions: Test vs Reversal vs Oddball')

box off
grid on


%% FXS: Expert vs Reversal vs Oddball

close all
clear all
clc

load ("FMTFXSEL.mat")
load ("FMTFXSAL.mat")
load ("FMTFXSRL.mat")

shadetype = 'SEM';
sm_val = 20;
Time = -1:0.01:4; % Ensure this matches your actual data resolution

figure
hold on

% Shaded stimulus backgrounds
y_limits = [0 0.25];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D
fill([2 3 3 2],         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Response

% Plot smoothed traces with SEM
stdshade_smooth3(FMTFXSEL, Time, 0.5, 'r', sm_val, shadetype); % Expert
stdshade_smooth3(FMTFXSAL, Time, 0.5, [0.45 0.00 0.70], sm_val, shadetype); % Reversal AL
stdshade_smooth3(FMTFXSRL, Time, 0.5, [0.90 0.70 0.10], sm_val, shadetype); % Reversal RL

% Labels for each stimulus
text(2.5, 0.245, 'RW', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(-0.40, 0.110, {'Standard', '(ABCD)'}, 'Color', 'r', 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(1, 0.04, {'Oddball', '(AAAA)'}, 'Color', [0.90 0.70 0.10], 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(1.2, 0.185, {'Reversal', '(DCBA)'}, 'Color', [0.45 0.00 0.70], 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')

xlim([-1 4])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Lick Probability')
title('FMR1 KO - Licking Across Sessions: Test vs Reversal vs Oddball')

box off
grid on


%% FXS Decode

% Naive vs Learned FXS Comparison with SEM

close all
clear all
clc

load ("FMTFXSNL.mat")
load ("FMTFXSEL.mat")

shadetype = 'SEM';
sm_val = 20;
Time = -1:0.01:4; % Ensure Time is defined; update if needed

figure
hold on

% Shaded stimulus backgrounds
y_limits = [0 0.25];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D
fill([2 3 3 2],         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Response

% Plot data with SEM
stdshade_smooth3(FMTFXSNL, Time, 0.5, 'k', sm_val, shadetype); % Naive
stdshade_smooth3(FMTFXSEL, Time, 0.5, 'r', sm_val, shadetype); % Expert

% Stimulus labels
text(0.15, 0.245, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.245, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.245, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.245, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(2.5, 0.245, 'RW', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(-0.8, 0.11, 'Expert', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(-0.8, 0.03, 'Naive', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')

text(-0.65, 0.21, '$$\mathbf{p^{cond} < 0.01^*}$$', 'Interpreter', 'latex', ...
    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.65, 0.19, '$$\mathbf{p^{stim} > 0.05}$$', 'Interpreter', 'latex', ...
    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.65, 0.17, '$$\mathbf{p^{inter} > 0.05}$$', 'Interpreter', 'latex', ...
    'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

xlim([-1 4])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Lick Probability')

box off
grid on

title('FMR1 KO - Licking Across Sessions: Naive vs Expert')

%% WT vs FXS

% Wild type vs FXS


close all
clear all
clc

load ("FMTEL.mat")
load ("FMTFXSEL.mat")

shadetype = 'SEM';
sm_val = 20;
Time = -1:0.01:4; % Ensure Time is defined; update if needed

figure
hold on

% Shaded stimulus backgrounds
y_limits = [0 0.25];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D
fill([2 3 3 2],         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Response

% Plot data with SEM
stdshade_smooth3(FMTEL, Time, 0.5, 'r', sm_val, shadetype); % Naive
stdshade_smooth3(FMTFXSEL, Time, 0.5, 'g', sm_val, shadetype); % Expert

% Stimulus labels
text(0.15, 0.245, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.245, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.245, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.245, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(2.5, 0.245, 'RW', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(-0.8, 0.12, 'Wild Type', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(-0.6, 0.058, 'FMR1 KO', 'Color', [0 0.5 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')


xlim([-1 4])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Lick Probability')

box off
grid on

title('Expert Mice Licking Across Sessions: WT vs FMR1 KO')



%% ------ Imaging ---- Sensory Cortex

close all
clear all
clc

load ("FXSMTEIS.mat")
load ("FXSMTNIS.mat")

Time = (0.1:0.1:3)-1;

% Naive vs Expert FXS (Sensory) - Avg. with SEM

figure
hold on

shadetype = 'SEM';
sm_val = 1;  % Smoothing value for sensory

% Shaded stimulus backgrounds
y_limits = [0 0.025];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D

% Main traces
stdshade_smooth3(FXSMTNIS, Time, 0.5, 'k', sm_val, shadetype); % Naive
stdshade_smooth3(FXSMTEIS, Time, 0.5, 'r', sm_val, shadetype); % Expert

% Stimulus text labels
text(0.15, 0.0245, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.0245, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.0245, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.0245, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(1, 0.005, 'Expert', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(0.3, 0.017, 'Naive',  'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')

% LaTeX-style stats
text(-0.4, 0.02, '$$\mathbf{p^{cond} = 0.07}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.4, 0.018, '$$\mathbf{p^{stim} < 0.01^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.4, 0.016, '$$\mathbf{p^{inter} = 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

% Final formatting
xlim([-0.5 2])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Neural Response')
title('FMR1 KO Sensory Cortex: Naive vs Expert')
box off
grid on



%% Wild Type: Sensory Cortex


close all
clear all
clc

load ("MTNSI.mat")
load ("MTESI.mat")

Time = (0.1:0.1:3)-1;

% Naive vs Expert FXS (Sensory) - Avg. with SEM

figure
hold on

shadetype = 'SEM';
sm_val = 1;  % Smoothing value for sensory

% Shaded stimulus backgrounds
y_limits = [0 0.025];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D

% Main traces
stdshade_smooth3(MTNSI, Time, 0.5, 'k', 1, shadetype);  
stdshade_smooth3(MTESI, Time, 0.5, 'r', 1, shadetype);

% Stimulus text labels
text(0.15, 0.0245, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.0245, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.0245, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.0245, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(0.55, 0.009, 'Expert', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(0.3, 0.017, 'Naive',  'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')

% LaTeX-style stats
text(-0.4, 0.02, '$$\mathbf{p^{cond} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.4, 0.018, '$$\mathbf{p^{stim} < 0.01^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.4, 0.016, '$$\mathbf{p^{inter} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

% Final formatting
xlim([-0.5 2])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Neural Response')
title('Wild Type Sensory Cortex: Naive vs Expert')
box off
grid on
%% Naive: Wild Type vs FXS Sensory Cortex


close all
clear all
clc

load ("MTNSI.mat")
load ("FXSMTNIS.mat")

Time = (0.1:0.1:3)-1;

% Naive vs Expert FXS (Sensory) - Avg. with SEM

figure
hold on

shadetype = 'SEM';
sm_val = 1;  % Smoothing value for sensory

% Shaded stimulus backgrounds
y_limits = [0 0.025];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D

% Main traces
stdshade_smooth3(MTNSI, Time, 0.5, 'b', 1, shadetype);  
stdshade_smooth3(FXSMTNIS, Time, 0.5, 'g', 1, shadetype);

% Stimulus text labels
text(0.15, 0.0245, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.0245, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.0245, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.0245, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(0.55, 0.01, 'FMR1 KO', 'Color', [0 0.5 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(0.3, 0.017, 'WT',  'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')

% LaTeX-style stats
text(-0.4, 0.02, '$$\mathbf{p^{geno} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.4, 0.018, '$$\mathbf{p^{stim} < 0.01^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.4, 0.016, '$$\mathbf{p^{inter} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

% Final formatting
xlim([-0.5 2])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Neural Response')
title('Naive Sensory Cortex: Wild Type vs FXS')
box off
grid on


%% Expert: Wild Type vs FXS Sensory Cortex


close all
clear all
clc

load ("MTESI.mat")
load ("FXSMTEIS.mat")

Time = (0.1:0.1:3)-1;

% Naive vs Expert FXS (Sensory) - Avg. with SEM

figure
hold on

shadetype = 'SEM';
sm_val = 1;  % Smoothing value for sensory

% Shaded stimulus backgrounds
y_limits = [0 0.025];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D

% Main traces
stdshade_smooth3(MTESI, Time, 0.5, 'b', 1, shadetype);  
stdshade_smooth3(FXSMTEIS, Time, 0.5, 'g', 1, shadetype);

% Stimulus text labels
text(0.15, 0.0245, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.0245, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.0245, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.0245, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(0.7, 0.011, 'FMR1 KO', 'Color', [0 0.5 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(0.35, 0.014, 'WT',  'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')

% LaTeX-style stats
text(-0.4, 0.02, '$$\mathbf{p^{geno} < 0.05^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.4, 0.018, '$$\mathbf{p^{stim} < 0.01^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.4, 0.016, '$$\mathbf{p^{inter} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

% Final formatting
xlim([-0.5 2])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Neural Response')
title('Expert Sensory Cortex: Wild Type vs FMR1 KO')
box off
grid on


%% Frontal -------- 

% Wild Type Frontal Cortex: Naive vs Expert

close all
clear all
clc

load ("MTNFI.mat")
load ("MTEFI.mat")

Time = (0.1:0.1:3)-1;

figure
hold on

shadetype = 'SEM';
sm_val = 1;  % Smoothing value for sensory

% Shaded stimulus backgrounds
y_limits = [0 0.01];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D

% Main traces
stdshade_smooth3(MTNFI, Time, 0.5, 'k', 1, shadetype);  
stdshade_smooth3(MTEFI, Time, 0.5, 'r', 1, shadetype);

% Stimulus text labels
text(0.15, 0.0098, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.0098, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.0098, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.0098, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(0.55, 0.006, 'Expert', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(0.55, 0.002, 'Naive',  'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')

% LaTeX-style stats
text(-0.3, 0.007, '$$\mathbf{p^{cond} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.3, 0.006, '$$\mathbf{p^{stim} < 0.01^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.3, 0.005, '$$\mathbf{p^{inter} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

% Final formatting
xlim([-0.5 2])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Neural Response')
title('Wild Type Frontal Cortex: Naive vs Expert')
box off
grid on

%% Frontal


close all
clear all
clc

load ("FXSMTNIF.mat")
load ("FXSMTEIF.mat")

Time = (0.1:0.1:3)-1;

figure
hold on

shadetype = 'SEM';
sm_val = 1;  % Smoothing value for sensory

% Shaded stimulus backgrounds
y_limits = [0 0.01];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D

% Main traces
stdshade_smooth3(FXSMTNIF, Time, 0.5, 'k', 1, shadetype);  
stdshade_smooth3(FXSMTEIF, Time, 0.5, 'r', 1, shadetype);

% Stimulus text labels
text(0.15, 0.0098, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.0098, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.0098, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.0098, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(0.5, 0.003, 'Expert', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(0.1, 0.003, 'Naive',  'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')

% LaTeX-style stats
text(-0.3, 0.007, '$$\mathbf{p^{cond} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.3, 0.006, '$$\mathbf{p^{stim} < 0.01^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.3, 0.005, '$$\mathbf{p^{inter} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

% Final formatting
xlim([-0.5 2])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Neural Response')
title('FMR1 KO Frontal Cortex: Naive vs Expert')
box off
grid on

%% Naive: Wild Type vs FXS Frontal Cortex


close all
clear all
clc

load ("MTNFI.mat")
load ("FXSMTNIF.mat")

Time = (0.1:0.1:3)-1;

% Naive vs Expert FXS (Sensory) - Avg. with SEM

figure
hold on

shadetype = 'SEM';
sm_val = 1;  % Smoothing value for sensory

% Shaded stimulus backgrounds
y_limits = [0 0.01];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D

% Main traces
stdshade_smooth3(MTNFI, Time, 0.5, 'b', 1, shadetype);  
stdshade_smooth3(FXSMTNIF, Time, 0.5, 'g', 1, shadetype);

% Stimulus text labels
text(0.15, 0.0098, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.0098, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.0098, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.0098, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(0.5, 0.002, 'FMR1 KO', 'Color', [0 0.5 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(0.6, 0.0045, 'WT',  'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')

% LaTeX-style stats
text(-0.3, 0.007, '$$\mathbf{p^{geno} < 0.01^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.3, 0.006, '$$\mathbf{p^{stim} < 0.01^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.3, 0.005, '$$\mathbf{p^{inter} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')


% Final formatting
xlim([-0.5 2])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Neural Response')
title('Naive Frontal Cortex: Wild Type vs FMR1 KO')
box off
grid on


%% Expert: Wild Type vs FXS Frontal Cortex



close all
clear all
clc

load ("MTEFI.mat")
load ("FXSMTEIF.mat")

Time = (0.1:0.1:3)-1;

% Naive vs Expert FXS (Sensory) - Avg. with SEM

figure
hold on

shadetype = 'SEM';
sm_val = 1;  % Smoothing value for sensory

% Shaded stimulus backgrounds
y_limits = [0 0.01];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D

% Main traces
stdshade_smooth3(MTEFI, Time, 0.5, 'b', 1, shadetype);  
stdshade_smooth3(FXSMTEIF, Time, 0.5, 'g', 1, shadetype);

% Stimulus text labels
text(0.15, 0.0098, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.0098, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.0098, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.0098, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(0.4, 0.0025, 'FMR1 KO', 'Color', [0 0.5 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')
text(0.5, 0.0048, 'WT',  'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left')

% LaTeX-style stats
text(-0.3, 0.007, '$$\mathbf{p^{geno} < 0.05^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.3, 0.006, '$$\mathbf{p^{stim} < 0.01^*}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

text(-0.3, 0.005, '$$\mathbf{p^{inter} > 0.05}$$', 'Interpreter', 'latex', ...
     'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left')

% Final formatting
xlim([-0.5 2])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Neural Response')
title('Expert Frontal Cortex: Wild Type vs FMR1 KO')
box off
grid on

%% reverse WT vs reverse FXS

close all
clear all
clc

load ("FMTRL.mat")
load ("FMTFXSRL.mat")

shadetype = 'SEM';
sm_val = 20;
Time = -1:0.01:4; % Ensure this matches your actual data resolution

figure
hold on

% Shaded stimulus backgrounds
y_limits = [0 0.25];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D
fill([2 3 3 2],         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Response

% Plot smoothed traces with SEM
stdshade_smooth3(FMTFXSRL, Time, 0.5, 'g', sm_val, 'SEM'); % Expert
stdshade_smooth3(FMTRL, Time, 0.5, 'r', sm_val, shadetype); % Reversal RL

% Labels for each stimulus
text(0.15, 0.245, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.245, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.245, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.245, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(2.5, 0.245, 'RW', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(1.2, 0.05, {'FMR1 KO', '(DCBA)'}, 'Color', [0 0.5 0], 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(1.2, 0.215, {'WT', '(DCBA)'}, 'Color', 'r', 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')

xlim([-1 4])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Lick Probability')
title('Reversal Licking Across Sessions: WT vs FMR1 KO')

box off
grid on
%% oddball WT vs reverse FXS

close all
clear all
clc

load ("FMTOL.mat")
load ("FMTFXSAL.mat")

shadetype = 'SEM';
sm_val = 20;
Time = -1:0.01:4; % Ensure this matches your actual data resolution

figure
hold on

% Shaded stimulus backgrounds
y_limits = [0 0.25];
fill([0 0.3 0.3 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % A
fill([0.5 0.8 0.8 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor', 'none', 'FaceAlpha', 0.3); % B
fill([1 1.3 1.3 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor', 'none', 'FaceAlpha', 0.3); % C
fill([1.5 1.8 1.8 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor', 'none', 'FaceAlpha', 0.3); % D
fill([2 3 3 2],         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Response

% Plot smoothed traces with SEM
stdshade_smooth3(FMTFXSAL, Time, 0.5, 'g', sm_val, 'SEM'); % Expert
stdshade_smooth3(FMTOL, Time, 0.5, 'r', sm_val, shadetype); % Reversal RL

% Labels for each stimulus
text(0.15, 0.245, 'A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.65, 0.245, 'B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.15, 0.245, 'C', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.65, 0.245, 'D', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(2.5, 0.245, 'RW', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

% Group labels
text(1.2, 0.2, {'FMR1 KO', '(AAAA)'}, 'Color', [0 0.5 0], 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')
text(1.2, 0.04, {'WT', '(AAAA)'}, 'Color', 'r', 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center')

xlim([-1 4])
ylim(y_limits)
xlabel('Time (s)')
ylabel('Lick Probability')
title('Oddball Licking Across Sessions: WT vs FMR1 KO')

box off
grid on



