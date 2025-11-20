%% ===== CREATE FAKE "UNREAL" DATASETS (everything high & similar) =====
close all; clearvars; clc

% --- Time axis and windows (same as your GLM script) ---
res  = 0.01;
bins = -1:res:4;

BL_win = [-0.1  0.0];
A_win  = [ 0.4  0.5];
B_win  = [ 0.9  1.0];
C_win  = [ 1.4  1.5];
D_win  = [ 1.9  2.0];

Ai  = bins>=A_win(1) & bins< A_win(2);
Bi  = bins>=B_win(1) & bins< B_win(2);
Ci  = bins>=C_win(1) & bins< C_win(2);
Di  = bins>=D_win(1) & bins< D_win(2);

% --- Sizes (match typical OE20) ---
n_ABCD = 130;
n_DCBA = 30;
n_AAAA = 40;

% --- Baseline + element strengths ---
p_base = 0.002;  % low background activity outside windows
pA = 0.10;       % 
pB = 0.30;       % 
pC = 0.40;       % 
pD = 0.50;       % 

noise_sd = 0.01; % mild trial noise

% --- Helper to create one trial row given per-window levels [A,B,C,D] ---
make_trial = @(levels) max(0, min(0.999, ...
    p_base*ones(1, numel(bins)) + ...
    levels(1)*Ai + levels(2)*Bi + levels(3)*Ci + levels(4)*Di + ...
    noise_sd*randn(1, numel(bins)) ));

%% ===== ABCD ===== (positions map to A,B, C, D)
OEFAKE_FMR1FAKE_ABCD_trialbytrial_licking = zeros(n_ABCD, numel(bins));
for t = 1:n_ABCD
    OEFAKE_FMR1FAKE_ABCD_trialbytrial_licking(t,:) = make_trial([pA, pB, pC, pD]);
end
save('OEFAKE_FMR1FAKE_ABCD_trialbytrial_licking.mat', ...
     'OEFAKE_FMR1FAKE_ABCD_trialbytrial_licking');

%% ===== DCBA ===== (positions map to D, C, B, A)
OEFAKE_FMR1FAKE_DCBA_trialbytrial_licking = zeros(n_DCBA, numel(bins));
for t = 1:n_DCBA
    OEFAKE_FMR1FAKE_DCBA_trialbytrial_licking(t,:) = make_trial([pD, pC, pB, pA]);
end
save('OEFAKE_FMR1FAKE_DCBA_trialbytrial_licking.mat', ...
     'OEFAKE_FMR1FAKE_DCBA_trialbytrial_licking');

%% ===== AAAA ===== (all positions high like A/B/C/D)
OEFAKE_FMR1FAKE_AAAA_trialbytrial_licking = zeros(n_AAAA, numel(bins));
for t = 1:n_AAAA
    OEFAKE_FMR1FAKE_AAAA_trialbytrial_licking(t,:) = make_trial([pA, pB, pC, pD]);
end
save('OEFAKE_FMR1FAKE_AAAA_trialbytrial_licking.mat', ...
     'OEFAKE_FMR1FAKE_AAAA_trialbytrial_licking');
