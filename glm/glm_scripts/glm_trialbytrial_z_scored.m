%% GLM Trial-by-Trial on testing day using z-scores – Wild Type Mice (Gaussian GLM, no baseline covariate)
close all
clearvars
clc

% ===== Load z-scored trial-by-trial licking data =====
load('OE18_ABCD_trialbytrial_licking_z.mat');  
load('OE18_DCBA_trialbytrial_licking_z.mat');   
load('OE18_AAAA_trialbytrial_licking_z.mat'); 

load('OE20_ABCD_trialbytrial_licking_z.mat');  
load('OE20_DCBA_trialbytrial_licking_z.mat');   
load('OE20_AAAA_trialbytrial_licking_z.mat');  

load('OE21_ABCD_trialbytrial_licking_z.mat');  
load('OE21_DCBA_trialbytrial_licking_z.mat');   
load('OE21_AAAA_trialbytrial_licking_z.mat');  

load('OE24_ABCD_trialbytrial_licking_z.mat');  
load('OE24_DCBA_trialbytrial_licking_z.mat');   
load('OE24_AAAA_trialbytrial_licking_z.mat');  

load('OE32_ABCD_trialbytrial_licking_z.mat');  
load('OE32_DCBA_trialbytrial_licking_z.mat');   
load('OE32_AAAA_trialbytrial_licking_z.mat');  

load('OE33_ABCD_trialbytrial_licking_z.mat');  
load('OE33_DCBA_trialbytrial_licking_z.mat');   
load('OE33_AAAA_trialbytrial_licking_z.mat');  

load('OE34_ABCD_trialbytrial_licking_z.mat');  
load('OE34_DCBA_trialbytrial_licking_z.mat');   
load('OE34_AAAA_trialbytrial_licking_z.mat');  

load('OE35_ABCD_trialbytrial_licking_z.mat');  
load('OE35_DCBA_trialbytrial_licking_z.mat');   
load('OE35_AAAA_trialbytrial_licking_z.mat');  

load('OE39_ABCD_trialbytrial_licking_z.mat');  
load('OE39_DCBA_trialbytrial_licking_z.mat');   
load('OE39_AAAA_trialbytrial_licking_z.mat');  

load('OE40_ABCD_trialbytrial_licking_z.mat');  
load('OE40_DCBA_trialbytrial_licking_z.mat');   
load('OE40_AAAA_trialbytrial_licking_z.mat');  

load('OE45_ABCD_trialbytrial_licking_z.mat');  
load('OE45_DCBA_trialbytrial_licking_z.mat');   
load('OE45_AAAA_trialbytrial_licking_z.mat');  

load('OE46_ABCD_trialbytrial_licking_z.mat');  
load('OE46_DCBA_trialbytrial_licking_z.mat');   
load('OE46_AAAA_trialbytrial_licking_z.mat');  

load('OE47_ABCD_trialbytrial_licking_z.mat');  
load('OE47_DCBA_trialbytrial_licking_z.mat');   
load('OE47_AAAA_trialbytrial_licking_z.mat');

% ===== Combine all WT mice =====
WT_ABCD_only_z_data = [
    OE18_ABCD_trialbytrial_licking_z;
    OE20_ABCD_trialbytrial_licking_z;
    OE21_ABCD_trialbytrial_licking_z;
    OE24_ABCD_trialbytrial_licking_z;
    OE32_ABCD_trialbytrial_licking_z;
    OE33_ABCD_trialbytrial_licking_z;
    OE34_ABCD_trialbytrial_licking_z;
    OE35_ABCD_trialbytrial_licking_z;
    OE39_ABCD_trialbytrial_licking_z;
    OE40_ABCD_trialbytrial_licking_z;
    OE45_ABCD_trialbytrial_licking_z;
    OE46_ABCD_trialbytrial_licking_z;
    OE47_ABCD_trialbytrial_licking_z
];

WT_AAAA_only_z_data = [
    OE18_AAAA_trialbytrial_licking_z;
    OE20_AAAA_trialbytrial_licking_z;
    OE21_AAAA_trialbytrial_licking_z;
    OE24_AAAA_trialbytrial_licking_z;
    OE32_AAAA_trialbytrial_licking_z;
    OE33_AAAA_trialbytrial_licking_z;
    OE34_AAAA_trialbytrial_licking_z;
    OE35_AAAA_trialbytrial_licking_z;
    OE39_AAAA_trialbytrial_licking_z;
    OE40_AAAA_trialbytrial_licking_z;
    OE45_AAAA_trialbytrial_licking_z;
    OE46_AAAA_trialbytrial_licking_z;
    OE47_AAAA_trialbytrial_licking_z
];

WT_DCBA_only_z_data = [
    OE18_DCBA_trialbytrial_licking_z;
    OE20_DCBA_trialbytrial_licking_z;
    OE21_DCBA_trialbytrial_licking_z;
    OE24_DCBA_trialbytrial_licking_z;
    OE32_DCBA_trialbytrial_licking_z;
    OE33_DCBA_trialbytrial_licking_z;
    OE34_DCBA_trialbytrial_licking_z;
    OE35_DCBA_trialbytrial_licking_z;
    OE39_DCBA_trialbytrial_licking_z;
    OE40_DCBA_trialbytrial_licking_z;
    OE45_DCBA_trialbytrial_licking_z;
    OE46_DCBA_trialbytrial_licking_z;
    OE47_DCBA_trialbytrial_licking_z
];

% ===== Define analysis windows =====
res  = 0.01;                    
bins = -1:res:4;               

BL_win = [-0.1  0.0];
A_win  = [ 0.4  0.5];
B_win  = [ 0.9  1.0];
C_win  = [ 1.4  1.5];
D_win  = [ 1.9  2.0];

BLi = bins>=BL_win(1) & bins< BL_win(2);
Ai  = bins>=A_win(1)  & bins< A_win(2);
Bi  = bins>=B_win(1)  & bins< B_win(2);
Ci  = bins>=C_win(1)  & bins< C_win(2);
Di  = bins>=D_win(1)  & bins< D_win(2);

seqs     = {'ABCD','DCBA','AAAA'};
datasets = {WT_ABCD_only_z_data, WT_DCBA_only_z_data, WT_AAAA_only_z_data};
winIdx   = {Ai, Bi, Ci, Di;   ... % ABCD: A,B,C,D
            Di, Ci, Bi, Ai;   ... % DCBA: D,C,B,A
            Ai, Ai, Ai, Ai};      % AAAA: A,A,A,A
elems    = {['A';'B';'C';'D'], ...
            ['D';'C';'B';'A'], ...
            ['A';'A';'A';'A']};

% ===== Build trial-by-trial table (response = mean z per window) =====
nPerSet = cellfun(@(X) size(X,1), datasets);
rows    = sum(nPerSet) * 4;

Mouse   = strings(rows,1);
Session = strings(rows,1);
TrialNum= zeros(rows,1);
Seq     = strings(rows,1);
Elem    = strings(rows,1);
Pos     = zeros(rows,1);

Count   = zeros(rows,1);   % mean z-score in the window

row = 0;
for s = 1:numel(seqs)
    X = datasets{s};             
    nT = size(X,1);
    for t = 1:nT
        tr = X(t,:);  % one trial (z-scored time series)

        % four stimulus windows per trial
        for pos = 1:4
            row = row + 1;

            idx = winIdx{s,pos};

            Mouse(row)   = sprintf("OE20_T%04d", t);
            Session(row) = "S1";
            TrialNum(row)= t;
            Seq(row)     = seqs{s};
            Elem(row)    = elems{s}(pos);
            Pos(row)     = pos;

            % mean z-score in the window (continuous response)
            Count(row) = mean(tr(idx));
        end
    end
end

% ===== Assemble table =====
T = table( categorical(Mouse), categorical(Session), TrialNum, ...
           categorical(Seq), categorical(Elem), Pos, Count, ...
           'VariableNames', {'Mouse','Session','TrialNum','Seq','Elem','Pos','Count'});

% Derived covariate: optional within-sequence drift (keep or drop)
T.TrialC = zscore(double(T.TrialNum));

% Factor coding
T.ElemCat = categorical(T.Elem);
T.ElemCat = reordercats(T.ElemCat, {'A','B','C','D'}); 

T.PosCat  = categorical(string(T.Pos));
T.PosCat  = reordercats(T.PosCat, {'1','2','3','4'});    

% ===== Fit Gaussian GLMs on z-scores (no baseline covariate) =====
mdl       = fitglm(T,'Count ~ 1 + ElemCat + PosCat + TrialC', 'Distribution','normal','Link','identity');


mdl_noPos = fitglm(T,'Count ~ 1 + ElemCat + TrialC',          'Distribution','normal','Link','identity');
mdl_noElem= fitglm(T,'Count ~ 1 + PosCat  + TrialC',          'Distribution','normal','Link','identity');

% ===== LR-style comparisons =====
dDev_pos  = max(mdl_noPos.Deviance  - mdl.Deviance,  0);
dDev_elem = max(mdl_noElem.Deviance - mdl.Deviance,  0);

df_pos  = max(height(mdl.Coefficients) - height(mdl_noPos.Coefficients),1);
df_elem = max(height(mdl.Coefficients) - height(mdl_noElem.Coefficients),1);

p_pos   = 1 - chi2cdf(dDev_pos,  df_pos);
p_elem  = 1 - chi2cdf(dDev_elem, df_elem);

n       = height(T);
R2_pos  = 1 - exp(-dDev_pos  / n);
R2_elem = 1 - exp(-dDev_elem / n);
w_pos   = R2_pos / max(R2_pos + R2_elem, eps);
w_elem  = 1 - w_pos;

fprintf('p_pos = %.3g (df=%d), p_elem = %.3g (df=%d) | R2_pos = %.3f, R2_elem = %.3f | w_pos = %.2f, w_elem = %.2f\n', ...
        p_pos, df_pos, p_elem, df_elem, R2_pos, R2_elem, w_pos, w_elem);

disp(mdl.Coefficients)


%% licking probability (average PSTHs + legend fix)

load('OE18_testing_ABCD_lick_z.mat');
load('OE18_testing_AAAA_lick_z.mat');
load('OE18_testing_DCBA_lick_z.mat');

load('OE20_testing_ABCD_lick_z.mat');
load('OE20_testing_AAAA_lick_z.mat');
load('OE20_testing_DCBA_lick_z.mat');

load('OE21_testing_ABCD_lick_z.mat');
load('OE21_testing_AAAA_lick_z.mat');
load('OE21_testing_DCBA_lick_z.mat');

load('OE24_testing_ABCD_lick_z.mat');
load('OE24_testing_AAAA_lick_z.mat');
load('OE24_testing_DCBA_lick_z.mat');

load('OE32_testing_ABCD_lick_z.mat');
load('OE32_testing_AAAA_lick_z.mat');
load('OE32_testing_DCBA_lick_z.mat');

load('OE33_testing_ABCD_lick_z.mat');
load('OE33_testing_AAAA_lick_z.mat');
load('OE33_testing_DCBA_lick_z.mat');

load('OE34_testing_ABCD_lick_z.mat');
load('OE34_testing_AAAA_lick_z.mat');
load('OE34_testing_DCBA_lick_z.mat');

load('OE35_testing_ABCD_lick_z.mat');
load('OE35_testing_AAAA_lick_z.mat');
load('OE35_testing_DCBA_lick_z.mat');

load('OE39_testing_ABCD_lick_z.mat');
load('OE39_testing_AAAA_lick_z.mat');
load('OE39_testing_DCBA_lick_z.mat');

load('OE40_testing_ABCD_lick_z.mat');
load('OE40_testing_AAAA_lick_z.mat');
load('OE40_testing_DCBA_lick_z.mat');

load('OE45_testing_ABCD_lick_z.mat');
load('OE45_testing_AAAA_lick_z.mat');
load('OE45_testing_DCBA_lick_z.mat');

load('OE46_testing_ABCD_lick_z.mat');
load('OE46_testing_AAAA_lick_z.mat');
load('OE46_testing_DCBA_lick_z.mat');

load('OE47_testing_ABCD_lick_z.mat');
load('OE47_testing_AAAA_lick_z.mat');
load('OE47_testing_DCBA_lick_z.mat');

% --- Combine ALL mice (WT) for ABCD / AAAA / DCBA (z-scored testing averages) ---

WT_ABCD_all_z = [
    OE18_testing_ABCD_lick_z;
    OE20_testing_ABCD_lick_z;
    OE21_testing_ABCD_lick_z;
    OE24_testing_ABCD_lick_z;
    OE32_testing_ABCD_lick_z;
    OE33_testing_ABCD_lick_z;
    OE34_testing_ABCD_lick_z;
    OE35_testing_ABCD_lick_z;
    OE39_testing_ABCD_lick_z;
    OE40_testing_ABCD_lick_z;
    OE45_testing_ABCD_lick_z;
    OE46_testing_ABCD_lick_z;
    OE47_testing_ABCD_lick_z
];

WT_AAAA_all_z = [
    OE18_testing_AAAA_lick_z;
    OE20_testing_AAAA_lick_z;
    OE21_testing_AAAA_lick_z;
    OE24_testing_AAAA_lick_z;
    OE32_testing_AAAA_lick_z;
    OE33_testing_AAAA_lick_z;
    OE34_testing_AAAA_lick_z;
    OE35_testing_AAAA_lick_z;
    OE39_testing_AAAA_lick_z;
    OE40_testing_AAAA_lick_z;
    OE45_testing_AAAA_lick_z;
    OE46_testing_AAAA_lick_z;
    OE47_testing_AAAA_lick_z
];

WT_DCBA_all_z = [
    OE18_testing_DCBA_lick_z;
    OE20_testing_DCBA_lick_z;
    OE21_testing_DCBA_lick_z;
    OE24_testing_DCBA_lick_z;
    OE32_testing_DCBA_lick_z;
    OE33_testing_DCBA_lick_z;
    OE34_testing_DCBA_lick_z;
    OE35_testing_DCBA_lick_z;
    OE39_testing_DCBA_lick_z;
    OE40_testing_DCBA_lick_z;
    OE45_testing_DCBA_lick_z;
    OE46_testing_DCBA_lick_z;
    OE47_testing_DCBA_lick_z
];


shadetype = 'SEM';
sm_val = 20;
Time = -1:0.01:4; % Ensure this matches your actual data resolution

figure
hold on

% Shaded stimulus backgrounds (excluded from legend)
y_limits = [0 0.5];
hA = fill([0 0.5 0.5 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor','none','FaceAlpha',0.3); set(hA,'HandleVisibility','off');
hB = fill([0.5 1 1 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor','none','FaceAlpha',0.3); set(hB,'HandleVisibility','off');
hC = fill([1 1.5 1.5 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor','none','FaceAlpha',0.3); set(hC,'HandleVisibility','off');
hD = fill([1.5 2.0 2.0 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor','none','FaceAlpha',0.3); set(hD,'HandleVisibility','off');
hRW= fill([2 3 3 2],         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.3); set(hRW,'HandleVisibility','off');

% Plot smoothed traces with SEM
stdshade_smooth3(WT_ABCD_all_z, Time, 0.5, 'r',       sm_val, shadetype);     % <--- changed
stdshade_smooth3(WT_AAAA_all_z, Time, 0.5, [0 0.5 0], sm_val, shadetype);     % <--- changed
stdshade_smooth3(WT_DCBA_all_z, Time, 0.5, 'b',       sm_val, shadetype);     % <--- changed

% Labels for each stimulus
text(0.25, 0.475, 'A',  'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.75, 0.475, 'B',  'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.25, 0.475, 'C',  'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.75, 0.475, 'D',  'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(2.5,  0.475, 'RW', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

xlim([-1 4])
ylim(y_limits);
xlabel('Time (s)')
ylabel('Lick Probability')
title('Wild Type - Licking Across Sessions: Test vs Reversal vs Oddball')
box off
grid on

% Legend (force line colors)
hABCD = plot(nan,nan,'r-','LineWidth',1.5);
hAAAA = plot(nan,nan,'-','Color',[0 0.5 0],'LineWidth',1.5);
hDCBA = plot(nan,nan,'b-','LineWidth',1.5);
legend([hABCD hAAAA hDCBA], {'ABCD','AAAA','DCBA'}, 'Location','northeast');


%% GLM Trial-by-Trial on testing day - FMR1KO Mice
close all
clearvars
clc

load('OE43_ABCD_trialbytrial_licking_z.mat');   
load('OE43_DCBA_trialbytrial_licking_z.mat');   
load('OE43_AAAA_trialbytrial_licking_z.mat');   

load('OE44_ABCD_trialbytrial_licking_z.mat');   
load('OE44_DCBA_trialbytrial_licking_z.mat');   
load('OE44_AAAA_trialbytrial_licking_z.mat');   

load('OE48_ABCD_trialbytrial_licking_z.mat');   
load('OE48_DCBA_trialbytrial_licking_z.mat');   
load('OE48_AAAA_trialbytrial_licking_z.mat');   

load('OE49_ABCD_trialbytrial_licking_z.mat');  
load('OE49_DCBA_trialbytrial_licking_z.mat');   
load('OE49_AAAA_trialbytrial_licking_z.mat');  

load('OE50_ABCD_trialbytrial_licking_z.mat');   
load('OE50_DCBA_trialbytrial_licking_z.mat');   
load('OE50_AAAA_trialbytrial_licking_z.mat');   

% --- Combine (following your original selection: OE48–OE50) ---
FMR1KO_ABCD_only_z_data = [
    OE43_ABCD_trialbytrial_licking_z;
    OE44_ABCD_trialbytrial_licking_z;
    OE48_ABCD_trialbytrial_licking_z;
    OE49_ABCD_trialbytrial_licking_z;
    OE50_ABCD_trialbytrial_licking_z
];

FMR1KO_DCBA_only_z_data = [
    OE43_DCBA_trialbytrial_licking_z;
    OE44_DCBA_trialbytrial_licking_z;
    OE48_DCBA_trialbytrial_licking_z;
    OE49_DCBA_trialbytrial_licking_z;
    OE50_DCBA_trialbytrial_licking_z
];

FMR1KO_AAAA_only_z_data = [
    OE43_AAAA_trialbytrial_licking_z;
    OE44_AAAA_trialbytrial_licking_z;
    OE48_AAAA_trialbytrial_licking_z;
    OE49_AAAA_trialbytrial_licking_z;
    OE50_AAAA_trialbytrial_licking_z
];

res  = 0.01;                    
bins = -1:res:4;               

BL_win = [-0.1  0.0];
A_win  = [ 0.4  0.5];
B_win  = [ 0.9  1.0];
C_win  = [ 1.4  1.5];
D_win  = [ 1.9  2.0];

BLi = bins>=BL_win(1) & bins< BL_win(2);
Ai  = bins>=A_win(1)  & bins< A_win(2);
Bi  = bins>=B_win(1)  & bins< B_win(2);
Ci  = bins>=C_win(1)  & bins< C_win(2);
Di  = bins>=D_win(1)  & bins< D_win(2);

seqs     = {'ABCD','DCBA','AAAA'};
datasets = {FMR1KO_ABCD_only_z_data, FMR1KO_DCBA_only_z_data, FMR1KO_AAAA_only_z_data};
winIdx   = {Ai, Bi, Ci, Di;   ... % ABCD: A,B,C,D
            Di, Ci, Bi, Ai;   ... % DCBA: D,C,B,A
            Ai, Ai, Ai, Ai};      % AAAA: A,A,A,A
elems    = {['A';'B';'C';'D'], ...
            ['D';'C';'B';'A'], ...
            ['A';'A';'A';'A']};

% ===== Build trial-by-trial table (response = mean z per window) =====
nPerSet = cellfun(@(X) size(X,1), datasets);
rows    = sum(nPerSet) * 4;

Mouse   = strings(rows,1);
Session = strings(rows,1);
TrialNum= zeros(rows,1);
Seq     = strings(rows,1);
Elem    = strings(rows,1);
Pos     = zeros(rows,1);

Count   = zeros(rows,1);   % mean z-score in the window

row = 0;
for s = 1:numel(seqs)
    X = datasets{s};             
    nT = size(X,1);
    for t = 1:nT
        tr = X(t,:);  % one trial (z-scored time series)

        % four stimulus windows per trial
        for pos = 1:4
            row = row + 1;

            idx = winIdx{s,pos};

            Mouse(row)   = sprintf("OE20_T%04d", t);
            Session(row) = "S1";
            TrialNum(row)= t;
            Seq(row)     = seqs{s};
            Elem(row)    = elems{s}(pos);
            Pos(row)     = pos;

            % mean z-score in the window (continuous response)
            Count(row) = mean(tr(idx));
        end
    end
end

% ===== Assemble table =====
T = table( categorical(Mouse), categorical(Session), TrialNum, ...
           categorical(Seq), categorical(Elem), Pos, Count, ...
           'VariableNames', {'Mouse','Session','TrialNum','Seq','Elem','Pos','Count'});

% Derived covariate: optional within-sequence drift (keep or drop)
T.TrialC = zscore(double(T.TrialNum));

% Factor coding
T.ElemCat = categorical(T.Elem);
T.ElemCat = reordercats(T.ElemCat, {'A','B','C','D'}); 

T.PosCat  = categorical(string(T.Pos));
T.PosCat  = reordercats(T.PosCat, {'1','2','3','4'});    

% ===== Fit Gaussian GLMs on z-scores (no baseline covariate) =====
mdl       = fitglm(T,'Count ~ 1 + ElemCat + PosCat + TrialC', 'Distribution','normal','Link','identity');


mdl_noPos = fitglm(T,'Count ~ 1 + ElemCat + TrialC',          'Distribution','normal','Link','identity');
mdl_noElem= fitglm(T,'Count ~ 1 + PosCat  + TrialC',          'Distribution','normal','Link','identity');

% ===== LR-style comparisons =====
dDev_pos  = max(mdl_noPos.Deviance  - mdl.Deviance,  0);
dDev_elem = max(mdl_noElem.Deviance - mdl.Deviance,  0);

df_pos  = max(height(mdl.Coefficients) - height(mdl_noPos.Coefficients),1);
df_elem = max(height(mdl.Coefficients) - height(mdl_noElem.Coefficients),1);

p_pos   = 1 - chi2cdf(dDev_pos,  df_pos);
p_elem  = 1 - chi2cdf(dDev_elem, df_elem);

n       = height(T);
R2_pos  = 1 - exp(-dDev_pos  / n);
R2_elem = 1 - exp(-dDev_elem / n);
w_pos   = R2_pos / max(R2_pos + R2_elem, eps);
w_elem  = 1 - w_pos;

fprintf('p_pos = %.3g (df=%d), p_elem = %.3g (df=%d) | R2_pos = %.3f, R2_elem = %.3f | w_pos = %.2f, w_elem = %.2f\n', ...
        p_pos, df_pos, p_elem, df_elem, R2_pos, R2_elem, w_pos, w_elem);

disp(mdl.Coefficients)


%% licking plotting


load('OE43_testing_ABCD_lick_z.mat');
load('OE43_testing_AAAA_lick_z.mat');
load('OE43_testing_DCBA_lick_z.mat');

load('OE44_testing_ABCD_lick_z.mat');
load('OE44_testing_AAAA_lick_z.mat');
load('OE44_testing_DCBA_lick_z.mat');

load('OE48_testing_ABCD_lick_z.mat');
load('OE48_testing_AAAA_lick_z.mat');
load('OE48_testing_DCBA_lick_z.mat');

load('OE49_testing_ABCD_lick_z.mat');
load('OE49_testing_AAAA_lick_z.mat');
load('OE49_testing_DCBA_lick_z.mat');

load('OE50_testing_ABCD_lick_z.mat');
load('OE50_testing_AAAA_lick_z.mat');
load('OE50_testing_DCBA_lick_z.mat');

% --- Combine all KO mice for each sequence (z-scored testing averages) ---
FMR1KO_ABCD_all_z = [
    OE43_testing_ABCD_lick_z;
    OE44_testing_ABCD_lick_z;
    OE48_testing_ABCD_lick_z;
    OE49_testing_ABCD_lick_z;
    OE50_testing_ABCD_lick_z
];

FMR1KO_AAAA_all_z = [
    OE43_testing_AAAA_lick_z;
    OE44_testing_AAAA_lick_z;
    OE48_testing_AAAA_lick_z;
    OE49_testing_AAAA_lick_z;
    OE50_testing_AAAA_lick_z
];

FMR1KO_DCBA_all_z = [
    OE43_testing_DCBA_lick_z;
    OE44_testing_DCBA_lick_z;
    OE48_testing_DCBA_lick_z;
    OE49_testing_DCBA_lick_z;
    OE50_testing_DCBA_lick_z
];

shadetype = 'SEM';
sm_val = 20;
Time = -1:0.01:4; % Ensure this matches your actual data resolution

figure
hold on

% Shaded stimulus backgrounds (excluded from legend)
y_limits = [0 0.5];
hA = fill([0 0.5 0.5 0],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.6 0.8 1],     'EdgeColor','none','FaceAlpha',0.3); set(hA,'HandleVisibility','off');
hB = fill([0.5 1 1 0.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.6 0.6],     'EdgeColor','none','FaceAlpha',0.3); set(hB,'HandleVisibility','off');
hC = fill([1 1.5 1.5 1],     [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 1 0.5],       'EdgeColor','none','FaceAlpha',0.3); set(hC,'HandleVisibility','off');
hD = fill([1.5 2.0 2.0 1.5], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [1 0.75 0.95],   'EdgeColor','none','FaceAlpha',0.3); set(hD,'HandleVisibility','off');
hRW= fill([2 3 3 2],         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.3); set(hRW,'HandleVisibility','off');

% Plot smoothed traces with SEM
stdshade_smooth3(FMR1KO_ABCD_all_z, Time, 0.5, 'r',       sm_val, shadetype);     % <--- changed
stdshade_smooth3(FMR1KO_AAAA_all_z, Time, 0.5, [0 0.5 0], sm_val, shadetype);     % <--- changed
stdshade_smooth3(FMR1KO_DCBA_all_z, Time, 0.5, 'b',       sm_val, shadetype);     % <--- changed

% Labels for each stimulus
text(0.25, 0.475, 'A',  'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(0.75, 0.475, 'B',  'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.25, 0.475, 'C',  'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(1.75, 0.475, 'D',  'HorizontalAlignment', 'center', 'FontWeight', 'bold')
text(2.5,  0.475, 'RW', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')

xlim([-1 4])
ylim(y_limits);
xlabel('Time (s)')
ylabel('Lick Probability')
title('FMR1KO Mice - Licking Across Sessions: Test vs Reversal vs Oddball')
box off
grid on

% Legend (force line colors)
hABCD = plot(nan,nan,'r-','LineWidth',1.5);
hAAAA = plot(nan,nan,'-','Color',[0 0.5 0],'LineWidth',1.5);
hDCBA = plot(nan,nan,'b-','LineWidth',1.5);
legend([hABCD hAAAA hDCBA], {'ABCD','AAAA','DCBA'}, 'Location','northeast');


%% t-test
clear; clc;

load('WT_w_pos.mat');
WT_w_pos = w_pos;
      
load('FMR1KO_w_pos.mat');  
FMR1KO_w_pos = w_pos;

WT = WT_w_pos(:);
KO = FMR1KO_w_pos(:);

[h, p, ci, stats] = ttest2(WT_w_pos, FMR1KO_w_pos, 'Vartype', 'unequal');

disp(['Hypothesis test result (h): ', num2str(h)]);
disp(['p-value (p): ', num2str(p, '%.6f')]);      
disp(['Confidence interval (ci): ', num2str(ci')]);
disp(['Test statistic (t): ', num2str(stats.tstat)]);
disp(['Degrees of freedom (df): ', num2str(stats.df)]);
