%% GLM Trial-by-Trial on testing day - Wild Type Mice
close all
clearvars
clc

mouseid = 'OE47';  

seqs = {'ABCD','DCBA','AAAA'};
for i = 1:numel(seqs)
    seq   = seqs{i};
    fname = sprintf('%s_%s_trialbytrial_licking.mat', mouseid, seq);
    L     = load(fname);
    vname = sprintf('%s_%s_trialbytrial_licking', mouseid, seq);
    switch seq
        case 'ABCD'
            ABCD_only_data = L.(vname);
        case 'DCBA'
            DCBA_only_data = L.(vname);
        case 'AAAA'
            AAAA_only_data = L.(vname);
    end
end

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
datasets = {ABCD_only_data, DCBA_only_data, AAAA_only_data};

% Positions are fixed windows across all sequences
winIdx   = {Ai, Bi, Ci, Di;   ... % ABCD positions 1..4
            Ai, Bi, Ci, Di;   ... % DCBA positions 1..4
            Ai, Bi, Ci, Di};      % AAAA positions 1..4
% Element identity per position (sequence-specific)
elems    = {['A';'B';'C';'D'], ...
            ['D';'C';'B';'A'], ...
            ['A';'A';'A';'A']};

nPerSet = cellfun(@(X) size(X,1), datasets);       % #trials in each sequence
rows    = sum(nPerSet) * 4;                        % 4 windows (A–D) per trial (no +1 for BL)

Mouse   = strings(rows,1);
Session = strings(rows,1);
TrialNum= zeros(rows,1);
Seq     = strings(rows,1);
Elem    = strings(rows,1);
Pos     = zeros(rows,1);    % positions 1..4

Count   = zeros(rows,1);
WindowSec = zeros(rows,1);
PreTrialCount      = zeros(rows,1);
PreTrialWindowSec  = zeros(rows,1);

row = 0;
for s = 1:numel(seqs)
    X = datasets{s};                 % trials x time
    nT = size(X,1);
    for t = 1:nT
        tr = X(t,:);                 % 1 x 501 trial trace (probability per bin)

        % per-trial baseline (used as covariate, but NOT added as a row)
        bl_cnt  = sum(tr(BLi));
        nBinsBL = sum(BLi);

        % --- A/B/C/D rows (sequence-specific windows) PER TRIAL ---
        for pos = 1:4
            row = row + 1;

            idx      = winIdx{s,pos};
            nBinsWin = sum(idx);

            Mouse(row) = sprintf("T%04d", t);
            Session(row) = "S1";
            TrialNum(row)= t;
            Seq(row)     = seqs{s};
            Elem(row)    = elems{s}(pos);
            Pos(row)     = pos;

            % sum of per-bin probabilities ~ expected lick count in window
            Count(row)     = sum(tr(idx));
            WindowSec(row) = nBinsWin * res;

            % carry the per-trial baseline into covariates
            PreTrialCount(row)     = bl_cnt;
            PreTrialWindowSec(row) = nBinsBL * res;
        end
    end
end

% Assemble Table
T = table( categorical(Mouse), categorical(Session), TrialNum, ...
           categorical(Seq), categorical(Elem), Pos, ...
           Count, WindowSec, PreTrialCount, PreTrialWindowSec, ...
    'VariableNames', {'Mouse','Session','TrialNum','Seq','Elem','Pos', ...
                      'Count','WindowSec','PreTrialCount','PreTrialWindowSec'});

% Derived Covariates
% T.Offset = log(T.WindowSec);                          % exposure = seconds per window
BL_rate  = T.PreTrialCount ./ max(T.PreTrialWindowSec, eps);
T.PreReg = zscore(log1p(BL_rate));                    % baseline rate per trial (z-scored)
T.TrialC = zscore(double(T.TrialNum));                % optional within-seq drift

% SINGLE GLM: Elem ref = A, Pos ref = 1 
T.ElemCat = categorical(T.Elem);
T.ElemCat = reordercats(T.ElemCat, {'A','B','C','D'});   % A is reference

T.PosCat  = categorical(string(T.Pos));
T.PosCat  = reordercats(T.PosCat, {'1','2','3','4'});    % 1 is reference

% Full model: splits variance between elements & positions
mdl = fitglm(T, 'Count ~ 1 + ElemCat + PosCat + PreReg + TrialC', ...
             'Distribution','poisson','Link','log','Offset',T.Offset);

% Reduced models (drop one factor at a time) for variance partitioning
mdl_noPos  = fitglm(T, 'Count ~ 1 + ElemCat + PreReg + TrialC', ...
                    'Distribution','poisson','Link','log','Offset',T.Offset);
mdl_noElem = fitglm(T, 'Count ~ 1 + PosCat  + PreReg + TrialC', ...
                    'Distribution','poisson','Link','log','Offset',T.Offset);

% LR stats (unique contributions)
dDev_pos  = max(mdl_noPos.Deviance  - mdl.Deviance,  0);
dDev_elem = max(mdl_noElem.Deviance - mdl.Deviance,  0);

df_pos  = max(height(mdl.Coefficients) - height(mdl_noPos.Coefficients),  1);
df_elem = max(height(mdl.Coefficients) - height(mdl_noElem.Coefficients), 1);

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


WT_GLM_summary = struct( ...
    'MouseID', mouseid, ...
    'R2_elem', R2_elem, ...
    'R2_pos',  R2_pos, ...
    'w_elem',  w_elem, ...
    'w_pos',   w_pos, ...
    'p_elem',  p_elem, ...
    'p_pos',   p_pos, ...
    'Coeffs',  mdl.Coefficients);

% savefile = sprintf('%s_WT_GLM_summary.mat', mouseid);
% save(savefile, 'WT_GLM_summary');

%% GLM Trial-by-Trial on testing day - Fmr1 KO mice
close all
clearvars
clc

mouseid = 'OE50';  

seqs = {'ABCD','DCBA','AAAA'};
for i = 1:numel(seqs)
    seq   = seqs{i};
    fname = sprintf('%s_%s_trialbytrial_licking.mat', mouseid, seq);
    L     = load(fname);
    vname = sprintf('%s_%s_trialbytrial_licking', mouseid, seq);
    switch seq
        case 'ABCD'
            ABCD_only_data = L.(vname);
        case 'DCBA'
            DCBA_only_data = L.(vname);
        case 'AAAA'
            AAAA_only_data = L.(vname);
    end
end

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
datasets = {ABCD_only_data, DCBA_only_data, AAAA_only_data};

% Positions are fixed windows across all sequences
winIdx   = {Ai, Bi, Ci, Di;   ... % ABCD positions 1..4
            Ai, Bi, Ci, Di;   ... % DCBA positions 1..4
            Ai, Bi, Ci, Di};      % AAAA positions 1..4
% Element identity per position (sequence-specific)
elems    = {['A';'B';'C';'D'], ...
            ['D';'C';'B';'A'], ...
            ['A';'A';'A';'A']};

% We DO NOT create a BL row. Only A/B/C/D positions per trial.
nPerSet = cellfun(@(X) size(X,1), datasets);       % #trials in each sequence
rows    = sum(nPerSet) * 4;                        % 4 windows (A–D) per trial (no +1 for BL)

Mouse   = strings(rows,1);
Session = strings(rows,1);
TrialNum= zeros(rows,1);
Seq     = strings(rows,1);
Elem    = strings(rows,1);
Pos     = zeros(rows,1);    % positions 1..4

Count   = zeros(rows,1);
WindowSec = zeros(rows,1);
PreTrialCount      = zeros(rows,1);
PreTrialWindowSec  = zeros(rows,1);

row = 0;
for s = 1:numel(seqs)
    X = datasets{s};                 % trials x time
    nT = size(X,1);
    for t = 1:nT
        tr = X(t,:);                 % 1 x 501 trial trace (probability per bin)

        % per-trial baseline (used as covariate, but NOT added as a row)
        bl_cnt  = sum(tr(BLi));
        nBinsBL = sum(BLi);

        % --- A/B/C/D rows (sequence-specific windows) PER TRIAL ---
        for pos = 1:4
            row = row + 1;

            idx      = winIdx{s,pos};
            nBinsWin = sum(idx);

            Mouse(row) = sprintf("T%04d", t);
            Session(row) = "S1";
            TrialNum(row)= t;
            Seq(row)     = seqs{s};
            Elem(row)    = elems{s}(pos);
            Pos(row)     = pos;

            % sum of per-bin probabilities ~ expected lick count in window
            Count(row)     = sum(tr(idx));
            WindowSec(row) = nBinsWin * res;

            % carry the per-trial baseline into covariates
            PreTrialCount(row)     = bl_cnt;
            PreTrialWindowSec(row) = nBinsBL * res;
        end
    end
end

% Assemble Table
T = table( categorical(Mouse), categorical(Session), TrialNum, ...
           categorical(Seq), categorical(Elem), Pos, ...
           Count, WindowSec, PreTrialCount, PreTrialWindowSec, ...
    'VariableNames', {'Mouse','Session','TrialNum','Seq','Elem','Pos', ...
                      'Count','WindowSec','PreTrialCount','PreTrialWindowSec'});

% ---- Derived covariates (PER TRIAL) ----
T.Offset = log(T.WindowSec);                          % exposure = seconds per window
BL_rate  = T.PreTrialCount ./ max(T.PreTrialWindowSec, eps);
T.PreReg = zscore(log1p(BL_rate));                    % baseline rate per trial (z-scored)
T.TrialC = zscore(double(T.TrialNum));                % optional within-seq drift

% SINGLE GLM: Elem ref = A, Pos ref = 1 
T.ElemCat = categorical(T.Elem);
T.ElemCat = reordercats(T.ElemCat, {'A','B','C','D'});   % A is reference

T.PosCat  = categorical(string(T.Pos));
T.PosCat  = reordercats(T.PosCat, {'1','2','3','4'});    % 1 is reference

% Full model: splits variance between elements & positions
mdl = fitglm(T, 'Count ~ 1 + ElemCat + PosCat + PreReg + TrialC', ...
             'Distribution','poisson','Link','log','Offset',T.Offset);

% Reduced models (drop one factor at a time) for variance partitioning
mdl_noPos  = fitglm(T, 'Count ~ 1 + ElemCat + PreReg + TrialC', ...
                    'Distribution','poisson','Link','log','Offset',T.Offset);
mdl_noElem = fitglm(T, 'Count ~ 1 + PosCat  + PreReg + TrialC', ...
                    'Distribution','poisson','Link','log','Offset',T.Offset);

% LR stats (unique contributions)
dDev_pos  = max(mdl_noPos.Deviance  - mdl.Deviance,  0);
dDev_elem = max(mdl_noElem.Deviance - mdl.Deviance,  0);

df_pos  = max(height(mdl.Coefficients) - height(mdl_noPos.Coefficients),  1);
df_elem = max(height(mdl.Coefficients) - height(mdl_noElem.Coefficients), 1);

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


FMR1KO_GLM_summary = struct( ...
    'MouseID', mouseid, ...
    'R2_elem', R2_elem, ...
    'R2_pos',  R2_pos, ...
    'w_elem',  w_elem, ...
    'w_pos',   w_pos, ...
    'p_elem',  p_elem, ...
    'p_pos',   p_pos, ...
    'Coeffs',  mdl.Coefficients);
% 
% savefile = sprintf('%s_FMR1KO_GLM_summary.mat', mouseid);
% save(savefile, 'FMR1KO_GLM_summary');

