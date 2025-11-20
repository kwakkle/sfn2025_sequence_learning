%% GLM Combine all single trials 


% Wild Type combine

clear all 
clearvars
close all

files = dir('*_WT_GLM_summary.mat'); 
all = [];

for f = 1:length(files)
    S = load(files(f).name);
    all = [all; S.WT_GLM_summary];
end

R2_elem_WT = [all.R2_elem];
R2_pos_WT  = [all.R2_pos];
w_elem_WT  = [all.w_elem];
w_pos_WT   = [all.w_pos];

% Group means and SEM
mean_R2_elem = mean(R2_elem_WT);
sem_R2_elem  = std(R2_elem_WT) / sqrt(numel(R2_elem_WT));

mean_R2_pos  = mean(R2_pos_WT);
sem_R2_pos   = std(R2_pos_WT)  / sqrt(numel(R2_pos_WT));

mean_w_elem  = mean(w_elem_WT);
sem_w_elem   = std(w_elem_WT)  / sqrt(numel(w_elem_WT));

mean_w_pos   = mean(w_pos_WT);
sem_w_pos    = std(w_pos_WT)   / sqrt(numel(w_pos_WT));

fprintf('R2_pos  = %.3f ± %.3f\n',  mean_R2_pos,  sem_R2_pos);
fprintf('R2_elem = %.3f ± %.3f\n',  mean_R2_elem, sem_R2_elem);
fprintf('w_pos   = %.2f ± %.2f\n',  mean_w_pos,   sem_w_pos);
fprintf('w_elem  = %.2f ± %.2f\n',  mean_w_elem,  sem_w_elem);

% save('WT_R2_pos.mat','R2_pos_WT');
% save('WT_R2_elem.mat','R2_elem_WT');

%% fmr1 combine

clear all 
clearvars
close all

files = dir('*_FMR1KO_GLM_summary.mat');
all = [];

for f = 1:length(files)
    S = load(files(f).name);
    all = [all; S.FMR1KO_GLM_summary];
end

R2_elem_KO = [all.R2_elem];
R2_pos_KO  = [all.R2_pos];
w_elem_KO  = [all.w_elem];
w_pos_KO   = [all.w_pos];

% Group means and SEM
mean_R2_elem = mean(R2_elem_KO);
sem_R2_elem  = std(R2_elem_KO) / sqrt(numel(R2_elem_KO));

mean_R2_pos  = mean(R2_pos_KO);
sem_R2_pos   = std(R2_pos_KO)  / sqrt(numel(R2_pos_KO));

mean_w_elem  = mean(w_elem_KO);
sem_w_elem   = std(w_elem_KO)  / sqrt(numel(w_elem_KO));

mean_w_pos   = mean(w_pos_KO);
sem_w_pos    = std(w_pos_KO)   / sqrt(numel(w_pos_KO));

fprintf('R2_pos  = %.3f ± %.3f\n',  mean_R2_pos,  sem_R2_pos);
fprintf('R2_elem = %.3f ± %.3f\n',  mean_R2_elem, sem_R2_elem);
fprintf('w_pos   = %.2f ± %.2f\n',  mean_w_pos,   sem_w_pos);
fprintf('w_elem  = %.2f ± %.2f\n',  mean_w_elem,  sem_w_elem);

% save('FMR1KO_R2_pos.mat','R2_pos_KO');
% save('FMR1KO_R2_elem.mat','R2_elem_KO');

%% Two-way anova for R2 variance 

clear all
clc
close all

load('WT_R2_pos.mat');      
load('WT_R2_elem.mat');    
load('FMR1KO_R2_pos.mat');  
load('FMR1KO_R2_elem.mat'); 

R2_values = [R2_pos_WT, R2_elem_WT, R2_pos_KO, R2_elem_KO]';
Group     = [ ...
    repmat("WT",     numel(R2_pos_WT),1); ...
    repmat("WT",     numel(R2_elem_WT),1); ...
    repmat("FMR1KO", numel(R2_pos_KO),1); ...
    repmat("FMR1KO", numel(R2_elem_KO),1) ];

Predictor = [ ...
    repmat("Position", numel(R2_pos_WT),1); ...
    repmat("Element",  numel(R2_elem_WT),1); ...
    repmat("Position", numel(R2_pos_KO),1); ...
    repmat("Element",  numel(R2_elem_KO),1) ];

T = table(R2_values, Group, Predictor);

[p,tbl,stats] = anovan(T.R2_values, {T.Group, T.Predictor}, ...
    'model','interaction', ...
    'varnames', {'Genotype','Predictor'});

fprintf('Genotype main effect p = %.4f\n', p(1));
fprintf('Predictor main effect p = %.4f\n', p(2));
fprintf('Interaction p = %.4f\n', p(3));


%% bar plot for GLM 

clear all
clearvars
close all

% Load saved R² data
load('WT_R2_pos.mat');      % R2_pos_WT
load('WT_R2_elem.mat');     % R2_elem_WT
load('FMR1KO_R2_pos.mat');  % R2_pos_KO
load('FMR1KO_R2_elem.mat'); % R2_elem_KO

% Compute means and SEMs
mElem = [mean(R2_elem_WT), mean(R2_elem_KO)];
sElem = [std(R2_elem_WT)/sqrt(numel(R2_elem_WT)), std(R2_elem_KO)/sqrt(numel(R2_elem_KO))];

mPos  = [mean(R2_pos_WT),  mean(R2_pos_KO)];
sPos  = [std(R2_pos_WT)/sqrt(numel(R2_pos_WT)),   std(R2_pos_KO)/sqrt(numel(R2_pos_KO))];

means = [mElem; mPos];   % Rows = [Elem; Pos], Cols = [WT, KO]
sems  = [sElem; sPos];

figure;
b = bar(means, 'grouped', 'FaceColor', 'flat');
hold on;

b(1).CData = [0 0 0.7];   % dark blue
b(2).CData = [0 0.6 0];   % dark green
    
% Error bars
ngroups = size(means,1);   % 2 (Elem, Pos)
nseries = size(means,2);   % 2 (WT, KO)

for j = 1:nseries
    x = b(j).XEndPoints;
    errorbar(x, means(:,j), sems(:,j), 'k', 'linestyle','none', 'linewidth',1);
end

set(gca, 'XTickLabel', {'Element','Position'});
ylabel('R^2');
legend({'Wild Type','FMR1KO'}, 'Location','northeast');
title('GLM Unique Variance (R^2) by Element vs Position');
box on;
set(gcf, 'Color', 'w');

yl = ylim;   % get y-axis limits
ytop = yl(2)*0.95;   % place near top of graph

text(1.2, 0.08, '$$\mathbf{p^{geno} = 0.0139^*}$$', 'Interpreter','latex', ...
    'FontSize',12, 'Color','k', 'HorizontalAlignment','left');
text(1.2, 0.07, '$$\mathbf{p^{pred} = 0.0598}$$', 'Interpreter','latex', ...
    'FontSize',12, 'Color','k', 'HorizontalAlignment','left');
text(1.2, 0.06, '$$\mathbf{p^{inter} = 0.1472}$$', 'Interpreter','latex', ...
    'FontSize',12, 'Color','k', 'HorizontalAlignment','left');


%% Combined barplots: Element & Position coefficients (WT vs Fmr1KO)
clear all
clearvars
close all

% element coefficients 
files_WT = dir('*_WT_GLM_summary.mat');
allElem_WT = [];

for f = 1:length(files_WT)
    S   = load(files_WT(f).name);
    mdl = S.WT_GLM_summary.Coeffs;

    erows  = contains(mdl.Properties.RowNames, 'ElemCat_');
    enames = mdl.Properties.RowNames(erows);
    evals  = mdl.Estimate(erows);

    % enforce order B, C, D
    [~, ord] = ismember({'ElemCat_B','ElemCat_C','ElemCat_D'}, enames);
    allElem_WT = [allElem_WT; evals(ord)'];   % 1x3
end

files_KO = dir('*_FMR1KO_GLM_summary.mat');
allElem_KO = [];

for f = 1:length(files_KO)
    S   = load(files_KO(f).name);
    mdl = S.FMR1KO_GLM_summary.Coeffs;

    erows  = contains(mdl.Properties.RowNames, 'ElemCat_');
    enames = mdl.Properties.RowNames(erows);
    evals  = mdl.Estimate(erows);

    [~, ord] = ismember({'ElemCat_B','ElemCat_C','ElemCat_D'}, enames);
    allElem_KO = [allElem_KO; evals(ord)'];
end

mWT_elem = mean(allElem_WT, 1);  sWT_elem = std(allElem_WT, [], 1)/sqrt(size(allElem_WT,1));
mKO_elem = mean(allElem_KO, 1);  sKO_elem = std(allElem_KO, [], 1)/sqrt(size(allElem_KO,1));
meansElem = [mWT_elem; mKO_elem];
semsElem  = [sWT_elem; sKO_elem];

% ---------- POSITION COEFFICIENTS (relative to pos 1) ----------
allPos_WT = []; allPos_KO = [];

for f = 1:length(files_WT)
    S   = load(files_WT(f).name);
    mdl = S.WT_GLM_summary.Coeffs;

    prows  = contains(mdl.Properties.RowNames, 'PosCat_');
    pnames = mdl.Properties.RowNames(prows);
    pvals  = mdl.Estimate(prows);

    [~, ord] = ismember({'PosCat_2','PosCat_3','PosCat_4'}, pnames);
    allPos_WT = [allPos_WT; pvals(ord)'];
end

for f = 1:length(files_KO)
    S   = load(files_KO(f).name);
    mdl = S.FMR1KO_GLM_summary.Coeffs;

    prows  = contains(mdl.Properties.RowNames, 'PosCat_');
    pnames = mdl.Properties.RowNames(prows);
    pvals  = mdl.Estimate(prows);

    [~, ord] = ismember({'PosCat_2','PosCat_3','PosCat_4'}, pnames);
    allPos_KO = [allPos_KO; pvals(ord)'];
end

mWT_pos = mean(allPos_WT, 1);  sWT_pos = std(allPos_WT, [], 1)/sqrt(size(allPos_WT,1));
mKO_pos = mean(allPos_KO, 1);  sKO_pos = std(allPos_KO, [], 1)/sqrt(size(allPos_KO,1));
meansPos = [mWT_pos; mKO_pos];
semsPos  = [sWT_pos; sKO_pos];

% ---------- PLOT (two panels) ----------
WTcolor = [0 0 0.7];
KOcolor = [0 0.6 0];

figure;
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% ===== Left: Elements =====
ax1 = nexttile;

% Data in the same shape as bar() expects: groups x series (3x2)
Melem = meansElem.';   % 3x2  (groups = B,C,D ; series = WT, KO)
Selem = semsElem.';    % 3x2

b = bar(Melem, 'grouped', 'FaceColor','flat'); hold on
b(1).CData = WTcolor; 
b(2).CData = KOcolor;

for j = 1:2
    x = b(j).XEndPoints;              % 3 x-positions for this series
    y = Melem(:,j);                   % 3 means
    e = Selem(:,j);                   % 3 SEMs
    errorbar(x, y, e, 'k','linestyle','none','linewidth',1);
end
set(gca,'XTickLabel',{'Elem B','Elem C','Elem D'});
ylabel('Coefficient (\beta)');
legend({'WT','Fmr1KO'},'Location','northeast'); legend boxoff
title('Element Coefficients (vs A)');
box on

% ===== Right: Positions =====
ax2 = nexttile;

Mpos = meansPos.';     % 3x2  (groups = 2,3,4 ; series = WT, KO)
Spos = semsPos.';      % 3x2

b = bar(Mpos, 'grouped', 'FaceColor','flat'); hold on
b(1).CData = WTcolor; 
b(2).CData = KOcolor;

for j = 1:2
    x = b(j).XEndPoints;
    y = Mpos(:,j);
    e = Spos(:,j);
    errorbar(x, y, e, 'k','linestyle','none','linewidth',1);
end
set(gca,'XTickLabel',{'Pos 2','Pos 3','Pos 4'});
title('Position Coefficients (vs 1)');
box on

% --- Add p-values for ELEMENTS (left subplot)
axes(ax1);
text(2.1, .8, '$$\mathbf{p^{geno} = 0.0017^{**}}$$', 'Interpreter','latex', ...
    'FontSize',10, 'Color','k', 'HorizontalAlignment','left');
text(2.1, .7, '$$\mathbf{p^{elem} = 0.0918}$$', 'Interpreter','latex', ...
    'FontSize',10, 'Color','k', 'HorizontalAlignment','left');
text(2.1, .6, '$$\mathbf{p^{inter} = 0.3119}$$', 'Interpreter','latex', ...
    'FontSize',10, 'Color','k', 'HorizontalAlignment','left');

% --- Add p-values for POSITIONS (right subplot)
axes(ax2);
text(2.2, -.3, '$$\mathbf{p^{geno} = 0.1308}$$', 'Interpreter','latex', ...
    'FontSize',10, 'Color','k', 'HorizontalAlignment','left');
text(2.2, -.4, '$$\mathbf{p^{pos} = 0.0254^{*}}$$', 'Interpreter','latex', ...
    'FontSize',10, 'Color','k', 'HorizontalAlignment','left');
text(2.2, -.5, '$$\mathbf{p^{inter} = 0.1106}$$', 'Interpreter','latex', ...
    'FontSize',10, 'Color','k', 'HorizontalAlignment','left');
%% Save coefficient for element data 
 
clear; clc; close all

files_WT = dir('*_WT_GLM_summary.mat');
files_KO = dir('*_FMR1KO_GLM_summary.mat');

want = {'ElemCat_B','ElemCat_C','ElemCat_D'};

% WT
ElemCoef_WT = [];   % rows = mice, cols = [B C D]
for f = 1:numel(files_WT)
    S = load(files_WT(f).name);
    C = S.WT_GLM_summary.Coeffs;  % table with RowNames
    names = string(C.Properties.RowNames);
    row = nan(1,3);
    for k = 1:3
        idx = find(names == want{k}, 1);
        if ~isempty(idx), row(k) = C.Estimate(idx); end
    end
    if any(isnan(row))
        fprintf('SKIP WT %s (missing %s)\n', files_WT(f).name, strjoin(want(isnan(row)), ', '));
    else
        ElemCoef_WT(end+1,:) = row; 
    end
end

% KO

ElemCoef_KO = [];
for f = 1:numel(files_KO)
    S = load(files_KO(f).name);
    C = S.FMR1KO_GLM_summary.Coeffs;
    names = string(C.Properties.RowNames);
    row = nan(1,3);
    for k = 1:3
        idx = find(names == want{k}, 1);
        if ~isempty(idx), row(k) = C.Estimate(idx); end
    end
    if any(isnan(row))
        fprintf('SKIP KO %s (missing %s)\n', files_KO(f).name, strjoin(want(isnan(row)), ', '));
    else
        ElemCoef_KO(end+1,:) = row; 
    end
end

Elem_B_WT = ElemCoef_WT(:,1); Elem_C_WT = ElemCoef_WT(:,2); Elem_D_WT = ElemCoef_WT(:,3);
Elem_B_KO = ElemCoef_KO(:,1); Elem_C_KO = ElemCoef_KO(:,2); Elem_D_KO = ElemCoef_KO(:,3);

% save('WT_ElementCoefs.mat',  'ElemCoef_WT','Elem_B_WT','Elem_C_WT','Elem_D_WT');
% save('FMR1KO_ElementCoefs.mat','ElemCoef_KO','Elem_B_KO','Elem_C_KO','Elem_D_KO');

%% Two-way ANOVA on element coefficients (Genotype x Element)

clear; clc; close all

load('WT_ElementCoefs.mat');      % Elem_B_WT, Elem_C_WT, Elem_D_WT
load('FMR1KO_ElementCoefs.mat');  % Elem_B_KO, Elem_C_KO, Elem_D_KO

Elem_values = [ ...
    Elem_B_WT; Elem_C_WT; Elem_D_WT; ...
    Elem_B_KO; Elem_C_KO; Elem_D_KO ];

Genotype = [ ...
    repmat("WT",     numel(Elem_B_WT),1); ...
    repmat("WT",     numel(Elem_C_WT),1); ...
    repmat("WT",     numel(Elem_D_WT),1); ...
    repmat("FMR1KO", numel(Elem_B_KO),1); ...
    repmat("FMR1KO", numel(Elem_C_KO),1); ...
    repmat("FMR1KO", numel(Elem_D_KO),1) ];

Element = [ ...
    repmat("B", numel(Elem_B_WT),1); ...
    repmat("C", numel(Elem_C_WT),1); ...
    repmat("D", numel(Elem_D_WT),1); ...
    repmat("B", numel(Elem_B_KO),1); ...
    repmat("C", numel(Elem_C_KO),1); ...
    repmat("D", numel(Elem_D_KO),1) ];

[p,tbl,stats] = anovan(Elem_values, {categorical(Genotype), categorical(Element)}, ...
    'model','interaction', 'varnames', {'Genotype','Element'});

fprintf('Genotype main effect p = %.4f\n', p(1));
fprintf('Element  main effect p = %.4f\n', p(2));
fprintf('Interaction p          = %.4f\n', p(3));

%% Save coefficient for positions

clear; clc; close all

files_WT = dir('*_WT_GLM_summary.mat');
files_KO = dir('*_FMR1KO_GLM_summary.mat');

want = {'PosCat_2','PosCat_3','PosCat_4'};   

PosCoef_WT = [];   % rows = mice, cols = [2 3 4]
for f = 1:numel(files_WT)
    S = load(files_WT(f).name);
    C = S.WT_GLM_summary.Coeffs;           % table with RowNames
    names = string(C.Properties.RowNames);
    row = nan(1,3);
    for k = 1:3
        idx = find(names == want{k}, 1);
        if ~isempty(idx), row(k) = C.Estimate(idx); end
    end
    if any(isnan(row))
        fprintf('SKIP WT %s (missing %s)\n', files_WT(f).name, strjoin(want(isnan(row)), ', '));
    else
        PosCoef_WT(end+1,:) = row; 
    end
end

PosCoef_KO = [];
for f = 1:numel(files_KO)
    S = load(files_KO(f).name);
    C = S.FMR1KO_GLM_summary.Coeffs;
    names = string(C.Properties.RowNames);
    row = nan(1,3);
    for k = 1:3
        idx = find(names == want{k}, 1);
        if ~isempty(idx), row(k) = C.Estimate(idx); end
    end
    if any(isnan(row))
        fprintf('SKIP KO %s (missing %s)\n', files_KO(f).name, strjoin(want(isnan(row)), ', '));
    else
        PosCoef_KO(end+1,:) = row;
    end
end

Pos_2_WT = PosCoef_WT(:,1); Pos_3_WT = PosCoef_WT(:,2); Pos_4_WT = PosCoef_WT(:,3);
Pos_2_KO = PosCoef_KO(:,1); Pos_3_KO = PosCoef_KO(:,2); Pos_4_KO = PosCoef_KO(:,3);

% save('WT_PositionCoefs.mat',   'PosCoef_WT','Pos_2_WT','Pos_3_WT','Pos_4_WT');
% save('FMR1KO_PositionCoefs.mat','PosCoef_KO','Pos_2_KO','Pos_3_KO','Pos_4_KO');

%% Two way anova on position coefficients 

clear; clc; close all

load('WT_PositionCoefs.mat');      % Pos_2_WT, Pos_3_WT, Pos_4_WT
load('FMR1KO_PositionCoefs.mat');  % Pos_2_KO, Pos_3_KO, Pos_4_KO

Pos_values = [ ...
    Pos_2_WT; Pos_3_WT; Pos_4_WT; ...
    Pos_2_KO; Pos_3_KO; Pos_4_KO ];

Genotype = [ ...
    repmat("WT",     numel(Pos_2_WT),1); ...
    repmat("WT",     numel(Pos_3_WT),1); ...
    repmat("WT",     numel(Pos_4_WT),1); ...
    repmat("FMR1KO", numel(Pos_2_KO),1); ...
    repmat("FMR1KO", numel(Pos_3_KO),1); ...
    repmat("FMR1KO", numel(Pos_4_KO),1) ];

Position = [ ...
    repmat("2", numel(Pos_2_WT),1); ...
    repmat("3", numel(Pos_3_WT),1); ...
    repmat("4", numel(Pos_4_WT),1); ...
    repmat("2", numel(Pos_2_KO),1); ...
    repmat("3", numel(Pos_3_KO),1); ...
    repmat("4", numel(Pos_4_KO),1) ];

[p,tbl,stats] = anovan(Pos_values, {categorical(Genotype), categorical(Position)}, ...
    'model','interaction', 'varnames', {'Genotype','Position'});

fprintf('Genotype main effect p = %.4f\n', p(1));
fprintf('Position main effect p = %.4f\n', p(2));
fprintf('Interaction p          = %.4f\n', p(3));
